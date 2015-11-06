#include<iostream>
#include <fstream>
#include<stdio.h>
#include<cstdlib>
#include <sstream>
#include <limits.h>
#include<string>
#include<qapio.h>
#include<curand_kernel.h>



using namespace std;

#define width 4950
#define TABULIST 10000
#define MAXTHREADS 1024 

#define AODMIN 0.1
#define AODMAX 0.33


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true) {

  if (code != cudaSuccess) {
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

/*
 * use basic formula to compute cost of a permuation 
 */
__device__ int initCost(int *d_flows, int *d_dist, int *d_sol, int nsize, int tidx) {
  int calcost = 0;
  int index = tidx * nsize;
  for(int i = 0; i < nsize - 1; i++) {
    for(int j = i + 1; j < nsize; j++)
      calcost = calcost 
	+ (d_flows[(d_sol[index + i] - 1) * nsize + (d_sol[index + j] - 1)]) 
	* d_dist[i * nsize  + j];
  }
  for(int k = 1; k < nsize; k++) {
    for(int l = 0; l < k;l++)
      calcost = calcost 
	+ d_flows[(d_sol[index + k] - 1) * nsize + (d_sol[index + l]- 1)] 
	* d_dist[k * nsize + l];	
  }
  return calcost;
}

/*
 * compute cost of a permutation based on Burkard 
 */
__device__ int neighborCost(int *d_flows, int *d_dist, int *d_sol, int nsize, int tidx, int i, int j) {

  int offset = tidx * nsize;

  int iUnit = d_sol[offset + i];
  int jUnit = d_sol[offset + j];

  int ccost = 0, gcost = 0, hcost = 0;

  for(int k = 0; k < nsize; k++) {
    int kUnit = d_sol[offset + k];
    if (k != i && k != j)  {
      gcost = (d_dist[j * nsize + k] - d_dist[i * nsize + k]) *
	(d_flows[(iUnit - 1) * nsize + (kUnit - 1)] - d_flows[(jUnit-1) * nsize + (kUnit - 1)]); 
      hcost = (d_dist[k * nsize + j] - d_dist[k * nsize + i]) * 
	(d_flows[(kUnit - 1) * nsize + (iUnit - 1)] - d_flows[(kUnit - 1) * nsize + (jUnit - 1)]);
      ccost = ccost + (gcost + hcost);
    }
  }
  return ccost;
}

__device__ void copyFromFlattened(int *dest, int *src, int size, int tidx) {
  int offset = tidx * size;
  for(int i = 0; i < size; i++)
    dest[i] = src[offset + i];
}

/* 
 *   copy src permutation to dest permutation 
 */
__device__ void copy(int *dest, int *src, int size, int tidx) {
  int offset = tidx * size;
  for(int i = 0; i < size; i++)
    dest[offset + i] = src[offset + i];
  return;
}

void copyToFlattened(int *dest, int *src, int size, int tidx) {
  int offset = tidx * size;
  for(int i = 0; i < size; i++)
    dest[offset + i] = src[i];
}


/* 
 * swap units 
 */
__device__ void swap(int *a,int *b) {
  int temp=0;
  temp = *a;
  *a = *b;
  *b = temp;
}

__device__ int findmin(int *a, int start, int length) {

  int min = a[start];
  int minIndex = start;
  for(int i = start; i < length; i++) {
    if (a[i] < min) {
      min = a[i];
      minIndex = i;
    }
  }
  swap(&a[start],&a[minIndex]);
  return min;
}


__device__ void sort(int *arr,int length) {
  int key, i;
  for(int j = 1; j < length; j++) {
     key = arr[j];
     i = j-1;
     while(arr[i] > key && i >= 0) {
       arr[i+1]=arr[i];
       i--;
     }
     arr[i+1]=key;
  }
}


__device__ int findMinPos(int least, int *d_divresult, int r, int tidx) {
  for(int i = 0; i < r; i++)
    if (least == d_divresult[tidx * r + i])
      return i;
  return 0;
}

__global__ void child_kernel(int tidx, int nsize, int *d_dist, int *d_flows, int *d_sol, int baseCost, 
			     int *d_pos, int *d_divresult, int thrds) {

  int tx = threadIdx.x + blockDim.x * blockIdx.x;

  int ipos =  d_pos[tx * 2];
  int jpos =  d_pos[tx * 2 + 1];
  
  int offset = tidx * nsize; 
  int cost; 
  
  cost = (d_dist[jpos * nsize + ipos] - d_dist[ipos * nsize + jpos])
    * (d_flows[(d_sol[offset + ipos] - 1) * nsize + (d_sol[offset + jpos] - 1)] 
       - d_flows[(d_sol[offset + jpos] - 1) * nsize + (d_sol[offset + ipos] - 1)]);
  
  cost += (d_dist[jpos * nsize + jpos] - d_dist[ipos * nsize + ipos]) 
    * (d_flows[(d_sol[offset + ipos] - 1) * nsize + (d_sol[offset + ipos] - 1)] 
       - d_flows[(d_sol[offset + jpos] - 1) * nsize + (d_sol[offset + jpos] - 1)]);
  
  // delta 
  cost += neighborCost(d_flows, d_dist, d_sol, nsize, tidx, ipos, jpos); 
  cost += baseCost;
  
  d_divresult[tidx * thrds + tx] = cost;
  return;
}


__global__ void tabu(int *d_dist,int *d_flows,int *d_sol,int nsize, 
		     int *d_bestsofar,int *d_bestcostsofar,int *d_pos,int *d_newarray,
		     int *d_divresult, curandState *state, int hoods) {	


  int d_tsol[4950];

  int least=0, count=0,lespos=0;
  int t=0,nbk=0,ntr=0;

  short tabu[TABULIST];

  float AOD = 0;
  int min = nsize * AODMIN;
  int max = nsize * AODMAX;
  int n = (max - min) + 1;

  float farray[110];
  int array[110];

  for(int i = min; i <= max; i++)
      array[i - min] = i;
  float jk = 1;
  for(int i = 0; i < n; i++,jk++)
    if(jk <= n)
      farray[i] = (float) jk/n;

  int tidx = threadIdx.x + blockDim.x * blockIdx.x;

  int baseCost = initCost(d_flows,d_dist,d_sol,nsize,tidx);
  d_bestcostsofar[tidx] = baseCost;
  copy(d_bestsofar, d_sol, nsize, tidx);
  copy(d_newarray, d_sol, nsize, tidx);

  curand_init(17, tidx, 0, &state[tidx]);

  // determine child launch configuration 
  int neighbors =  ((nsize - 1) * nsize) / 2;
  int threadsPerBlock;
  int blocksPerGrid;
  
  if (neighbors < MAXTHREADS) {
    threadsPerBlock = neighbors;
    blocksPerGrid = 1;
  }
  else {
    threadsPerBlock = 32;
    blocksPerGrid = (neighbors + threadsPerBlock - 1) / threadsPerBlock;
  }

  // each child returns one result 
  int results = threadsPerBlock * blocksPerGrid;

  // ALGM parameter
  for(int l = 0 ; l < nsize * hoods; l++) {	
    count = 0;

    child_kernel<<<blocksPerGrid,threadsPerBlock>>>
      (tidx,nsize,d_dist,d_flows,d_sol,baseCost,d_pos,d_divresult, results);	
			
    if (cudaSuccess != cudaGetLastError())
      return;
    if (cudaSuccess != cudaDeviceSynchronize())
      return;
    
    copyFromFlattened(d_tsol, d_divresult, results, tidx);
    
    //    sort(d_tsol,r);
    //    least = d_tsol[count]; 
    least = findmin(d_tsol, count, results);
    lespos = findMinPos(least, d_divresult, results, tidx);

    do {

      int npos = tidx * nsize + d_pos[lespos*2];
      int nposNext = tidx * nsize + d_pos[lespos*2+1];

      if(d_newarray[npos] < d_newarray[nposNext]) {
	ntr = d_newarray[npos];
	nbk = d_newarray[nposNext];
      }
      else {
	ntr = d_newarray[nposNext];
	nbk = d_newarray[npos];
      }
      
      if (tabu[ntr * nsize + nbk] <= l) {
	AOD = curand_uniform(&state[tidx]);
	if(AOD < farray[0])
	  t = array[0];
	else {
	  for(int j = 0;j < n-1; j++) {
	    if(farray[j] < AOD && AOD < farray[j + 1])
	      t = array[j + 1];
	  }
	}
	tabu[ntr * nsize + nbk] = l + t;

	swap(&d_newarray[npos],&d_newarray[nposNext]);
	if(least < d_bestcostsofar[tidx]) {
	  d_bestcostsofar[tidx] = least;
	  for(int j = 0; j< nsize; j++)
	      d_bestsofar[tidx * nsize + j] = d_newarray[tidx*nsize+j];
	  tabu[nbk * nsize + ntr] = tabu[nbk * nsize + ntr] + 1;
	}
	break;
      }
      if (least < d_bestcostsofar[tidx]) {
	swap(&d_newarray[npos],&d_newarray[nposNext]);
      	d_bestcostsofar[tidx] = least;
	for(int j = 0; j < nsize; j++)
	  d_bestsofar[tidx * nsize + j] = d_newarray[tidx* nsize + j];
	break;
      }

      count++;
      // least = d_tsol[count]; 
      least = findmin(d_tsol, count, results);
      lespos = findMinPos(least, d_divresult, results, tidx);

    } while (true);

    for(int j = 0;j < nsize; j++)
      d_sol[tidx * nsize + j] = d_newarray[tidx * nsize + j];
    baseCost = initCost(d_flows, d_dist, d_sol, nsize, tidx);
  }

  free(tabu);
  free(d_tsol);
}


void perm(int *init_sol,int nsize) {
  for(int i = nsize-1; i > 0; i--) {
    int j= rand() % (i+1);
    int temp = init_sol[i];
    init_sol[i] = init_sol[j];
    init_sol[j]=temp;
  }
}

bool unique(int *h_sol,int *init_sol,int nsize,int k) {
  int ntr = 0;
  for(int lk=0;lk<k;lk++) {
    for(int l=0;l<nsize;l++) {
      if(h_sol[lk*nsize+l] != init_sol[l]) {
	ntr++;
	break;
      }
    }
    
  }
  if(ntr == k)
    return true;
  else
    return false;
}



int main(int argc,char *argv[]) {

  curandState* devStates;

  clock_t cpu_startTime, cpu_endTime;
  double cpu_ElapseTime=0;
  cpu_startTime = clock();

  cudaEvent_t start , stop;
  float ctime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  
  if (argc != 5) {
    cout << "usage: " << endl;
    cout << "\t./tabu_dyn datafile solns randseed hoods" << endl;
    exit(1);
  }

  string filename = argv[1];
  int solns = atoi(argv[2]);
  int iseed = atoi(argv[3]);
  int hoods = atoi(argv[4]);


  // read data from file 
  int size;
  int **array;
  readData(filename, &array, &size);
 
  // split and flatten matrix
  int *h_flows, *h_dist;
  splitAndFlattenInt(&h_flows, &h_dist, array, size);
    
  srand(iseed);

  // set up initial solution 
  int *init_sol = (int *) malloc(size * sizeof(int));
  for(int i = 0; i < size; i++)
    init_sol[i] = i + 1;

  // create set of initial solutions via random permutations
  int *h_sol = (int *) malloc( sizeof(int) * (solns) * size); 
  for(int k = 0; k < solns; k++) {
    do {
      perm(init_sol,size);
      copyToFlattened(h_sol,init_sol,size,k);
    } while(!unique(h_sol,init_sol,size,k));
  }


  int *h_pos = (int *) malloc(sizeof(int) * (size*(size-1))/2 * 2);
  int q = 0;
  for(int i = 0;i < size-1; i++) {
    for(int j = i + 1;j < size; j++) {
      h_pos[q * 2] = i;
      h_pos[q * 2 + 1] = j;
      q++;
    }
  }
  
  int *h_bestsofar = (int *)malloc(solns * size * sizeof(int));
  int *h_bestcostsofar = (int *)malloc((solns) * sizeof(int));

  cudaEventRecord(start,0);
  
  int *d_dist = NULL,*d_flows = NULL ,*d_sol = NULL,*d_pos=NULL;
  int *d_bestsofar=NULL,*d_bestcostsofar=NULL,*d_newarray=NULL;
  int *d_divresult;

  gpuErrchk( cudaMalloc(&devStates, solns * sizeof( curandState )) );
  gpuErrchk( cudaMalloc((void **)&d_bestcostsofar, solns * sizeof(int)) );
  gpuErrchk( cudaMalloc((void **)&d_bestsofar, solns * size * sizeof(int)) );
  gpuErrchk( cudaMalloc((void **)&d_divresult,(solns*(size*(size-1))/2) * sizeof(int)) );
  gpuErrchk( cudaMalloc((void **)&d_newarray,(solns) * size * sizeof(int)));
  gpuErrchk( cudaMalloc((void **)&d_dist,size*size*sizeof(int)) );
  gpuErrchk( cudaMalloc((void **)&d_flows,size*size*sizeof(int)) );
  gpuErrchk( cudaMalloc((void **)&d_pos,(size*(size-1))/2 * 2 * sizeof(int)));
  gpuErrchk( cudaMalloc((void **)&d_sol,solns*size*sizeof(int)) );
  
  // copy from host to device 	
  gpuErrchk( cudaMemcpy(d_dist,h_dist,size*size*sizeof(int),cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(d_flows,h_flows,size*size*sizeof(int),cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(d_sol,h_sol,solns*size*sizeof(int),cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(d_pos,h_pos,(size*(size-1))/2 * 2 * sizeof(int),cudaMemcpyHostToDevice) );
  
  int threadsPerBlock=32;
  int blockPerGrid = (solns + threadsPerBlock - 1) / threadsPerBlock;
  
#ifdef DEBUG
  cout<<"input file name:" << filename << endl; 
  cout<<"size of the array:"<<size;						
  cout<<"number of initial solutions:"<<solns<<endl;
  cout<<"number of blocks:"<<blockPerGrid<<" "<<endl;
  cout<<"number of threads:"<<threadsPerBlock<<" "<<endl;
#endif

  tabu<<<blockPerGrid,threadsPerBlock>>>(d_dist,d_flows,d_sol,size,
					 d_bestsofar,d_bestcostsofar,d_pos,
					 d_newarray,d_divresult,devStates,hoods);
  
  gpuErrchk( cudaPeekAtLastError() );
  if (cudaSuccess != cudaGetLastError())
    return 1;
  if (cudaSuccess != cudaDeviceSynchronize())
    return 2;
  
  // copy from device to host
  gpuErrchk( cudaMemcpy(h_bestcostsofar,d_bestcostsofar,solns * sizeof(int),cudaMemcpyDeviceToHost) );
  gpuErrchk( cudaMemcpy(h_bestsofar, d_bestsofar , solns * size* sizeof(int),cudaMemcpyDeviceToHost )); 
  
  cudaEventRecord(stop,0);
  cudaEventSynchronize(start);
  cudaEventSynchronize(stop);
    
  int bestsoln = 0;
  int bestcost = h_bestcostsofar[0];
  for(int i = 1;i < solns; i++) {
    if(h_bestcostsofar[i] < bestcost) {
      bestcost = h_bestcostsofar[i];
      bestsoln = i;
    }
  }

  cout << "problem size: " << size << endl;
  cout << "best cost: " << bestcost << endl;
  cout << "best solution: ";
  for(int j = 0; j < size; j++)
    cout<< h_bestsofar[bestsoln * size + j] << " ";
  cout << endl;

  cpu_endTime = clock();
  cpu_ElapseTime= ((cpu_endTime - cpu_startTime) /(double) CLOCKS_PER_SEC);
 
  fprintf(stdout, "total exec time (s):\t %2.2f\n", cpu_ElapseTime); 
  cudaEventElapsedTime(&ctime, start , stop);
  fprintf(stdout, "kernel exec:\t %2.2f\n", (ctime / 1000));

  cudaEventDestroy(start);
  cudaEventDestroy(stop);


  free(array);
  free(h_dist);
  free(h_flows);
  free(init_sol);
  free(h_sol);
  free(h_pos);
  free(h_bestsofar);
  free(h_bestcostsofar);
  cudaFree(d_dist);
  cudaFree(d_flows);
  cudaFree(d_sol);
  cudaFree(d_bestcostsofar);
  cudaFree(d_bestsofar);
  cudaFree(d_pos);

  return 0;
}
