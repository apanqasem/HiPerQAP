#include<iostream>
#include <fstream>
#include<stdio.h>
#include<cstdlib>
#include <sstream>
#include<string>
#include<qapio.h>
#include<2optlib.h>

using namespace std;

#define TILE_WIDTH 100

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true) {
  if (code != cudaSuccess) {
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
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

/* 
 * swap units 
 */
__device__ void swap(int *a,int *b) {
  int temp=0;
  temp = *a;
  *a = *b;
  *b = temp;
}

/*
 * use basic formula to compute cost of a permuation 
 */
__device__ int initCost(short *d_flows, short *d_dist, int *d_sol, int nsize, int tidx) {
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
__device__ int neighborCost(short *d_flows, short *d_dist, int *d_sol, int nsize, int tidx, int i, int j) {

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

  
__global__ void twoOptShared(short *d_dist,short *d_flows,int *d_sol,int nsize,int row,int *d_result,int *d_bestsofar,int *d_bestcostsofar,int *d_newarray)
{	int totalcost =0,tcost;
  int ipos,jpos;
  int ti = threadIdx.x;
  int b0 = blockDim.x * blockIdx.x;
  int tidx = ti + b0;

  __shared__ short int a[TILE_WIDTH*TILE_WIDTH];
  __shared__ short int b[TILE_WIDTH*TILE_WIDTH];

  if(ti<nsize)
    {	
      for(int j=0;j<nsize;j++)
	{
	  a[ti*nsize+j] = d_dist[ti*nsize+j];
	  b[ti*nsize+j] = d_flows[ti*nsize+j];
	}
			
    }
  __syncthreads();

  if(tidx < row)
    {	
      totalcost=initCost(b,a,d_sol,nsize,tidx);
      d_result[tidx] = totalcost;
      __syncthreads();
      d_bestcostsofar[tidx]=d_result[tidx];
      copy(d_bestsofar,d_sol,nsize,tidx);
    }
  for(int l=0;l<nsize;l++)
    {	if(tidx<row)
	{	for(int k=0;k<nsize;k++)
	    { 				
	      for(int j=k+1;j<nsize;j++)
		{				ipos =0;
		  jpos=0;
		  copy(d_newarray,d_sol,nsize,tidx);
		  swap(&d_newarray[tidx*nsize+k],&d_newarray[tidx*nsize+j]);
		  ipos =  k;
		  jpos =  j;
		  int dcost=0,ecost=0,fcost=0;
		  dcost = (a[jpos*nsize+ipos] - a[ipos*nsize+jpos])*(b[(d_sol[(tidx * nsize)+ipos]-1) * nsize + (d_sol[(tidx * nsize)+jpos] - 1)] - b[(d_sol[(tidx * nsize)+jpos]-1) * nsize + (d_sol[(tidx * nsize)+ipos] - 1)]);
		  ecost = (a[jpos*nsize+jpos] - a[ipos*nsize+ipos])*(b[(d_sol[(tidx * nsize)+ipos]-1) * nsize + (d_sol[(tidx * nsize)+ipos] - 1)] - b[(d_sol[(tidx * nsize)+jpos]-1) * nsize + (d_sol[(tidx * nsize)+jpos] - 1)]);
		  fcost = dcost + ecost;
		  int totcost=0,delta=0;
		  totcost = neighborCost(b,a,d_sol,nsize,tidx,ipos,jpos);
		  delta = fcost + totcost;
		  tcost = d_result[tidx] + delta;
		  if(tcost<d_bestcostsofar[tidx])
		    {	
		      d_bestcostsofar[tidx]=tcost;
		      for(int j=0;j<nsize;j++)
			d_bestsofar[tidx * nsize + j] = d_newarray[tidx * nsize + j];
		    }
		}
	    }
			
		
	  for(int jc=0;jc<nsize;jc++)
	    {
	      d_sol[tidx*nsize+jc]= d_bestsofar[tidx*nsize+jc];
	    }
	  d_result[tidx] = d_bestcostsofar[tidx];
	}
    }
}



int main(int argc,char *argv[]) {

#ifdef PROFILE
  clock_t cpu_startTime, cpu_endTime; 
  double cpu_ElapseTime=0;

  cpu_startTime = clock();
#endif

  if (argc != 4) {
    cout << "usage: " << endl;
    cout << "\t./2opt datafile solns iterations randseed" << endl;
    exit(1);
  }

#ifdef DEBUG 
  cout <<"input file name: "<< argv[1] << endl; 
  cout <<"initial solutions: " << argv[2] << endl;
  cout <<"seed value: " << argv[3] << endl;
#endif 

  string filename = argv[1];
  int solns = atoi(argv[2]);
  int iseed = atoi(argv[3]);

#ifdef PROFILE
  cudaEvent_t start , stop;
  float ctime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
#endif 


  // read data from file 
  int size;
  int **array;
  readData(filename, &array, &size);

  // // Split and flatten matrix
   short *h_flows, *h_dist;
   splitAndFlatten(&h_flows, &h_dist, array, size);

#ifdef DEBUG
  cout << "problem size:" << size << endl;						
  printFlattenedArray(h_dist, size, size, "distances:");
  printFlattenedArray(h_flows, size, size, "flows:");
#endif

  if (iseed == 0) 
    srand(time(NULL));
  else
    srand(iseed);

  int *h_sol;
  genInitSolutions(&h_sol, size, solns);

#ifdef PROFILE
  cudaEventRecord(start,0);
#endif

  // allocate GPU memory
  int *d_result = NULL;
  int *d_bestsofar = NULL, *d_bestcostsofar = NULL, *d_newarray = NULL;
  short *d_dist = NULL,*d_flows = NULL;
  int *d_sol = NULL;

  gpuErrchk(cudaMalloc((void **) &d_bestcostsofar, solns * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_bestsofar, solns * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_newarray, solns * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_result, solns * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_dist, size * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_flows, size * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_sol, solns * size * sizeof(int)));

  // allocated CPU memory
  int *h_result = (int *) malloc(solns * sizeof(int));
  int *h_bestsofar = (int *) malloc(solns * size * sizeof(int));
  int *h_bestcostsofar = (int *) malloc(solns * sizeof(int));

  // copy host data to device 	
  gpuErrchk(cudaMemcpy(d_dist,h_dist,size*size*sizeof(int),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_flows,h_flows,size*size*sizeof(int),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_sol,h_sol,solns*size*sizeof(int),cudaMemcpyHostToDevice));

  // determine thread configurations 
  int threadsPerBlock = 128;
  int blockPerGrid = (solns + threadsPerBlock - 1) / threadsPerBlock;

#if DEBUG 
  cout << "blocks: " << blockPerGrid << " " << endl;
  cout << "threads: " << threadsPerBlock << " " << endl;
#endif

  twoOptShared<<< blockPerGrid, threadsPerBlock >>>(d_dist, d_flows, d_sol, size, solns, 
					   d_result, d_bestsofar, d_bestcostsofar, d_newarray);

  // copy device data to host
  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaMemcpy(h_bestcostsofar, d_bestcostsofar, solns * sizeof(int),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(h_bestsofar, d_bestsofar, solns * size * sizeof(int),cudaMemcpyDeviceToHost)); 

#ifdef PROFILE
  cudaEventRecord(stop,0);
  cudaEventSynchronize(start);
  cudaEventSynchronize(stop);
#endif 

#ifdef DEBUG  
  printFlattenedArray(h_bestsofar, solns, size, "Best permutations");
  printRegularArray(h_bestcostsofar, solns, "Best costs:");
#endif

  // find the best among all the solutions resturneed 
  int minIndex = 0;
  int minCost = findMin(h_bestcostsofar, solns, &minIndex);

  //  print results
  cout << "problem size: " << size << endl;
  cout << "best cost: " << minCost << endl;
  cout << "best solution: " ;
  for(int i = 0; i < size; i++)
    cout << h_bestsofar[minIndex * size + i] << " ";
  cout << endl;

#ifdef PROFILE
  cpu_endTime = clock();
  cpu_ElapseTime= ((cpu_endTime - cpu_startTime) / (double) CLOCKS_PER_SEC);

  fprintf(stdout, "total exec time (s):\t %2.2f\n", cpu_ElapseTime); 
  cudaEventElapsedTime(&ctime, start , stop);
  fprintf(stdout, "kernel exec:\t %2.2f\n", (ctime / 1000));

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
#endif


  free(array);
  free(h_dist);
  free(h_flows);
  free(h_sol);
  free(h_result);
  free(h_bestsofar);
  free(h_bestcostsofar);
  cudaFree(d_dist);
  cudaFree(d_flows);
  cudaFree(d_sol);
  cudaFree(d_bestcostsofar);
  cudaFree(d_newarray);
  cudaFree(d_result);
  cudaFree(d_bestsofar);
  return 0;
}
