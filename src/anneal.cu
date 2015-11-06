#include<iostream>
#include<fstream>
#include<stdio.h>
#include<cstdlib>
#include<qapio.h>
#include<2optlib.h>

#include<curand_kernel.h>

using namespace std;

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

__device__ float accProb(int cost, int newCost, float t) {
  float diff = (newCost - cost)/cost;
  float prob = exp(diff/t);
  if (prob > 1) 
    return 1.0;
  else 
    return prob;
}

__global__ void anneal(int *d_dist,int *d_flows,int *d_sol, int nsize, int row,
		       int *d_bestsofar,int *d_bestcostsofar,int *d_newarray) {	


  int tidx = threadIdx.x + blockDim.x * blockIdx.x;

  // number of inital solutions should equal number of threads 
  // if tidx > number of solutions, then something is wrong 
  if (tidx >= row) 
    return; 

  curandState_t state;
  curand_init(1023, tidx, 0, &state);


  float t = 1.0;
  float tMin = 0.001;
  int tMax = t/tMin;
  float alpha = 0.9;
  int neighbors = nsize;

  
  int index = tidx * nsize;
  int dcost = 0, ecost = 0;
  int delta = 0;
  int tcost;
  
  // calculate cost of initial solution
  d_bestcostsofar[tidx]  = initCost(d_flows, d_dist, d_sol, nsize, tidx);
  copy(d_bestsofar,d_sol,nsize,tidx);

  // search until temp threashold is reached 
  for(int n = 0; n < tMax; n++) {	
    for (int i = 0; i < neighbors; i++) {     
      // generate neighbor permuation by swapping a pair of units 
      copy(d_newarray, d_sol, nsize, tidx);

      int k = curand(&state) % nsize;
      int j = curand(&state) % nsize;
      
      swap(&d_newarray[tidx * nsize + k], &d_newarray[tidx * nsize + j]);

      // calculate cost of neighbor
      int kUnit = d_sol[index + k] - 1;
      int jUnit = d_sol[index + j] - 1;
    
      dcost = (d_dist[j * nsize + k] - d_dist[k * nsize + j]) 
	* (d_flows[kUnit * nsize + jUnit] - d_flows[jUnit * nsize + kUnit]);
    
      ecost = (d_dist[j * nsize + j] - d_dist[k * nsize + k])
	* (d_flows[kUnit * nsize + kUnit] - d_flows[jUnit * nsize + jUnit]);
    
      delta = dcost + ecost + neighborCost(d_flows, d_dist, d_sol, nsize, tidx, k, j);
      tcost = d_bestcostsofar[tidx] + delta;

      float rand = curand_uniform(&state);
      // update results if a permuation with lower cost is found
      if (accProb(d_bestcostsofar[tidx], tcost, t) > rand) {
	d_bestcostsofar[tidx] = tcost;
	copy(d_bestsofar, d_newarray, nsize, tidx);
	copy(d_sol, d_newarray, nsize, tidx);
      }
    }
    t = t * alpha;
  }	    
	   
  return;
}


int main(int argc,char *argv[]) {

#ifdef PROFILE
  clock_t cpu_startTime, cpu_endTime; 
  double cpu_ElapseTime=0;

  cpu_startTime = clock();
#endif

  if (argc != 4) {
    cout << "usage: " << endl;
    cout << "\t./2opt datafile solns randseed" << endl;
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

  // read data from file 
  int size;
  int **array;
  readData(filename, &array, &size);
 
  // split and flatten matrix
  int *h_flows, *h_dist;
  splitAndFlattenInt(&h_flows, &h_dist, array, size);

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

#ifdef DEBUG
  printFlattenedArray(h_sol, size, size, "initial solutions:");
#endif

#ifdef PROFILE
  cudaEvent_t start , stop;
  float ctime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,0);
#endif

  // allocate GPU memory
  int *d_bestsofar = NULL, *d_bestcostsofar = NULL, *d_newarray = NULL;
  int *d_dist = NULL,*d_flows = NULL ,*d_sol = NULL;

  gpuErrchk(cudaMalloc((void **) &d_bestcostsofar, solns * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_bestsofar, solns * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_newarray, solns * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_dist, size * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_flows, size * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_sol, solns * size * sizeof(int)));

  // allocated CPU memory
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

  anneal<<< blockPerGrid, threadsPerBlock >>>(d_dist, d_flows, d_sol, size, solns, 
					      d_bestsofar, d_bestcostsofar, d_newarray);

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
  free(h_bestsofar);
  free(h_bestcostsofar);
  cudaFree(d_dist);
  cudaFree(d_flows);
  cudaFree(d_sol);
  cudaFree(d_bestcostsofar);
  cudaFree(d_newarray);
  cudaFree(d_bestsofar);

  return 0;
}
