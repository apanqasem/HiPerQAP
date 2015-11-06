#include<iostream>
#include<fstream>
#include<stdio.h>
#include<cstdlib>

#include<qapio.h>
#include<2optlib.h>

using namespace std;
#define width 4950

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

__device__ void copyNoDestOffset(int *dest, int *src, int size, int tidx) {
  int offset = tidx * size;
  for(int i = 0; i < size; i++)
    dest[i] = src[offset + i];
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


__global__ void child_kernel(int tidx,int nsize,int *d_dist,int *d_flows,int *d_sol,int *d_bestcostsofar,int *d_bestsofar,int *d_result,int row,int *d_pos)
{
  int ipos=0,jpos=0,ir=0,ik=0;
  int tx = threadIdx.x;
  int xj = ((nsize *(nsize-1))/2)/(nsize-1);
  int d_tmpsol[100];
  __shared__ int pos[width * 2];
  if(tx<(nsize-1))
    {	ik =(nsize-2)*tx;
      for(int j=0;j<xj;j++)
	{
	  pos[(tx*2)+ik] = d_pos[(tx*2)+ik];
	  pos[(tx*2)+(ik+1)] = d_pos[(tx*2)+(ik+1)];
	  ik = ik + 2;
	}
    }
  __syncthreads();
  if(tx<nsize-1)
    {
      ir = (nsize-2)*tx;
      for(int j=0;j<xj;j++)
	{	
	  ipos=0;
	  jpos=0;
	  copyNoDestOffset(d_tmpsol,d_sol,nsize,tidx);
	  //swap(&d_tmpsol[pos[(tx*2)+ir]],&d_tmpsol[pos[(tx*2)+(ir+1)]]);
	  ipos =  pos[(tx*2)+ir];
	  jpos =  pos[(tx*2)+(ir+1)];
	  //printf("parent id%d\t child id %d\t iposition and jpos to swap %d %d \n",tidx,tx,ipos,jpos);

	  int dcost=0,ecost=0,fcost=0;
	  dcost = (d_dist[jpos*nsize+ipos] - d_dist[ipos*nsize+jpos])*(d_flows[(d_sol[(tidx * nsize)+ipos]-1) * nsize + (d_sol[(tidx * nsize)+jpos] - 1)] - d_flows[(d_sol[(tidx * nsize)+jpos]-1) * nsize + (d_sol[(tidx * nsize)+ipos] - 1)]);
	  ecost = (d_dist[jpos*nsize+jpos] - d_dist[ipos*nsize+ipos])*(d_flows[(d_sol[(tidx * nsize)+ipos]-1) * nsize + (d_sol[(tidx * nsize)+ipos] - 1)] - d_flows[(d_sol[(tidx * nsize)+jpos]-1) * nsize + (d_sol[(tidx * nsize)+jpos] - 1)]);
	  fcost = dcost + ecost;
	  int totcost=0,delta=0,tcost=0;
	  totcost = neighborCost(d_flows,d_dist,d_sol,nsize,tidx,ipos,jpos);
	  //__syncthreads();
	  //tcost=calculating(d_flows,d_dist,d_newarray,nsize,tidx);									
	  delta = fcost + totcost;
	  tcost = d_result[tidx] + delta;
	  if(tcost<d_bestcostsofar[tidx])
	    {	
	      d_bestcostsofar[tidx]=tcost;
	      swap(&d_tmpsol[ipos],&d_tmpsol[jpos]);
	      for(int j=0;j<nsize;j++)
		{
		  d_bestsofar[tidx * nsize + j] = d_tmpsol[j];
		}
	    }
	  ir = ir +2;	
	}//end of j

    }//end of k
  free(d_tmpsol);		
}

__device__ void diversification(int *d_sol,int *d_bestsofar,int nsize,int tidx,int row)
{
  int offset=25,pos,istart;
  //for(int i=0;i<row;i++)
  //{
  pos=0;
  istart =0;
  for(int start=offset;start>=0;start--)
    {
      istart = start;
      while(istart<nsize)
	{
	  d_sol[tidx*nsize+pos] = d_bestsofar[tidx*nsize+istart];
	  pos=pos+1;
	  if(istart!=0)
	    istart = istart + offset;
	  else
	    break;
	}

    }
  //}

}

__global__ void max(int *d_dist,int *d_flows,int *d_sol,int nsize,int row,int *d_result,int *d_bestsofar,int *d_bestcostsofar,int *d_pos)
{	int totalcost =0;
  int tidx = threadIdx.x+blockDim.x * blockIdx.x;
  if(tidx < row)
    {	
      totalcost=initCost(d_flows,d_dist,d_sol,nsize,tidx);
      d_result[tidx] = totalcost;
      d_bestcostsofar[tidx]=d_result[tidx];
      //      copy(d_sol,d_bestsofar,nsize,tidx);
      copy(d_bestsofar,d_sol,nsize,tidx);
    }
  __syncthreads();
  for(int l=0;l<2;l++)
    {	
      if(tidx<row)
	{
	  int threadsPerBlock= (nsize * (nsize - 1))/nsize;
	  //printf("child kernel number of threads %d:",threadsPerBlock);
	  child_kernel<<<1,threadsPerBlock>>>(tidx,nsize,d_dist,d_flows,d_sol,d_bestcostsofar,d_bestsofar,d_result,row,d_pos);			
	  if (cudaSuccess != cudaGetLastError())
	    {
	      return;
	    }
	  //cudaDeviceSynchronize();
	  //cudaError_t err = cudaGetLastError();
	  //if (err != cudaSuccess) printf("!");
	  // wait for child to complete
	  if (cudaSuccess != cudaDeviceSynchronize()) {
	    return;
	  }
	  for(int jc=0;jc<nsize;jc++)
	    {
	      d_sol[tidx*nsize+jc]= d_bestsofar[tidx*nsize+jc];
	    }
	  d_result[tidx] = d_bestcostsofar[tidx];
	  //diversification(d_sol,d_bestsofar,nsize,tidx,row);
	  //totalcost=calculating(d_flows,d_dist,d_sol,nsize,tidx);
	  //d_result[tidx] = totalcost;


	}
    }
  /*__syncthreads();
    if(tidx<row)
    {	
    diversification(d_sol,d_bestsofar,nsize,tidx,row);
    totalcost=calculating(d_flows,d_dist,d_sol,nsize,tidx);
    d_result[tidx] = totalcost;
    }
    __syncthreads();
    for(int l=0;l<nsize;l++)
    {	
    if(tidx<row)
    {
    int threadsPerBlock= (nsize * (nsize - 1))/nsize;
    //printf("child kernel number of threads %d:",threadsPerBlock);
    child_kernel<<<1,threadsPerBlock>>>(tidx,nsize,d_dist,d_flows,d_sol,d_bestcostsofar,d_bestsofar,d_result,row,d_pos);			
    if (cudaSuccess != cudaGetLastError())
    {
    return;
    }
    //cudaDeviceSynchronize();
    //cudaError_t err = cudaGetLastError();
    //if (err != cudaSuccess) printf("!");
    // wait for child to complete
    if (cudaSuccess != cudaDeviceSynchronize()) {
    return;
    }
    for(int jc=0;jc<nsize;jc++)
    {
    d_sol[tidx*nsize+jc]= d_bestsofar[tidx*nsize+jc];
    }
    d_result[tidx] = d_bestcostsofar[tidx];
    }
    }
  */	

  /*for(int i=1;i<row;i*=2)
    {
    if(tidx % (2 * i) == 0)
    {
    if(d_bestcostsofar[tidx]>d_bestcostsofar[tidx+i])
    {
    d_bestcostsofar[tidx] = d_bestcostsofar[tidx+i];
    }
    __syncthreads();
    }
    }

    printf("best cost %d:",d_bestcostsofar[0]);
  */
}

#if 0
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
  int size_B1 = (size*(size-1))/2 * 2;
  int mem_size_B1 = sizeof(int) * size_B1;
  int *h_pos = (int *)malloc(mem_size_B1);
  int l=0;
  for(int i=0;i<size-1;i++)
    {
      int q =l;
      for(int j=i+1;j<size;j++)
	{
	  int ipos=i;
	  int jpos=j;
	  h_pos[q*2+0] = ipos;
	  h_pos[q*2+1] = jpos;
	  q++;
	}
      l=q;
    }

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
  int *d_result = NULL;
  int *d_bestsofar = NULL, *d_bestcostsofar = NULL, *d_newarray = NULL;
  int *d_dist = NULL,*d_flows = NULL ,*d_sol = NULL,*d_pos=NULL;

  gpuErrchk(cudaMalloc((void **) &d_bestcostsofar, solns * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_bestsofar, solns * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_newarray, solns * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_result, solns * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_dist, size * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_flows, size * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_sol, solns * size * sizeof(int)));
  gpuErrchk(cudaMalloc((void **) &d_pos,(size * (size - 1))/2 * 2 * sizeof(int)));

  // allocated CPU memory
  int *h_result = (int *) malloc(solns * sizeof(int));
  int *h_bestsofar = (int *) malloc(solns * size * sizeof(int));
  int *h_bestcostsofar = (int *) malloc(solns * sizeof(int));

  // copy host data to device 	
  gpuErrchk(cudaMemcpy(d_dist,h_dist,size*size*sizeof(int),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_flows,h_flows,size*size*sizeof(int),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_sol,h_sol,solns*size*sizeof(int),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_pos,h_pos,(size*(size-1))/2 * 2 * sizeof(int),cudaMemcpyHostToDevice) );

  // determine thread configurations 
  int threadsPerBlock = 1024;
  int blockPerGrid = (solns + threadsPerBlock - 1) / threadsPerBlock;

#if DEBUG 
  cout << "blocks: " << blockPerGrid << " " << endl;
  cout << "threads: " << threadsPerBlock << " " << endl;
#endif

  max<<< blockPerGrid, threadsPerBlock >>>(d_dist, d_flows, d_sol, size, solns, 
					   d_result, d_bestsofar, d_bestcostsofar, d_newarray);

  // copy device data to host
  gpuErrchk(cudaPeekAtLastError());
  if (cudaSuccess != cudaGetLastError()) {
    return 1;
  }
  // wait for parent to complete
  if (cudaSuccess != cudaDeviceSynchronize()) {
    return 2;
  }

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
  fprintf(stdout, "kernel exec time (s):\t %2.2f\n", (ctime / 1000));

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

#endif

int main(int argc,char *argv[])
{

  int arraySizeX,arraySizeY,size,num,seed,a=0,b=0;	
  clock_t cpu_startTime, cpu_endTime;
  double cpu_ElapseTime=0;
  cpu_startTime = clock();
  ifstream input;
  cudaEvent_t start , stop;
  float ctime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  int solns = atoi(argv[2]);

  //  cout<<"input file name:"<<argv[1]<<endl; 

  int iseed = atoi(argv[3]);
  input.open(argv[1],ios::in);

  if(!input.is_open())
    {
      cout<<"error opening file";

    }
  //reading the size,seed from the file
  input>>size>>seed;

  int nsize;
  nsize=size;
  arraySizeX=2*nsize;
  arraySizeY=nsize;
  //declaring array to copy the matrix from file into array    	
  int** array;
  array = (int**) malloc(arraySizeX*sizeof(int*));
  for (int i = 0; i < arraySizeX; i++)
    array[i] = (int*) malloc(arraySizeY*sizeof(int));
  for(int row=0;row<(arraySizeX);row++)
    {
      for(int col=0;col<arraySizeY;col++){
	array[row][col]=0;
      }//end col for
    }//end row for
  //flatten array dist nd flows declarations
  int size_A = nsize * nsize;
  int mem_size_A = sizeof(int) * size_A;
  int *h_dist = (int *)malloc(mem_size_A);
  int *h_flows = (int *)malloc(mem_size_A);

  int *h_a = (int *)malloc(mem_size_A);

  for(int i=0;i<nsize;i++)
    {
      for(int j=0;j<nsize;j++)
	{
	  h_dist[i *nsize +j] = 0;
	  h_flows[i* nsize +j] = 0;
	}
    }
  while(!input.eof())
    {
      input>>num;
      if(b==nsize)
	{
	  a++;
	  b=0;
	}
      //a->row,b->col
      if(a!=(nsize*2) && b!=nsize)
	{
	  array[a][b]=num;
	  b++;
	}//end if
    }// end-while
  input.close();
  for(int row=0;row<nsize;row++)
    {
      for(int col=0;col<nsize;col++)
	{
	  h_flows[row *nsize + col]=array[row][col];
	}
    }

  //storing in dist_dup array
  int irow=0;
  for(int row=nsize;row<nsize*2;row++)
    {
      int icol=0;
      for(int col=0;col<nsize;col++)
	{
	  h_dist[irow *nsize + icol]=array[row][col];
	  icol++;

	}
      irow++;
    }

  srand(iseed);
  int *init_sol,j; //,solns=6144;
  init_sol = (int *)malloc(nsize * sizeof(int));
  for(int i=0;i<nsize;i++)
    init_sol[i] = i+1;
  int size_B = (solns) * nsize;
  int mem_size_B = sizeof(int) * size_B;
  int *h_sol = (int *)malloc(mem_size_B);
  for(int k=0;k<solns;k++)
    {
      for(int i=nsize-1;i>0;i--)
	{
	  j= rand() % (i+1);
	  int temp = init_sol[i];
	  init_sol[i] = init_sol[j];
	  init_sol[j]=temp;
	}
      for(int l=0;l<nsize;l++)
	{
	  h_sol[k *nsize + l] = init_sol[l];
	}
    }
				
  int size_B1 = (nsize*(nsize-1))/2 * 2;
  int mem_size_B1 = sizeof(int) * size_B1;
  int *h_pos = (int *)malloc(mem_size_B1);
  int l=0;
  for(int i=0;i<nsize-1;i++)
    {
      int q =l;
      for(int j=i+1;j<nsize;j++)
	{
	  int ipos=i;
	  int jpos=j;
	  h_pos[q*2+0] = ipos;
	  h_pos[q*2+1] = jpos;
	  q++;
	}
      l=q;
    }
  int *h_result;			
  int *d_result;
  int *h_bestsofar,*h_bestcostsofar;
  h_result = (int *)malloc(solns * sizeof(int));
  h_bestsofar = (int *)malloc(solns * nsize * sizeof(int));
  h_bestcostsofar = (int *)malloc(solns * sizeof(int));
  cudaEventRecord(start,0);
  int *d_bestsofar=NULL,*d_bestcostsofar=NULL;
  gpuErrchk( cudaMalloc((void **)&d_bestcostsofar,solns * sizeof(int)) );
  gpuErrchk( cudaMalloc((void **)&d_bestsofar,solns * nsize * sizeof(int)) );
  gpuErrchk( cudaMalloc((void **)&d_result,solns * sizeof(int)) );

  // declaring device array and allocating memory on gpu
  int *d_dist = NULL,*d_flows = NULL ,*d_sol = NULL,*d_pos=NULL;
  gpuErrchk( cudaMalloc((void **)&d_dist,nsize*nsize*sizeof(int)) );
  gpuErrchk( cudaMalloc((void **)&d_flows,nsize*nsize*sizeof(int)) );
  gpuErrchk( cudaMalloc((void **)&d_pos,(nsize*(nsize-1))/2 * 2 * sizeof(int)));
  gpuErrchk( cudaMalloc((void **)&d_sol,solns*nsize*sizeof(int)) );
  //copying arrays from host to device 	
  gpuErrchk( cudaMemcpy(d_dist,h_dist,nsize*nsize*sizeof(int),cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(d_flows,h_flows,nsize*nsize*sizeof(int),cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(d_sol,h_sol,solns*nsize*sizeof(int),cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(d_pos,h_pos,(nsize*(nsize-1))/2 * 2 * sizeof(int),cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(d_bestcostsofar,h_bestcostsofar,solns * sizeof(int),cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(d_bestsofar,h_bestsofar,solns * nsize* sizeof(int), cudaMemcpyHostToDevice));
  gpuErrchk( cudaMemcpy(d_result,h_result,solns * sizeof(int), cudaMemcpyHostToDevice) );
  //cuda kernel call 
  int threadsPerBlock=256;
  int blockPerGrid = (solns + threadsPerBlock - 1) / threadsPerBlock;
  //int smSize= threadsPerBlock *nsize*nsize*sizeof(int);

#if 1
  cout<<"number of initial solutions: "<<solns<<endl;
  cout<<"number of blocks: "<<blockPerGrid<<" "<<endl;
  cout<<"number of threads: "<<threadsPerBlock<<" "<<endl;
#endif

  max<<<blockPerGrid,threadsPerBlock>>>(d_dist,d_flows,d_sol,nsize,solns,d_result,d_bestsofar,d_bestcostsofar,d_pos);
  gpuErrchk( cudaPeekAtLastError() );
  if (cudaSuccess != cudaGetLastError()) {
    return 1;
  }
  // wait for parent to complete
  if (cudaSuccess != cudaDeviceSynchronize()) {
    return 2;
  }
  gpuErrchk( cudaMemcpy(h_result,d_result,solns * sizeof(int),cudaMemcpyDeviceToHost) );
  gpuErrchk( cudaMemcpy(h_bestcostsofar,d_bestcostsofar,solns * sizeof(int),cudaMemcpyDeviceToHost) );
  gpuErrchk( cudaMemcpy(h_bestsofar, d_bestsofar , solns * nsize* sizeof(int),cudaMemcpyDeviceToHost )); 
  //gpuErrchk( cudaMemcpy(temp_sol, d_tmpsol , (solns*(nsize-1)) * nsize* sizeof(int),cudaMemcpyDeviceToHost ));
  gpuErrchk( cudaMemcpy(h_sol,d_sol , (solns) * nsize* sizeof(int),cudaMemcpyDeviceToHost ));
  cudaEventRecord(stop,0);
  cudaEventSynchronize(start);
  cudaEventSynchronize(stop);

#ifdef DEBUG
  cout<<"cost of best solutions sofar and best solutin array's:"<<endl;
  for(int i=0;i<solns;i++)
    {
      for(j=0;j<nsize;j++)
        {
	  cout<<h_bestsofar[i*nsize+j]<<" ";
        }
      cout<<h_bestcostsofar[i]<<" ";
      cout<<endl;
    }

  cout<<"cost of diversified initial sols:";
    for(int i=0;i<solns;i++)
    {
    for(j=0;j<nsize;j++)
    {
    cout<<h_sol[i*nsize+j]<<" ";
    }
    cout<<h_result[i]<<" ";
    cout<<endl;
    }

    int offset=5,pos,istart;
    for(int i=0;i<solns;i++)
    {
    pos=0;
    istart =0;
    for(int start=offset;start>=0;start--)
    {
    istart = start;
    while(istart<nsize)
    {
    h_newarray[i*nsize+pos] = h_bestsofar[i*nsize+istart];
    pos=pos+1;
    if(istart!=0)
    istart = istart + offset;
    else
    break;
    }

    }
    }
    cout<<"newarray sol:";
    cout<<endl;
    for(int i=0;i<(solns);i++)
    {
    for(j=0;j<nsize;j++)
    {
    cout<<h_newarray[i*nsize+j]<<" ";
    }
    cout<<endl;
    }
  cout<<endl;
#endif

  cout<<"best solution:"<<endl;
  int temp=0,lrow=0;
  temp = h_bestcostsofar[0];
  for(int i=1;i<solns;i++)
    {
      if(h_bestcostsofar[i]<temp)
	{
	  temp = h_bestcostsofar[i];
	  lrow=i;
	}
    }
  for(int j=0;j<nsize;j++)
    {
      cout<<h_bestsofar[lrow*nsize+j]<<" ";
    }
  cout<<endl;


  cout<< "best cost:\t" << temp << endl;

  cpu_endTime = clock();
  cpu_ElapseTime= ((cpu_endTime - cpu_startTime) / (double) CLOCKS_PER_SEC);

  fprintf(stdout, "total exec time (s):\t %2.2f\n", cpu_ElapseTime); 
  cudaEventElapsedTime(&ctime, start , stop);

  fprintf(stdout, "kernel exec time (s):\t %2.2f\n", (ctime / 1000));


  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  cout<<endl;
  free(array);
  free(h_dist);
  free(h_flows);
  free(init_sol);
  free(h_sol);
  free(h_pos);
  free(h_result);
  free(h_bestsofar);
  free(h_bestcostsofar);
  cudaFree(d_dist);
  cudaFree(d_flows);
  cudaFree(d_sol);
  cudaFree(d_bestcostsofar);
  cudaFree(d_result);
  cudaFree(d_bestsofar);
  cudaFree(d_pos);
  return 0;
}
