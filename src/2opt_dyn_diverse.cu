#include<iostream>
#include <fstream>
#include<stdio.h>
#include<cstdlib>
#include <sstream>
#include<string>
using namespace std;
#define width 4950
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

//device function
__device__ int calculating(int *d_flows,int *d_dist,int *d_sol,int nsize,int tidx)
{
	int calcost =0;
	for(int i=0;i<nsize-1;i++)
          {
                for(int j=i+1;j<nsize;j++)
                {
                        calcost = calcost + ( d_flows[ ( d_sol[(tidx* nsize)+i]-1) *nsize + (d_sol[(tidx* nsize)+j]-1)]) * d_dist[i*nsize +j];

                }
          }
 		for(int k=1;k<nsize;k++)
                {
                        for(int l=0;l<k;l++)
                        {

				calcost = calcost + d_flows[(d_sol[(tidx* nsize)+k]-1) *nsize + (d_sol[(tidx* nsize)+l]-1)] * d_dist[k *nsize + l];	
			}
		}
return calcost;
}

__device__ int calculate(int *d_flows,int *d_dist,int *d_sol,int nsize,int tidx,int i,int j)
{
	int ccost=0,gcost=0,hcost=0;
	for(int k=0;k<nsize;k++)
	{
		if(k!=i && k!=j)
		{
			gcost = (d_dist[j*nsize+k] - d_dist[i*nsize+k]) *(d_flows[(d_sol[(tidx * nsize)+i]-1) * nsize + (d_sol[(tidx * nsize)+k] - 1)] - d_flows[(d_sol[(tidx * nsize)+j]-1) * nsize + (d_sol[(tidx * nsize)+k] - 1)]); 
			hcost = (d_dist[k*nsize+j] - d_dist[k*nsize+i]) *(d_flows[(d_sol[(tidx * nsize)+k]-1) * nsize + (d_sol[(tidx * nsize)+i] - 1)] - d_flows[(d_sol[(tidx * nsize)+k]-1) * nsize + (d_sol[(tidx * nsize)+j] - 1)]);
			ccost = ccost + (gcost + hcost);
		}

	}
 return ccost;
}

__device__ void copy(int *d_sol,int *d_newarray,int nsize,int tidx)
{
	
		for(int j=0;j<nsize;j++)
		{
			d_newarray[tidx * nsize + j] = d_sol[tidx * nsize + j];
		}
	
}
__device__ void copy1(int *d_tmpsol,int *d_sol,int nsize,int row,int tidx)
{
                     for(int j=0;j<nsize;j++)
                     {
                          d_tmpsol[j] = d_sol[tidx * nsize + j];
                     }
}

__device__ void copy2(int *d_tmpsol,int *d_newt,int nsize)
{
	for(int j=0;j<nsize;j++)
	{
			//d_newarray[((tidx *(nsize-1))+tx)* nsize + j] = d_tmpsol[j];
			d_newt[j] = d_tmpsol[j];
	}
	
}
__device__ void swap(int *a,int *b)
{
		int temp=0;
		temp = *a;
		*a = *b;
		*b = temp;
	
}
__device__ void least(int *d_newarray,int nsize,int row,int tidx,int *d_divresult,int *d_pos)
{

int temp=0,lrow=0,r;
r = (nsize * (nsize-1))/2;
temp = d_divresult[tidx*r];
for(int i=1;i<r;i++)
{
	if(d_divresult[tidx*r+i]<temp)
	{
		temp = d_divresult[tidx*r+i];
		lrow=i;
	}
}
//printf("%d least position in array",lrow);
swap(&d_newarray[tidx*nsize+d_pos[lrow*2]],&d_newarray[tidx*nsize+d_pos[lrow*2+1]]);
}
__global__ void child_kernel(int tidx,int nsize,int *d_dist,int *d_flows,int *d_sol,int *d_bestcostsofar,int *d_bestsofar,int *d_result,int row,int *d_pos,int *d_newarray,int *d_divresult,int *d_frequency)
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
copy1(d_tmpsol,d_sol,nsize,row,tidx);
if(tx<nsize-1)
{
 ir = (nsize-2)*tx;
for(int j=0;j<xj;j++)
{	
	ipos=0;
	jpos=0;
	copy1(d_tmpsol,d_sol,nsize,row,tidx);
	//swap(&d_tmpsol[pos[(tx*2)+ir]],&d_tmpsol[pos[(tx*2)+(ir+1)]]);
	ipos =  pos[(tx*2)+ir];
    	jpos =  pos[(tx*2)+(ir+1)];
	//printf("parent id%d\t child id %d\t iposition and jpos to swap %d %d \n",tidx,tx,ipos,jpos);
	int dcost=0,ecost=0,fcost=0;
	dcost = (d_dist[jpos*nsize+ipos] - d_dist[ipos*nsize+jpos])*(d_flows[(d_sol[(tidx * nsize)+ipos]-1) * nsize + (d_sol[(tidx * nsize)+jpos] - 1)] - d_flows[(d_sol[(tidx * nsize)+jpos]-1) * nsize + (d_sol[(tidx * nsize)+ipos] - 1)]);
	ecost = (d_dist[jpos*nsize+jpos] - d_dist[ipos*nsize+ipos])*(d_flows[(d_sol[(tidx * nsize)+ipos]-1) * nsize + (d_sol[(tidx * nsize)+ipos] - 1)] - d_flows[(d_sol[(tidx * nsize)+jpos]-1) * nsize + (d_sol[(tidx * nsize)+jpos] - 1)]);
	fcost = dcost + ecost;
	int totcost=0,delta=0,tcost=0;
	totcost = calculate(d_flows,d_dist,d_sol,nsize,tidx,ipos,jpos);
	//__syncthreads();
	//tcost=calculating(d_flows,d_dist,d_newarray,nsize,tidx);									
	delta = fcost + totcost;
	tcost = d_result[tidx] + delta;
	d_divresult[tidx*(nsize*(nsize-1)/2)+(tx*xj+j)]=tcost;
	if(tcost<d_bestcostsofar[tidx])
	{
		d_bestcostsofar[tidx]=tcost;
		swap(&d_tmpsol[ipos],&d_tmpsol[jpos]);
		for(int j=0;j<nsize;j++)
		{
			d_bestsofar[tidx * nsize + j] = d_tmpsol[j];
			//d_newarray[tidx*nsize+j] = d_tmpsol[j];
		}
	}
	
ir = ir +2;	
}//end of j

}//end of k
//least(d_newarray,nsize,row,tidx,d_divresult,d_pos);
free(d_tmpsol);		
}

__device__ void diversification(int *d_sol,int *d_newarray,int nsize,int tidx,int l)
{
	/*if(l>nsize/2)
	{
		l = l/4*10;
	}*/
	int offset=9;
	int pos,istart;
	//for(int i=0;i<row;i++)
	//{
		pos=0;
		istart =0;
		for(int start=offset;start>=0;start--)
		{
			istart = start;
			while(istart<nsize)
			{
				d_sol[tidx*nsize+pos] = d_newarray[tidx*nsize+istart];
				pos=pos+1;
				if(istart!=0)
					istart = istart + offset;
				else
					break;
			}

		}
	//}

}

__global__ void max(int *d_dist,int *d_flows,int *d_sol,int nsize,int row,int *d_result,int *d_bestsofar,int *d_bestcostsofar,int *d_pos,int *d_newarray,int *d_divresult,int *d_frequency)
{	int totalcost =0,divcost=0;//lrow=0;
	 int tidx = threadIdx.x+blockDim.x * blockIdx.x;
	if(tidx < row)
	{	
		totalcost=calculating(d_flows,d_dist,d_sol,nsize,tidx);
	 	d_result[tidx] = totalcost;
		d_bestcostsofar[tidx]=d_result[tidx];
		copy(d_sol,d_bestsofar,nsize,tidx);
		copy(d_sol,d_newarray,nsize,tidx);
	}
__syncthreads();
	/*for(int l=0;l<nsize;l++)
	{	
	    if(tidx<row)
		{
			int threadsPerBlock= (nsize * (nsize - 1))/nsize;
			//printf("child kernel number of threads %d:",threadsPerBlock);
			child_kernel<<<1,threadsPerBlock>>>(tidx,nsize,d_dist,d_flows,d_sol,d_bestcostsofar,d_bestsofar,d_result,row,d_pos,d_newarray,d_divresult,d_frequency);			
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
			if(l<nsize/2)
			{
			for(int jc=0;jc<nsize;jc++)
                	{
                    		d_sol[tidx*nsize+jc]= d_bestsofar[tidx*nsize+jc];
                	}
                       	d_result[tidx] = d_bestcostsofar[tidx];
			}
			else
			{
			diversification(d_sol,d_bestsofar,nsize,tidx,l);
			divcost=calculating(d_flows,d_dist,d_sol,nsize,tidx);
			d_result[tidx] = divcost;
			}
				
		}
   	}*/
	
	/*__syncthreads();
	if(tidx<row)
	{	
		diversification(d_newarray,d_bestsofar,nsize,tidx);
		divcost=calculating(d_flows,d_dist,d_newarray,nsize,tidx);
	 	d_result[tidx] = divcost;
	}
	__syncthreads();
	*/
	for(int l=0;l<nsize*10;l++)
	{	
	    if(tidx<row)
		{
			int threadsPerBlock= (nsize * (nsize - 1))/nsize;
			//printf("child kernel number of threads %d:",threadsPerBlock);
			child_kernel<<<1,threadsPerBlock>>>(tidx,nsize,d_dist,d_flows,d_sol,d_bestcostsofar,d_bestsofar,d_result,row,d_pos,d_newarray,d_divresult,d_frequency);				
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
		
			if(l%5==0)
			//if(l<(nsize/2+(nsize/4)))
			{
				if(tidx>row/2)
				{
					for(int jc=0;jc<nsize;jc++)
					{
						d_sol[tidx*nsize+jc]= d_bestsofar[tidx*nsize+jc];
					}
                    			d_result[tidx] = d_bestcostsofar[tidx];
				}
				else
				{
					diversification(d_sol,d_newarray,nsize,tidx,l);
					divcost=calculating(d_flows,d_dist,d_sol,nsize,tidx);
					d_result[tidx] = divcost;
				}
			}
			else
			{	
				copy(d_sol,d_newarray,nsize,tidx);
				least(d_newarray,nsize,row,tidx,d_divresult,d_pos);
				/*if(l%2==0)
				{
					diversification(d_sol,d_newarray,nsize,tidx,l);
					divcost=calculating(d_flows,d_dist,d_sol,nsize,tidx);
					d_result[tidx] = divcost;
				
				}*/
				//else
				//{
						for(int jc=0;jc<nsize;jc++)
                		{
                    			d_sol[tidx*nsize+jc]= d_newarray[tidx*nsize+jc];
                		}
                               divcost=calculating(d_flows,d_dist,d_sol,nsize,tidx);
								d_result[tidx] = divcost;
				//}
			}
	}
   }	
/*
temp = d_bestcostsofar[0];
//for(int i=1;i<row;i++)
if(tidx<row)
{
        if(d_bestcostsofar[tidx]<temp)
        {
		temp = d_bestcostsofar[tidx];
                lrow=tidx;
        }
__syncthreads();
}
 for(int j=0;j<nsize;j++)
 {
                printf("%d \t",d_bestsofar[lrow*nsize+j]);
 }
printf("\n");
printf("best cost %d:",temp);
//printf("\n");
*/
}




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
	cout<<"input file name:"<<argv[1]<<endl; 
	int iseed = atoi(argv[2]);
	input.open(argv[1],ios::in);
     //input.open("tai25a.txt",ios::in);
		         if(!input.is_open())
		         {
		                      cout<<"error opening file";

		          }
        //reading the size,seed from the file
		          input>>size>>seed;
		          cout<<"array size:"<<size<<endl;
			cout<<"seed value:"<<argv[2]<<endl;
	//copying the size into nsize,seed into iseed variable
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

	cout<<"flatten distance:";
	for(int i=0;i<nsize;i++)
	{
		for(int j=0;j<nsize;j++)
		{
			cout<<h_dist[i * nsize+ j]<<" ";
		}
	cout<<endl;
	}
	cout<<endl;
	cout<<"flatten flows:";
	for(int i=0;i<nsize;i++)
	{
		for(int j=0;j<nsize;j++)
		{
			cout<<h_flows[i * nsize+j]<<" ";
		}
	cout<<endl;
	}
	cout<<endl;
	cout<<"size of the array:"<<nsize;						
	cout<<endl;				
				//srand(time(NULL));
				srand(iseed);
	 			int *init_sol,j,row=6144;
                init_sol = (int *)malloc(nsize * sizeof(int));
                for(int i=0;i<nsize;i++)
                init_sol[i] = i+1;
                int size_B = (row) * nsize;
				int mem_size_B = sizeof(int) * size_B;
				int *h_sol = (int *)malloc(mem_size_B);
				for(int k=0;k<row;k++)
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
				
				/*cout<<"initial solution array cpu:";
                			for(int i=0;i<row;i++)
					{
						for(int j=0;j<nsize;j++)
							cout<<h_sol[i *nsize + j]<<" ";
                                		cout<<endl;
					}
*/
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
/*cout<<"swapping locations of the array one:";				
for(int i=0;i<(nsize*(nsize-1))/2;i++)
{
	for(int j=0;j<2;j++)
	{
		cout<<h_pos[i*2+j];
	
	}
	cout<<endl;
}
cout<<endl;
int xj = ((nsize*(nsize-1))/2)/(nsize-1);
int x=0;
cout<<"swapping locations of the array two:";				
for(int i=0;i<(nsize-1);i++)
{int k=x;
for(int j=0;j<xj;j++)
{
	cout<<h_pos[(i*2)+k];
	cout<<h_pos[(i*2)+(k+1)];
	cout<<endl;
	k = k+2;
}
cout<<endl;
x=k-2;
cout<<"x value:"<<x;
cout<<endl;
}
*/
int *h_newarray;
h_newarray = (int *)malloc(row * nsize * sizeof(int));

/*cout<<"diversified array:";
for(int i=0;i<(row);i++)
{
        for(j=0;j<nsize;j++)
        {
                cout<<h_newarray[i*nsize+j]<<" ";
        }
cout<<endl;
}
*/
int *h_result,*h_frequency;			
int *d_result,*d_frequency;
int *h_divresult;			
int *d_divresult;
int *h_bestsofar,*h_bestcostsofar;
h_result = (int *)malloc(row * sizeof(int));
h_frequency = (int *)malloc(row * sizeof(int));
h_divresult = (int *)malloc((row*(nsize*(nsize-1))/2) * sizeof(int));
h_bestsofar = (int *)malloc(row * nsize * sizeof(int));
h_bestcostsofar = (int *)malloc((row) * sizeof(int));
for(int i=0;i<row;i++)
h_frequency[i] = (nsize * (nsize-1))/2;
cudaEventRecord(start,0);
int *d_bestsofar=NULL,*d_bestcostsofar=NULL,*d_newarray=NULL;
gpuErrchk( cudaMalloc((void **)&d_bestcostsofar,row * sizeof(int)) );
gpuErrchk( cudaMalloc((void **)&d_bestsofar,row * nsize * sizeof(int)) );
gpuErrchk( cudaMalloc((void **)&d_result,row * sizeof(int)) );
gpuErrchk( cudaMalloc((void **)&d_frequency,row * sizeof(int)) );
gpuErrchk( cudaMalloc((void **)&d_newarray,(row) * nsize * sizeof(int)));
gpuErrchk( cudaMalloc((void **)&d_divresult,(row*(nsize*(nsize-1))/2) * sizeof(int)) );
// declaring device array and allocating memory on gpu
int *d_dist = NULL,*d_flows = NULL ,*d_sol = NULL,*d_pos=NULL;
gpuErrchk( cudaMalloc((void **)&d_dist,nsize*nsize*sizeof(int)) );
gpuErrchk( cudaMalloc((void **)&d_flows,nsize*nsize*sizeof(int)) );
gpuErrchk( cudaMalloc((void **)&d_newarray,(row) * nsize * sizeof(int)));
gpuErrchk( cudaMalloc((void **)&d_pos,(nsize*(nsize-1))/2 * 2 * sizeof(int)));
gpuErrchk( cudaMalloc((void **)&d_sol,row*nsize*sizeof(int)) );
//copying arrays from host to device 	
gpuErrchk( cudaMemcpy(d_dist,h_dist,nsize*nsize*sizeof(int),cudaMemcpyHostToDevice) );
gpuErrchk( cudaMemcpy(d_flows,h_flows,nsize*nsize*sizeof(int),cudaMemcpyHostToDevice) );
gpuErrchk( cudaMemcpy(d_sol,h_sol,row*nsize*sizeof(int),cudaMemcpyHostToDevice) );
gpuErrchk( cudaMemcpy(d_pos,h_pos,(nsize*(nsize-1))/2 * 2 * sizeof(int),cudaMemcpyHostToDevice) );
gpuErrchk( cudaMemcpy(d_bestcostsofar,h_bestcostsofar,row * sizeof(int),cudaMemcpyHostToDevice) );
gpuErrchk( cudaMemcpy(d_bestsofar,h_bestsofar,row * nsize* sizeof(int), cudaMemcpyHostToDevice));
gpuErrchk( cudaMemcpy(d_newarray,h_newarray, row* nsize* sizeof(int), cudaMemcpyHostToDevice));
gpuErrchk( cudaMemcpy(d_result,h_result,row * sizeof(int), cudaMemcpyHostToDevice) );
gpuErrchk( cudaMemcpy(d_frequency,h_frequency,row * sizeof(int), cudaMemcpyHostToDevice) );
gpuErrchk( cudaMemcpy(d_divresult,h_divresult,(row*(nsize*(nsize-1))/2) * sizeof(int), cudaMemcpyHostToDevice) );
//cuda kernel call 
int threadsPerBlock=256;
int blockPerGrid = (row + threadsPerBlock - 1) / threadsPerBlock;
//int smSize= threadsPerBlock *nsize*nsize*sizeof(int);
cout<<"number of initial solutions:"<<row<<endl;
cout<<"number of blocks:"<<blockPerGrid<<" "<<endl;
cout<<"number of threads:"<<threadsPerBlock<<" "<<endl;
max<<<blockPerGrid,threadsPerBlock>>>(d_dist,d_flows,d_sol,nsize,row,d_result,d_bestsofar,d_bestcostsofar,d_pos,d_newarray,d_divresult,d_frequency);
gpuErrchk( cudaPeekAtLastError() );
if (cudaSuccess != cudaGetLastError()) {
return 1;
}
// wait for parent to complete
if (cudaSuccess != cudaDeviceSynchronize()) {
return 2;
}
gpuErrchk( cudaMemcpy(h_result,d_result,row * sizeof(int),cudaMemcpyDeviceToHost) );
gpuErrchk( cudaMemcpy(h_divresult,d_divresult,(row*(nsize*(nsize-1))/2) * sizeof(int),cudaMemcpyDeviceToHost) );
//gpuErrchk( cudaMemcpy(h_newarray, d_newarray , (row) * nsize* sizeof(int),cudaMemcpyDeviceToHost ));
gpuErrchk( cudaMemcpy(h_bestcostsofar,d_bestcostsofar,row * sizeof(int),cudaMemcpyDeviceToHost) );
gpuErrchk( cudaMemcpy(h_bestsofar, d_bestsofar , row * nsize* sizeof(int),cudaMemcpyDeviceToHost )); 
//gpuErrchk( cudaMemcpy(temp_sol, d_tmpsol , (row*(nsize-1)) * nsize* sizeof(int),cudaMemcpyDeviceToHost ));
gpuErrchk( cudaMemcpy(h_sol,d_sol , (row) * nsize* sizeof(int),cudaMemcpyDeviceToHost ));
cudaEventRecord(stop,0);
cudaEventSynchronize(start);
cudaEventSynchronize(stop);
/*
cout<<"last sols and their cost:";
for(int i=0;i<row;i++)
{
	for(int j=0;j<nsize;j++)
		cout<<h_sol[i *nsize + j]<<" ";
    cout<<h_result[i]<<" ";
	cout<<endl;
}
cout<<"cost of best solutions sofar and best solutin array's:"<<endl;
for(int i=0;i<row;i++)
{
        for(j=0;j<nsize;j++)
        {
                cout<<h_bestsofar[i*nsize+j]<<" ";
        }
cout<<h_bestcostsofar[i]<<" ";
cout<<endl;
}

cout<<"cost of pair wise swaps"<<endl;
for(int i=0;i<(row*(nsize*(nsize-1))/2);i++)
	cout<<h_divresult[i]<<" ";
cout<<endl;
*/
/*
int offset=5,pos,istart;
for(int i=0;i<row;i++)
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
cout<<"cpu newarray sol:";
cout<<endl;
for(int i=0;i<(row);i++)
{
        for(j=0;j<nsize;j++)
        {
                cout<<h_newarray[i*nsize+j]<<" ";
        }
cout<<endl;
}
cout<<"cost of gpu diversified initial sols:";
for(int i=0;i<row;i++)
{
	for(j=0;j<nsize;j++)
    {
                cout<<h_newarray[i*nsize+j]<<" ";
    }
cout<<h_result[i]<<" ";
cout<<endl;
}
cout<<endl;
*/
cout<<"best solutionsofar array:"<<endl;
int temp=0,lrow=0;
temp = h_bestcostsofar[0];
for(int i=1;i<row;i++)
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
cout<<"best cost:"<<temp<<endl;
//*/
cpu_endTime = clock();
cpu_ElapseTime= ((cpu_endTime - cpu_startTime) /(double) CLOCKS_PER_SEC);
cout<<"total execution time in seconds:"<<cpu_ElapseTime<<endl;
cudaEventElapsedTime(&ctime, start , stop);
cout<<"time for the kernel in milliseconds:"<<ctime<<endl;
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
