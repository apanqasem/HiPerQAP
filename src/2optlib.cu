#include<cstdlib>
#include<2optlib.h>

/*
 * generate initial solutions  
 */
void genInitSolutions(int **solns, int size, int num_solns) {

  // create a single sequential solution 
  // i.e., 0 -> 1, 1 -> 2, 2 -> 3 ... 
  int *seq_soln = (int *) malloc(size * sizeof(int));
  for(int i = 0; i < size; i++)
    seq_soln[i] = i + 1;

  int *loc_solns = (int *) malloc(sizeof(int) * num_solns * size);

  for(int j = 0; j < num_solns; j++) {
    for(int i = size - 1; i > 0; i--) {
      // swap a random pair 
      int randPos = rand() % (i + 1);
      int temp = seq_soln[i];
      seq_soln[i] = seq_soln[randPos];
      seq_soln[randPos] = temp;
    }
    
    // copy rand
    for(int i = 0; i < size; i++) {
      loc_solns[j * size + i] = seq_soln[i];
    }
  }
  (*solns) = loc_solns;
}

int findMin(int *values, int solns, int *index) {
  int minCost = 0, minIndex = 0;
  minCost = values[0];
  for(int i = 1; i < solns; i++)
    if (values[i] < minCost) {
      minCost = values[i];
      minIndex = i;
    }
  (*index) = minIndex;
  return minCost;  
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
  
