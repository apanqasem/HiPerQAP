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

