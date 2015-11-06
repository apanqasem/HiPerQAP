/*
 * @file        : costcheck.cpp
 * @author      : Abhilash Chaparala 
 * @description : Verify QAP solution. Compare cost of discovered solution with the one contained in the original data set 
 *
 * modification : Apan Qasem 
                  05/03/15 
 *
 */

#include<iostream>
#include<fstream>
#include<cstdlib>

using namespace std;

int calculating(int *h_flows, int *h_dist, int *sol, int size, int row) {
  int calcost;
  for(int i = 0; i < row; i++) {
    calcost = 0;
    for(int j = 0; j < size-1;j++) {
      for(int k= j + 1; k<= size-1; k++) {
	calcost = calcost + (h_flows[(sol[(i * size) + j] - 1) * size + (sol[(i* size)+k]-1)]) * h_dist[j* size+k];
      }
    }
    for(int k=1; k< size; k++) {
      for(int l=0; l < k; l++) {
	calcost = calcost + h_flows[(sol[(i * size) + k]- 1) * size + (sol[(i* size)+l] - 1)] * h_dist[k *size + l];
      }
    }
  }
  return calcost;
}

int main(int argc, char *argv[]) {

  int size, num, seed, a=0, b=0;

  int arraySizeX, arraySizeY;

  if (argc != 3) {
    cout << "Usage: " << endl;
    cout << "\t./costcheck datafile solutionfile" << endl;
    exit(1);
  }
  
  string dataFileName = argv[1];
  string solFileName = argv[2];
  
  ifstream dataFile, solFile;
  
  dataFile.open(dataFileName.c_str());
  if (!dataFile) {
    cout << "Could not open data file" << dataFileName << endl; 
    exit(1);
  }

  solFile.open(solFileName.c_str());
  if (!solFile) {
    cout << "Could not open solution file" << solFileName << endl; 
    exit(1);
  }
    
  dataFile >> size >> seed;

  arraySizeX = 2 * size;
  arraySizeY = size;

  int** array;
  array = (int**) malloc(arraySizeX * sizeof(int*));
  for (int i = 0; i < arraySizeX; i++)
    array[i] = (int*) malloc(arraySizeY*sizeof(int));

  for(int row=0; row < arraySizeX; row++) {
    for(int col=0; col < arraySizeY; col++){
      array[row][col]=0;
    } 
  }

  // flatten array dist nd flows declarations
  int size_A = size * size;
  int mem_size_A = sizeof(int) * size_A;
  int *h_dist = (int *) malloc(mem_size_A);
  int *h_flows = (int *) malloc(mem_size_A);

  for(int i=0; i< size; i++) {
    for(int j=0;j<size;j++) {
      h_dist[i * size + j] = 0;
      h_flows[i* size + j] = 0;
    }
  }
  
  while(!dataFile.eof()) {
    dataFile >> num;
    if(b == size) {
      a++;
      b=0;
    }
    if(a!= (size * 2) && b != size) {
      array[a][b]=num;
      b++;
    }
  }
  dataFile.close();
  
  for(int row = 0; row < size; row++) {
    for(int col = 0; col < size; col++) {
      h_flows[row * size + col] = array[row][col];
    }
  }


  // store in dist_dup array
  int irow = 0;
  for(int row = size; row < size * 2; row++) {
    int icol = 0;
    for(int col = 0; col < size; col++) {
      h_dist[irow *size + icol] = array[row][col];
      icol++;
    }
    irow++;
  }
  

  int row=1;
  int sol[size];
  int solSize, solCost;


  if (!(solFile >> solSize)) {
    cout << "solution file formatted incorrectly" << endl;
    exit(1);
  }

  if (solSize != size) {
    cout << "solution size does not match data set size: " << size << " " << solSize << endl;
    exit(1);
  }
  
  if (!(solFile >> solCost)) {
    cout << "solution file formatted incorrectly" << endl;
    exit(1);
  }

  for (int i = 0; i < solSize; i++) {
    if (!(solFile >> sol[i])) {
      cout << "solution file formatted incorrectly" << endl;
      exit(1);
    }
    
  }


  solFile.close();
  int cost = calculating(h_flows, h_dist, sol, size, row);

  if (cost != solCost) 
    cout << "FAIL" << endl;
  else 
    cout << "PASS" << endl;

  return 0;
}
