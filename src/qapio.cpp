#include<iostream>
#include<fstream>
#include<cstdlib>
#include<qapio.h>

using namespace std;

void printFlattenedArray(int *values, int size0, int size1, string msg) {
  cout << msg << endl;
  for(int i = 0; i < size0; i++) {
    for(int j = 0; j < size1; j++)
      cout << values[i * size1 + j] << " ";
    cout << endl;
  }
  cout << endl;
  return;
}

void printRegularArray(int *values, int size, string msg) {
  cout << msg << endl;
  for(int i = 0; i < size; i++) 
    cout << values[i] << " ";
  cout << endl;
}


void readData(string filename, int ***matrix, int *size) {

  ifstream input;
  input.open(filename.c_str(),ios::in);
  if(!input.is_open()) {
    cout <<"error opening file: " << filename << endl;
    exit(0);
  }

  int seed, probSize;
  // first line in input file specifies problem size 
  input >> probSize >> seed;
  (*size) = probSize;

  // declaring array to copy the matrix from file into array    	
  int** array;
  array = (int**) malloc(2 * probSize * sizeof(int*));
  for (int i = 0; i < 2 * probSize; i++)
    array[i] = (int*) malloc(probSize * sizeof(int));

  for(int i = 0; i < (2 * probSize); i++) {
    for(int j = 0; j < probSize; j++)
      array[i][j] = 0;
  }

  int num = 0, a = 0, b=0;	
  while(!input.eof()) {
    input >> num;
    if(b == probSize) {
      a++;
      b=0;
    }
    if(a != (probSize * 2) && b != probSize) {
      array[a][b] = num;
      b++;
    }
  }
  input.close();
  
  (*matrix) = array;
  return;
}

/* 
 * Split a matrix into two flattened single dimensional arrays
 * The arrays contain flows and distances 
 * Flatenned arrays makes it easier to allocate data to GPU threads 
 */ 
void splitAndFlatten(short **flow_out, short **dist_out, int **matrix, int size) {

  short *dist = (short *) malloc (sizeof(int) * size * size);
  short *flow = (short *) malloc (sizeof(int) * size * size);

  for(int i = 0; i < size * size; i++) {
      flow[i] = 0;
      dist[i] = 0;
  }

  for(int row = 0; row < size; row++) {
    for(int col = 0; col < size; col++)
      flow[row * size + col] = matrix[row][col];
  }

  int irow = 0;
  for(int row = size; row < size * 2; row++) {
    for(int col = 0; col < size; col++)
      dist[irow * size + col] = matrix[row][col];
    irow++;
  }

  *flow_out = flow;
  *dist_out = dist;
  return;
}
  

/* 
 * Split a matrix into two flattened single dimensional arrays
 * The arrays contain flows and distances 
 * Flatenned arrays makes it easier to allocate data to GPU threads 
 */ 
void splitAndFlattenInt(int **flow_out, int **dist_out, int **matrix, int size) {

  int *dist = (int *) malloc (sizeof(int) * size * size);
  int *flow = (int *) malloc (sizeof(int) * size * size);

  for(int i = 0; i < size * size; i++) {
      flow[i] = 0;
      dist[i] = 0;
  }

  for(int row = 0; row < size; row++) {
    for(int col = 0; col < size; col++)
      flow[row * size + col] = matrix[row][col];
  }

  int irow = 0;
  for(int row = size; row < size * 2; row++) {
    for(int col = 0; col < size; col++)
      dist[irow * size + col] = matrix[row][col];
    irow++;
  }

  *flow_out = flow;
  *dist_out = dist;
  return;
}


