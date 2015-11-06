using namespace std;

void printFlattenedArray(int *values, int size0, int size1, string msg);
void printRegularArray(int *values, int size, string msg);

void readData(string filename, int ***matrix, int *size);
void splitAndFlatten(short **flow_out, short **dist_out, int **matrix, int size);
void splitAndFlattenInt(int **flow_out, int **dist_out, int **matrix, int size);



