// Fill two square matrices with random numbers and multiply them together
#include <mpi.h>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <sys/time.h>

using namespace std;

double my_rand();
double trace(double**,int);
double **matrix_init(int row, int col);
void matrix_delete(double **A);
void symm_matrix_fill(double **A, int dim);
uint64_t time_musec();

void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
double* a, int lda, double* b, int ldb, double beta, double* c, int ldc);
 
int main(int argc, char *argv[])
{
  if(argc != 2) {
    printf("Matrix rank  required as argument.\n");
    return 1;
  }
  
   MPI_Init(&argc, &argv);
 
      int nproc, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  
  int dim = atoi(argv[1]); // matrix rank

  const char * path = argv[2]; 
	
  const char * plain = "plain";
  const char * mkl = "mkl";
  const char * omp = "omp";
 const char * mpi = "mpi";

double**Cf= matrix_init(dim, dim);


  double **A = matrix_init(dim, dim);
  double **B = matrix_init(dim, dim);
  double **C = matrix_init(dim, dim);

  uint64_t begin = time_musec();


if(rank==0){
  printf("Matrix dimension is %d.\n", dim);
   srand(time(0));
   symm_matrix_fill(A, dim);
  symm_matrix_fill(B, dim);
}

MPI_Bcast(*A,dim*dim, MPI_DOUBLE,0,MPI_COMM_WORLD); 
MPI_Bcast(*B,dim*dim, MPI_DOUBLE,0,MPI_COMM_WORLD); 


/*
if(rank==0){

  for(int i=0; i < dim; i++){
    for(int j=0; j < dim; j++) {
      //C[i][j] = 0.0;
      for(int k=0; k < dim; k++){
        C[i][j] += A[i][k] * B[j][k];
    }
    }
   }
   printf("\n Matrix Trace, just in case: %f \n", trace(C,dim));

for(int i=0;i<dim;i++){
for(int j=0;j<dim;j++){
C[i][j]=0;
}
}


}


if(strcmp(path,mkl)){
  C_DGEMM('n','t', dim, dim, dim, 1.0, A[0], dim, B[0], dim, 0.0, C[0], dim);
}

if(strcmp(path,omp)){
#pragma omp parallel for collapse(3)
 for(int i=0; i < dim; i++){
    for(int j=0; j < dim; j++) {
      //C[i][j] = 0.0;
      for(int k=0; k < dim; k++)
        C[i][j] += A[i][k] * B[j][k];
    }
  }
}

*/


int count =0;
 
 for(int i=0; i < dim; i++){
	for(int j=0; j < dim; j++) {
      //C[i][j] = 0.0;
    for(int k=0; k < dim; k++){
     if(rank==(count%nproc)){
	  C[i][j] += A[i][k] * B[j][k];

    }
     count++;
  
    }
    
    }
}

   if(MPI_Reduce(*C, *Cf, dim*dim, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD)!= MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, 3);


uint64_t end = time_musec();
if(rank==0){

  printf("Time for MM = %5.10f sec.\n", (double) (end-begin)/1000000.0);
  printf("\n Matrix Trace, just in case: %f \n", trace(Cf,dim));

}
matrix_delete(Cf);

  matrix_delete(A);
  matrix_delete(B);
  matrix_delete(C);

MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}

// Allocate a (row-wise) square matrix of dimension dim
// Rows are contiguous in memory
double **matrix_init(int row, int col)
{
  double **A = new double* [row];
  double *B = new double[row*col];
  memset(static_cast<void*>(B), 0, row*col*sizeof(double));
  for(int i=0; i < row; i++){
    A[i] = &(B[i*col]);
}
//let's fill the matrix with zeros here.
for(int i=0;i<row;i++){
for(int j=0;j<col;j++){
A[i][j]=0;
}
}

  return A;
}

// Delete a matrix A of dimension dim initialized with matrix_init()
void matrix_delete(double **A)
{
  delete[] A[0];
  delete[] A;
}

// Fill a symmetric matrix A of dimension dim with (somewhat) random numbers
void symm_matrix_fill(double **A, int dim)
{
  for(int i=0; i < dim; i++)
    for(int j=0; j <= i; j++)
      A[i][j] = A[j][i] = my_rand();
}

// Generates positive and negative small numbers
double my_rand()  // caller must seed
{
  return ((double) (rand() % RAND_MAX/16 - RAND_MAX/32))/((double) RAND_MAX);
}

// returns the number of microseconds since the start of epoch
uint64_t time_musec()
{
 struct timeval tv;
 gettimeofday(&tv, NULL);

 uint64_t ret = tv.tv_sec * 1000000 + tv.tv_usec;

 return ret;
}


double trace(double**C,int N){

double trace =0;
for(int i=0;i<N;i++){

trace+= C[i][i];
}

return trace;
}

extern "C" {
extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*,
double*, int*, double*, double*, int*);
}

void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
double* a, int lda, double* b, int ldb, double beta, double* c, int ldc)
{
    if(m == 0 || n == 0 || k == 0) return;
    dgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
}

