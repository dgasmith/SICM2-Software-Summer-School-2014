//
// authors: T. Daniel Crawford (crawdad@vt.edu) & Ed Valeev (eduard@valeyev.net)
// date  : July 8, 2014
// the use of this software is permitted under the conditions GNU General
// Public License (GPL) version 2
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cassert>

#include "diag.h"
#include "mmult.h"

#define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))

struct Atom {
    int Z;
    double x, y, z;
};

void read_geometry(const char*, std::vector<Atom>&);
double** read_1e_ints(const char* filename, int nao);
double* read_2e_ints(const char* filename, int nao);

void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
double* a, int lda, double* b, int ldb, double beta, double* c, int ldc);
int C_DSYEV(char jobz, char uplo, int n, double *a, int lda, double *w,
double *work, int lwork);
void C_DAXPY(int n, double a, double *x, int incx, double *y, int incy);

int main(int argc, char *argv[]) {

  try {
    double **X, **F, **Fp, **C, **D, **D_last, *eps;
    double **evecs, *evals, **TMP;

    /*** =========================== ***/
    /*** initialize integrals, etc.  ***/
    /*** =========================== ***/

    // read geometry from xyz file
    std::vector<Atom> atoms;
    read_geometry("geom.dat", atoms);

    // count the number of electrons
    int nelectron = 0;
    for (unsigned int i = 0; i < atoms.size(); i++) nelectron += atoms[i].Z;
    int ndocc = nelectron / 2;

    /* nuclear repulsion energy */
    double enuc = 0.0;
    for (unsigned int i = 0; i < atoms.size(); i++)
      for (unsigned int j = i + 1; j < atoms.size(); j++) {
        double xij = atoms[i].x - atoms[j].x;
        double yij = atoms[i].y - atoms[j].y;
        double zij = atoms[i].z - atoms[j].z;
        double r2 = xij*xij + yij*yij + zij*zij;
        double r = sqrt(r2);
        enuc += atoms[i].Z * atoms[j].Z / r;
      }
    printf("\tNuclear repulsion energy = %20.10lf\n", enuc);

    /* Have the user input some key data */
    int do_blas;
    printf("\n1 for BLAS 0 for plain: ");
    scanf("%d", &do_blas);
    
    int nao;
    printf("\nEnter the number of AOs: ");
    scanf("%d", &nao);

    /* overlap integrals */
    double **S = read_1e_ints("s.dat", nao);
    printf("\n\tOverlap Integrals:\n");
    print_mat(S, nao, nao, stdout);

    /* kinetic-energy integrals */
    double **H = read_1e_ints("t.dat", nao);//Store T in here right away for efficiency.
    printf("\n\tKinetic-Energy Integrals:\n");
    print_mat(H, nao, nao, stdout);

    /* nuclear-attraction integrals */
    double **V = read_1e_ints("v.dat", nao);
    printf("\n\tNuclear Attraction Integrals:\n");
    print_mat(V, nao, nao, stdout);

    /* Core Hamiltonian */
    
 
     if(do_blas){
       C_DAXPY(nao*nao, 1, V[0], 1,H[0], 1);
     }
        
     else{
       for (int i = 0; i < nao; i++)
         for (int j = 0; j < nao; j++)
           H[i][j] +=  V[i][j];
     }

    printf("\n\tCore Hamiltonian:\n");
    print_mat(H, nao, nao, stdout);

    //delete_matrix(T); //No longer need this vector, store T in H right away
    delete_matrix(V);

    /* two-electron integrals */
    double *TEI = read_2e_ints("eri.dat", nao);

    /* build the symmetric orthogonalizer X = S^(-1/2) */
    evals = init_array(nao);
    
    if(do_blas){
      evecs  = read_1e_ints("s.dat", nao);    
      double** work=init_matrix(nao,nao);
      C_DSYEV('V', 'U', nao, *evecs, nao, evals, *work, nao*nao);    
      delete[] work;
    }

    else{ 
     evecs = init_matrix(nao, nao);
     diag(nao, nao, S, evals, 1, evecs, 1e-13);
    }
     
    TMP = init_matrix(nao, nao);
     
    for (int i = 0; i < nao; i++) {
      for (int j = 0; j < nao; j++) {
        S[i][j] = 0.0;
      }
      S[i][i] = 1.0 / sqrt(evals[i]);
    }     
    
    X = init_matrix(nao, nao);
    
    if(do_blas){
    //Remember that evecs needs to be transposed here since it is stored column major
      C_DGEMM('t','n',nao,nao,nao,1, *evecs,nao,*S,nao,1,*TMP,nao);
      C_DGEMM('n','n',nao,nao,nao,1, *TMP,nao,*evecs,nao,1,*X,nao);
    }
    
    else{
      mmult(evecs, 0, S, 0, TMP, nao, nao, nao);
      mmult(TMP, 0, evecs, 1, X, nao, nao, nao);
    }
    
    delete_matrix(TMP);
    delete[] evals;
    delete_matrix(evecs);
    printf("\n\tS^-1/2 Matrix:\n");
    print_mat(X, nao, nao, stdout);


    /*** =========================== ***/
    /*** build initial-guess density ***/
    /*** =========================== ***/

    F = init_matrix(nao, nao);
    for (int i = 0; i < nao; i++)
      for (int j = 0; j < nao; j++)
        F[i][j] = H[i][j]; /* core Hamiltonian guess */

    TMP = init_matrix(nao, nao);
    Fp = init_matrix(nao, nao);
    
    if(do_blas){
      C_DGEMM('n','n',nao,nao,nao,1, *X,nao,*F,nao,1,*TMP,nao);
      C_DGEMM('n','n',nao,nao,nao,1, *TMP,nao,*X,nao,1,*Fp,nao);
	}
    
    else{
      mmult(X, 0, F, 0, TMP, nao, nao, nao);
      mmult(TMP, 0, X, 0, Fp, nao, nao, nao);
    }
    
    printf("\n\tInitial F' Matrix:\n");
    print_mat(Fp, nao, nao, stdout);

    eps = init_array(nao);
    if(do_blas){
      for(int i=0; i<nao; i++){
        for(int j=0; j<nao; j++){
          TMP[i][j] = Fp[i][j];
        }
      }
	  
	  double** work=init_matrix(nao,nao);
	  C_DSYEV('V', 'U', nao, *TMP, nao, eps, *work, nao*nao);    
	  delete[] work;
	}
	
	else{
      diag(nao, nao, Fp, eps, 1, TMP, 1e-13);
    }
    
    C = init_matrix(nao, nao);
    
    if(do_blas){
      C_DGEMM('n','t',nao,nao,nao,1, *X,nao,*TMP,nao,1,*C,nao);
    }
    else{
     mmult(X, 0, TMP, 0, C, nao, nao, nao);
    }
    
    printf("\n\tInitial C Matrix:\n");
    print_mat(C, nao, nao, stdout);

    D = init_matrix(nao, nao);
    if(do_blas){
      // now do the contraction to form the density matrix
      C_DGEMM('n','t',nao,nao,ndocc,1, *C,nao,*C,nao,1,*D,nao);
    }
    
    
    else{ //no blas
      for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++)
          for (int k = 0; k < ndocc; k++)
            D[i][j] += C[i][k] * C[j][k];
    }   
    printf("\n\tInitial Density Matrix:\n");
    print_mat(D, nao, nao, stdout);

    double escf = 0.0;
    for (int i = 0; i < nao; i++)
      for (int j = 0; j < nao; j++)
        escf += D[i][j] * (H[i][j] + F[i][j]); 
        //there are other blas routines, adds and traces, that would be useful here

    int iter = 0;
    int maxiter = 1000;

    printf(
        "\n\n Iter        E(elec)              E(tot)               Delta(E)             RMS(D)\n");
    printf(" %02d %20.12f %20.12f\n", iter, escf, escf + enuc);

    D_last = init_matrix(nao, nao);

    /*** =========================== ***/
    /*** main iterative loop ***/
    /*** =========================== ***/

    double ediff;
    double rmsd;
    double escf_last = 0.0;
    double conv = 1e-12;

    do {
      iter++;

      /* Save a copy of the energy and the density */
      escf_last = escf;
      for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++)
          D_last[i][j] = D[i][j];

      /* build a new Fock matrix */
      for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++) {
          F[i][j] = H[i][j];
          for (int k = 0; k < nao; k++)
            for (int l = 0; l < nao; l++) {
              int ij = INDEX(i, j);
              int kl = INDEX(k, l);
              int ijkl = INDEX(ij, kl);
              int ik = INDEX(i, k);
              int jl = INDEX(j, l);
              int ikjl = INDEX(ik, jl);

              F[i][j] += D[k][l] * (2.0 * TEI[ijkl] - TEI[ikjl]);
            }
        }

      if (iter == 1) {
        printf("\n\tFock Matrix:\n");
        print_mat(F, nao, nao, stdout);
      }

      /* Build new guess */
      zero_matrix(TMP, nao, nao);
      if(do_blas){
        C_DGEMM('n', 'n', nao, nao, nao, 1, *X, nao, *F, nao, 1, *TMP, nao);
      }
      
      else{
        mmult(X, 0, F, 0, TMP, nao, nao, nao);
      }
      zero_matrix(Fp, nao, nao);
      
      if(do_blas){
        C_DGEMM('n', 'n', nao, nao, nao, 1,  *TMP, nao, *X, nao, 1, *Fp, nao);
      }
      else{
        mmult(TMP, 0, X, 0, Fp, nao, nao, nao);
      }
      
      zero_matrix(TMP, nao, nao);
      zero_array(eps, nao);
      
      if(do_blas){
        double** work=init_matrix(nao,nao);
        C_DSYEV('V', 'U', nao, *Fp, nao, eps, *work, nao*nao);    
        delete[] work;
      }
      else{
        diag(nao, nao, Fp, eps, 1, TMP, 1e-13);
      }
      
      zero_matrix(C, nao, nao);
      
      if(do_blas){
        //remember to transpose the eigenvectors in TMP
        C_DGEMM('n', 't', nao, nao, nao, 1, *X, nao, *Fp, nao, 1, *C, nao);
      }
      else{
        mmult(X, 0, TMP, 0, C, nao, nao, nao);
      }
      zero_matrix(D, nao, nao);
      
      /* Build new density matrix */
      if(do_blas){
		// do the contraction to form the density matrix
		C_DGEMM('n', 't', nao, nao, ndocc, 1, *C, nao, *C, nao, 1, *D, nao);
      }
      
      else{
        for (int i = 0; i < nao; i++)
          for (int j = 0; j < nao; j++)
            for (int k = 0; k < ndocc; k++)
              D[i][j] += C[i][k] * C[j][k];
      }
      
      /* Compute new SCF energy */
      escf = 0.0;
      for(int i = 0; i < nao; i++)
        for(int j = 0; j < nao; j++)
          escf += D[i][j] * (H[i][j] + F[i][j]);

      /* Compute RMSD */
      ediff = escf - escf_last;
      rmsd = 0.0;
      for(int i = 0; i < nao; i++)
        for(int j = 0; j < nao; j++)
          rmsd += (D[i][j] - D_last[i][j]) * (D[i][j] - D_last[i][j]);

      printf(" %02d %20.12f %20.12f %20.12f %20.12f\n", iter, escf, escf + enuc,
             ediff, sqrt(rmsd));

    } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));

    delete_matrix(TMP);
    delete_matrix(D_last);
    delete_matrix(D);
    delete_matrix(C);
    delete[] eps;
    delete_matrix(Fp);
    delete_matrix(F);
    delete_matrix(X);
    delete_matrix(H);
    delete_matrix(S);
    delete[] TEI;
  } // end of try block

  catch (const char* ex) {
    std::cerr << "caught exception: " << ex << std::endl;
    return 1;
  }
  catch (std::string& ex) {
    std::cerr << "caught exception: " << ex << std::endl;
    return 1;
  }
  catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << "caught unknown exception\n";
    return 1;
  }

  return 0;
}

void read_geometry(const char *filename, std::vector<Atom>& atoms)
{
  std::ifstream is(filename);
  assert(is.good());

  size_t natom;
  is >> natom;

  atoms.resize(natom);
  for(unsigned int i = 0; i < natom; i++)
    is >> atoms[i].Z >> atoms[i].x >> atoms[i].y >> atoms[i].z;

  is.close();
}

double** read_1e_ints(const char* filename, int nao) {
  std::ifstream is(filename);
  assert(is.good());

  double** result = init_matrix(nao, nao);

  int i, j;
  double val;
  while (!is.eof()) {
    is >> i >> j >> val;
    result[i - 1][j - 1] = result[j - 1][i - 1] = val;
  }
  is.close();

  return result;
}

double* read_2e_ints(const char* filename, int nao) {
  double* result = init_array((nao*(nao+1)/2)*((nao*(nao+1)/2)+1)/2);
  std::ifstream is(filename);
  assert(is.good());

  int i, j, k, l;
  double val;
  while (!is.eof()) {
    is >> i >> j >> k >> l >> val;
    long ij = INDEX(i - 1, j - 1);
    long kl = INDEX(k - 1, l - 1);
    long ijkl = INDEX(ij, kl);
    result[ijkl] = val;
  }
  is.close();

  return result;
}

extern "C" {
extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*,
double*, int*, double*, double*, int*);
extern void dsyev_(char*, char*, int*, double*, int*, double*, double*,
int*, int*);
extern void daxpy_(int*, double*, double*, int*, double*, int*);
}

void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
double* a, int lda, double* b, int ldb, double beta, double* c, int ldc)
{
    if(m == 0 || n == 0 || k == 0) return;
    dgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta,
c, &ldc);
}

int C_DSYEV(char jobz, char uplo, int n, double *a, int lda, double *w,
double *work, int lwork){
    int info;
    dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

    return info;
}

void C_DAXPY(int n, double a, double *x, int incx, double *y, int incy)
{
   daxpy_(&n, &a, x, &incx, y, &incy);
}

