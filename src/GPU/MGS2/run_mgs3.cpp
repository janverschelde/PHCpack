// This is the main program for the modified Gram-Schmidt method 
// on a matrix of random complex numbers on the unit circle.
// There are four subroutines, to compute
// 1. an orthonormalization of a sequence of vectors;
// 2. a QR decomposition of a sequence of vector; and
// 3. the least squares solution of a linear system,
//    after a QR decomposition;
// 4. to run Newton's method on the chandra problem.
// For each of the three subroutines, there are three different modes:
// 0, 3, 6, 9: GPU only; 1, 4, 7, 10: CPU only; 2, 5, 8, 11: GPU+CPU with checks
// on the correctness and accuracy.
// Each of the three subroutines can be compiled for the three different
// levels of precision: double, double double, and quad double.

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "vector_functions.h"
#include "gqd_type.h"
#include "gqd_qd_util.h"
#include "gqd_qd_utilT.h"
#include "DefineType.h"
#include "mgs2_host.h"
#include "mgs2_kernels.h"
#include "chandra.h"

using namespace std;

void run_orthonormalization
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode );
/*
 * DESCRPTION :
 *   Runs the orthonormalization on the matrix v and its counterpart v_h
 *   (as stored on the host for acceleration on the device).
 * 
 * ON ENTRY :
 *   v        dim-by-dim matrix of complex numbers;
 *   v_h      dim*dim complex numbers stored on host for device;
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h;
 *   freq     frequency for timings in lower dimensions;
 *   mode     mode of execution: 0 (GPU only), 1 (CPU only),
 *            2 for both on GPU and CPU with tests on correctness. */

/*void run_QR_decomposition
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode );*/
/*
 * DESCRPTION :
 *   Runs the QR decomposition on the matrix v and its counterpart v_h
 *   (as stored on the host for acceleration on the device).
 * 
 * ON ENTRY :
 *   v        dim-by-dim matrix of complex numbers;
 *   v_h      dim*dim complex numbers stored on host for device;
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h;
 *   freq     frequency for timings in lower dimensions;
 *   mode     mode of execution: 3 (GPU only), 4 (CPU only),
 *            5 for both on GPU and CPU with tests on correctness. */

void run_QR_least_squares
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode );
/*
 * DESCRPTION :
 *   Solves a linear system in the least squares sense with
 *   the QR decomposition on the matrix v and its counterpart v_h
 *   (as stored on the host for acceleration on the device).
 * 
 * ON ENTRY :
 *   v        dim-by-(dim+1) matrix of complex numbers;
 *   v_h      dim*(dim+1) complex numbers stored on host for device;
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h,
 *            the right hand size vector is stored in column dim+1;
 *   freq     frequency for timings in lower dimensions;
 *   mode     mode of execution: 6 (GPU only), 7 (CPU only),
 *            8 for both on GPU and CPU with tests on correctness. */

void run_newton_on_chandra
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode );
/*
 * DESCRIPTION :
 *   Runs Newton's method on the system chandra of dimension dim,
 *   performing as many iterations as the value of freq.
 * 
 * ON ENTRY :
 *   v        dim-by-(dim+1) matrix of complex numbers;
 *   v_h      dim*(dim+1) complex numbers stored on host for device;
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h,
 *            the right hand size vector is stored in column dim+1;
 *   freq     number of Newton iterations,
 *   mode     mode of execution: 9 (GPU only), 10 (CPU only),
 *            11 for both on GPU and CPU with tests on correctness. */

int main ( int argc, char *argv[] )
{
   int BS,dim,freq,mode;  // initialization of the execution parameters
   if(parse_arguments(argc,argv,&BS,&dim,&freq,&mode) == 1) return 1;


   complexH<T1>** v; // the sequence of vectors on the host is stored in v
   complex<T>* v_h;  // v_h is stored on the host as one long array 
   v = new complexH<T1>*[dim+1];
   for(int i=0; i<dim+1; i++) v[i] = new complexH<T1>[dim];
   v_h = new complex<T>[dim*(dim+1)];

   int timevalue;
   timevalue = 1287178355; // fixed seed for timings
   srand(timevalue);
   run_newton_on_chandra(v,v_h,BS,dim,freq,mode);

   return 0;
}

void run_orthonormalization
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode )
{
   if(mode == 0 || mode == 2) GPU_mgs2(v_h,dim,dim,freq,BS);
   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<freq; i++) CPU_mgs2(v,dim,dim);
      // if(mode == 2) print_matrices(v,v_h,dim,dim); // only for small dim
      if(mode == 2)
      {
         print_difference(v,v_h,dim,dim);
         checkGPUnormal(v_h,dim,dim);
         checkCPUnormal(v,dim,dim);
         checkGPUorthogonal(v_h,dim,dim);
         checkCPUorthogonal(v,dim,dim);
      }
   }
}

/*void run_QR_decomposition
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode )
{
   int dimR = dim*(dim+1)/2;
   complex<T>* R_h = new complex<T>[dimR];
   complexH<T1>** R = new complexH<T1>*[dim];
   for(int i=0; i<dim; i++)
   {
      R[i] = new complexH<T1>[dim];
      for(int j=0; j<dim; j++) R[i][j].init(0.0,0.0);
   }
   complexH<T1>** w = new complexH<T1>*[dim];
   for(int i=0; i<dim; i++) w[i] = new complexH<T1>[dim];
   copy_matrices(v,w,dim,dim);

   if(mode == 3 || mode == 5) GPU_mgs2qr(v_h,R_h,dimR,dim,dim,freq,BS);
   if(mode == 4 || mode == 5)
   {
      for(int i=0; i<freq; i++) CPU_mgs2qr(v,R,dim,dim);
      if(mode == 5)
      {
         checkCPUnormal(v,dim,dim);
         checkGPUnormal(v_h,dim,dim);
         checkCPUorthogonal(v,dim,dim);
         checkGPUorthogonal(v_h,dim,dim);
         checkCPUdecomposition(w,v,R,dim,dim);
         checkGPUdecomposition(w,v_h,R_h,dimR,dim,dim);
      }
   }
}*/

void run_QR_least_squares
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode )
{
   int dimR = (dim+1)*(dim+2)/2;      // one extra column for R
   complex<T>* R_h = new complex<T>[dimR];
   complexH<T1>** R = new complexH<T1>*[dim+1];
   for(int i=0; i<dim+1; i++)
   {
      R[i] = new complexH<T1>[dim+1];
      for(int j=0; j<dim+1; j++) R[i][j].init(0.0,0.0);
   }
   complexH<T1>** w = new complexH<T1>*[dim+1];
   for(int i=0; i<dim+1; i++) w[i] = new complexH<T1>[dim];
   copy_matrices(v,w,dim,dim+1);
   complex<T>* sol_h = new complex<T>[dim];
   complexH<T1>* sol = new complexH<T1>[dim];

   if(mode == 6 || mode == 8)
      GPU_mgs2qrls(v_h,R_h,sol_h,dimR,dim,dim+1,freq,BS);

   if(mode == 7 || mode == 8)
   {
      for(int i=0; i<freq; i++) CPU_mgs2qrls(v,R,sol,dim,dim+1);
      if(mode == 8)
      {
         // print_matrices(v,v_h,dim,dim+1); // column dim+1 should be zero
         checkCPUnormal(v,dim,dim);
         checkGPUnormal(v_h,dim,dim);
         checkCPUorthogonal(v,dim,dim);
         checkGPUorthogonal(v_h,dim,dim);
         checkCPUdecomposition(w,v,R,dim,dim+1);
         checkGPUdecomposition(w,v_h,R_h,dimR,dim,dim+1);
         //cout << "the matrix R : " << endl; print_matrix(R,dim+1,dim+1);
         //cout << "the matrix R_h : " << endl; print_vector(R_h,dimR);
         checkCPUsolution(w,sol,dim,dim+1);
         //cout << "the CPU solution sol : " << endl; print_vector(sol,dim);
         checkGPUsolution(w,sol_h,dim,dim+1);
         //cout << "the GPU solution sol_h : " << endl; print_vector(sol_h,dim);
      }
   }
}

void run_newton_on_chandra
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode )
{
   complexH<T1> c_a;
   c_a.init(33.0,0.0);
   complexH<T1> c_b;
   c_b.init(64.0,0.0);
   complexH<T1> c = c_a/c_b; // avoid representation errors, unlike 0.51234
   complexH<T1>* x = new complexH<T1>[dim];
   complexH<T1>* sol = new complexH<T1>[dim];
   complex<T>* sol_h = new complex<T>[dim];
   complexH<T1>** R = new complexH<T1>*[dim+1];
   int dimR = (dim+1)*(dim+2)/2;
   complex<T>* R_h = new complex<T>[dim];

   for(int i=0; i<dim+1; i++)
   {
      R[i] = new complexH<T1>[dim+1];
      for(int j=0; j<dim+1; j++) R[i][j].init(0.0,0.0);
   }

   for(int i=0; i<dim; i++){
      //x[i].init(T1(1)/(i+2),T1(i+2)/(i+3));
      x[i].init(1.0, 0.0);
   }

   for(int L=0; L<dim; L++) comp1_qd2gqd(&x[L],&sol_h[L]);

   GPU_mgs2qrls(v_h,R_h,sol_h,dimR,dim,dim+1,freq,BS);

   for(int k=0; k<dim; k++) comp1_gqd2qd(&sol_h[k],&sol[k]);

   cout << sol[0];
}
