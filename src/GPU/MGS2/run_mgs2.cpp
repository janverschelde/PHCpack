// This is the main program for the modified Gram-Schmidt method 
// on a matrix of random complex numbers on the unit circle.
// Subroutines, governed by the fourth parameter "mode", compute
// 1. orthonormalization of a sequence of complex vectors
//    mode = 0: GPU only; = 1: CPU only; = 2: GPU+CPU with correctness checks;
// 2. QR decomposition of a sequence of complex vectors
//    mode = 3: GPU only; = 4: CPU only; = 5: GPU+CPU with correctness checks;
// 3. the least squares solution of a complex linear system, after a QR
//    mode = 6: GPU only; = 7: CPU only; = 8: GPU+CPU with correctness checks;
// 4. complex Newton's method on the H-Chandrasekhar problem
//    mode = 9: GPU only; = 10: CPU only; = 11: GPU+CPU with comparisons;
// 5. orthonormalization of a sequence of real vectors
//    mode = 12: GPU only; = 13: CPU only; = 14: GPU+CPU with comparisons;
// 6. QR decomposition of a sequence of real vectors
//    mode = 15: GPU only; = 16: CPU only; = 17: GPU+CPU with comparisons;
// 7. the least squares solution of a real linear system, after a QR
//    mode = 18: GPU only; = 19: CPU only; = 20: GPU+CPU with comparisons.
// 8. real Newton's method on the H-Chandrasekhar problem
//    mode = 21: GPU only; = 22: CPU only; = 23: GPU+CPU with comparisons;
// 9. complex Newton's method on the H-Chadrasekhar problem
//    mode = 24: GPU only and everything: eval, diff, qrls, update;
//    mode = 25: same as 24, with output at the end of the computations.
// The modes for GPU and CPU only are for timing purposes.
// When both GPU and CPU compute, comparisons are made for correctness
// and the accuracy of the results.
// Each of the three subroutines can be compiled for the three different
// levels of precision: double, double double, and quad double.

#include <iostream>
#include <iomanip>
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

void run_complex_orthonormalization
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

void run_real_orthonormalization
 ( T1** v, T* v_h, int BS, int dim, int freq, int mode );
/*
 * DESCRPTION :
 *   Runs the orthonormalization on the matrix v and its counterpart v_h
 *   (as stored on the host for acceleration on the device).
 * 
 * ON ENTRY :
 *   v        dim-by-dim matrix of real numbers;
 *   v_h      dim*dim real numbers stored on host for device;
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h;
 *   freq     frequency for timings in lower dimensions;
 *   mode     mode of execution: 12 (GPU only), 13 (CPU only),
 *            14 for both on GPU and CPU with tests on correctness. */

void run_complex_QR_decomposition
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode );
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

void run_real_QR_decomposition
 ( T1** v, T* v_h, int BS, int dim, int freq, int mode );
/*
 * DESCRPTION :
 *   Runs the QR decomposition on the matrix v and its counterpart v_h
 *   (as stored on the host for acceleration on the device).
 * 
 * ON ENTRY :
 *   v        dim-by-dim matrix of real numbers;
 *   v_h      dim*dim real numbers stored on host for device;
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h;
 *   freq     frequency for timings in lower dimensions;
 *   mode     mode of execution: 15 (GPU only), 16 (CPU only),
 *            17 for both on GPU and CPU with tests on correctness. */

void run_complex_QR_least_squares
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

void run_real_QR_least_squares
 ( T1** v, T* v_h, int BS, int dim, int freq, int mode );
/*
 * DESCRPTION :
 *   Solves a linear system in the least squares sense with
 *   the QR decomposition on the matrix v and its counterpart v_h
 *   (as stored on the host for acceleration on the device).
 * 
 * ON ENTRY :
 *   v        dim-by-(dim+1) matrix of real numbers;
 *   v_h      dim*(dim+1) real numbers stored on host for device;
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h,
 *            the right hand size vector is stored in column dim+1;
 *   freq     frequency for timings in lower dimensions;
 *   mode     mode of execution: 6 (GPU only), 7 (CPU only),
 *            8 for both on GPU and CPU with tests on correctness. */

void run_complex_Newton_on_chandra
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

void run_real_Newton_on_chandra
 ( T1** v, T* v_h, int BS, int dim, int freq, int mode );
/*
 * DESCRIPTION :
 *   Runs Newton's method on the system chandra of dimension dim,
 *   performing as many iterations as the value of freq.
 * 
 * ON ENTRY :
 *   v        dim-by-(dim+1) matrix of real numbers;
 *   v_h      dim*(dim+1) real numbers stored on host for device;
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h,
 *            the right hand size vector is stored in column dim+1;
 *   freq     number of Newton iterations,
 *   mode     mode of execution: 9 (GPU only), 10 (CPU only),
 *            11 for both on GPU and CPU with tests on correctness. */

void run_GPU_newton_on_chandra
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode );
/*
 * DESCRIPTION :
 *   Runs Newton's method on the discretization of the H-Chandrasekhar
 *   equation, entirely on the GPU.
 * 
 * ON ENTRY :
 *   v        dim-by-(dim+1) matrix of real numbers;
 *   v_h      dim*(dim+1) real numbers stored on host for device;
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h,
 *            the right hand size vector is stored in column dim+1;
 *   freq     number of Newton iterations,
 *   mode     mode of execution: 24 without output,
 *            25 with output at the end. */

int complex_main ( int BS, int dim, int freq, int mode, int seed );
/*
 * DESCRIPTION :
 *   This is the main program for complex data for mode < 12.
 *
 * ON ENTRY :
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h,
 *            the right hand size vector is stored in column dim+1;
 *   freq     number of Newton iterations,
 *   mode     mode of execution: 9 (GPU only), 10 (CPU only),
 *            11 for both on GPU and CPU with tests on correctness;
 *   seed     seed used to randomize sequence of random numbers. */

int real_main ( int BS, int dim, int freq, int mode, int seed );
/*
 * DESCRIPTION :
 *   This is the main program for real data for mode >= 12.
 *
 * ON ENTRY :
 *   BS       number of threads per block (block size) for GPU execution;
 *   dim      dimension of the matrices v and v_h,
 *            the right hand size vector is stored in column dim+1;
 *   freq     number of Newton iterations,
 *   mode     mode of execution: 9 (GPU only), 10 (CPU only),
 *            11 for both on GPU and CPU with tests on correctness.
 *   seed     seed used to randomize sequence of random numbers. */

int main ( int argc, char *argv[] )
{
   int BS,dim,freq,mode;  // initialization of the execution parameters
   if(parse_arguments(argc,argv,&BS,&dim,&freq,&mode) == 1) return 1;

   int timevalue;
   if(mode == 2 || mode == 5 || mode == 8)
      timevalue = time(NULL); // no fixed seed to verify correctness
      // timevalue = 1389380254;
   else
      timevalue = 1287178355; // fixed seed for timings
   srand(timevalue);

   if((mode < 12) || (mode > 23))
      return complex_main(BS,dim,freq,mode,timevalue);
   else
      return real_main(BS,dim,freq,mode,timevalue);
}

int complex_main ( int BS, int dim, int freq, int mode, int seed )
{
   complexH<T1>** v; // the sequence of vectors on the host is stored in v
   complex<T>* v_h;  // v_h is stored on the host as one long array 
   if(mode < 6)
   {
      v = new complexH<T1>*[dim];
      for(int i=0; i<dim; i++) v[i] = new complexH<T1>[dim];
      v_h = new complex<T>[dim*dim];
   }
   else              // allocate extra right hand side column
   {
      v = new complexH<T1>*[dim+1];
      for(int i=0; i<dim+1; i++) v[i] = new complexH<T1>[dim];
      v_h = new complex<T>[dim*(dim+1)];
   }
   if(mode < 9)      // generate random complex numbers on unit circle
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            double temp = random_double()*2*M_PI;
            v[i][j].init(cos(temp),sin(temp));
            v_h[i*dim+j].initH(cos(temp),sin(temp));
         }
      if(mode > 5)   // random numbers for extra right hand side column
      {
         for(int j=0; j<dim; j++)
         {
            double temp = random_double()*2*M_PI;
            v[dim][j].init(cos(temp),sin(temp));
            v_h[dim*dim+j].initH(cos(temp),sin(temp));
         }
      }
   }
   if(mode < 3)
      run_complex_orthonormalization(v,v_h,BS,dim,freq,mode);
   else if(mode < 6)
      run_complex_QR_decomposition(v,v_h,BS,dim,freq,mode);
   else if(mode < 9)
      run_complex_QR_least_squares(v,v_h,BS,dim,freq,mode);
   else if(mode < 12)
      run_complex_Newton_on_chandra(v,v_h,BS,dim,freq,mode);
   else
      run_GPU_newton_on_chandra(v,v_h,BS,dim,freq,mode);

   if(mode == 2 || mode == 5 || mode == 8)
      cout << "seed for random numbers : " << seed << endl;

   return 0;
}

int real_main ( int BS, int dim, int freq, int mode, int seed )
{
   T1** v; // sequence of vectors on the host is stored in v
   T* v_h; // v_h is stored on the host as one long array

   if(mode < 18)
   {
      v = new T1*[dim];
      for(int i=0; i<dim; i++) v[i] = new T1[dim];
      v_h = new T[dim*dim];
   }
   else               // allocate extra right hand side column
   {
      v = new T1*[dim+1];
      for(int i=0; i<dim+1; i++) v[i] = new T1[dim];
      v_h = new T[dim*(dim+1)];
   }
   if(mode < 21)      // generate random numbers in [-1, +1]
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            v[i][j] = 2.0*random_double() - 1.0; // [0, 1] becomes [-1, +1]
            qd2gqd(&v[i][j],&v_h[i*dim+j]);
         }
      if(mode > 18)   // random numbers in [-1, +1] for right hand side
      {
         for(int j=0; j<dim; j++)
         {
            v[dim][j] = 2.0*random_double() - 1.0;
            qd2gqd(&v[dim][j],&v_h[dim*dim+j]);
         }
      }
   }

   // print_real_matrix(v,dim,dim);

   if(mode < 15)
      run_real_orthonormalization(v,v_h,BS,dim,freq,mode);
   else if(mode < 18)
      run_real_QR_decomposition(v,v_h,BS,dim,freq,mode);
   else if(mode < 21)
      run_real_QR_least_squares(v,v_h,BS,dim,freq,mode);
   else
      run_real_Newton_on_chandra(v,v_h,BS,dim,freq,mode);
 
   if(mode == 14 || mode == 17 || mode == 20)
      cout << "seed for random numbers : " << seed << endl;

   return 0;
}

void run_complex_orthonormalization
 ( complexH<T1>** v, complex<T>* v_h, int BS, int dim, int freq, int mode )
{
   if(mode == 0 || mode == 2) GPU_complex_mgs2(v_h,dim,dim,freq,BS);
   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<freq; i++) CPU_complex_mgs2(v,dim,dim);
      // if(mode == 2) print_complex_matrices(v,v_h,dim,dim); // for small dim
      if(mode == 2)
      {
         complex_print_difference(v,v_h,dim,dim);
         complex_checkGPUnormal(v_h,dim,dim);
         complex_checkCPUnormal(v,dim,dim);
         complex_checkGPUorthogonal(v_h,dim,dim);
         complex_checkCPUorthogonal(v,dim,dim);
      }
   }
}

void run_real_orthonormalization
 ( T1** v, T* v_h, int BS, int dim, int freq, int mode )
{
   if(mode == 12 || mode == 14) GPU_real_mgs2(v_h,dim,dim,freq,BS);
   if(mode == 13 || mode == 14)
   {
      for(int i=0; i<freq; i++) CPU_real_mgs2(v,dim,dim);
      if(mode == 14)
      {
         real_print_difference(v,v_h,dim,dim);
         real_checkGPUnormal(v_h,dim,dim);
         real_checkCPUnormal(v,dim,dim);
         real_checkGPUorthogonal(v_h,dim,dim);
         real_checkCPUorthogonal(v,dim,dim);
      }
   }
}

void run_complex_QR_decomposition
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
   copy_complex_matrices(v,w,dim,dim);

   if(mode == 3 || mode == 5) GPU_complex_mgs2qr(v_h,R_h,dimR,dim,dim,freq,BS);
   if(mode == 4 || mode == 5)
   {
      for(int i=0; i<freq; i++) CPU_complex_mgs2qr(v,R,dim,dim);
      if(mode == 5)
      {
         complex_checkCPUnormal(v,dim,dim);
         complex_checkGPUnormal(v_h,dim,dim);
         complex_checkCPUorthogonal(v,dim,dim);
         complex_checkGPUorthogonal(v_h,dim,dim);
         complex_checkCPUdecomposition(w,v,R,dim,dim);
         complex_checkGPUdecomposition(w,v_h,R_h,dimR,dim,dim);
      }
   }
}

void run_real_QR_decomposition
 ( T1** v, T* v_h, int BS, int dim, int freq, int mode )
{
   int dimR = dim*(dim+1)/2;
   T* R_h = new T[dimR];
   T1** R = new T1*[dim];
   for(int i=0; i<dim; i++)
   {
      R[i] = new T1[dim];
      for(int j=0; j<dim; j++) R[i][j] = 0.0;
   }
   T1** w = new T1*[dim];
   for(int i=0; i<dim; i++) w[i] = new T1[dim];
   copy_real_matrices(v,w,dim,dim);

   if(mode == 15 || mode == 17) GPU_real_mgs2qr(v_h,R_h,dimR,dim,dim,freq,BS);
   if(mode == 16 || mode == 17)
   {
      for(int i=0; i<freq; i++) CPU_real_mgs2qr(v,R,dim,dim);
      if(mode == 17)
      {
         real_checkCPUnormal(v,dim,dim);
         real_checkGPUnormal(v_h,dim,dim);
         real_checkCPUorthogonal(v,dim,dim);
         real_checkGPUorthogonal(v_h,dim,dim);
         real_checkCPUdecomposition(w,v,R,dim,dim);
         real_checkGPUdecomposition(w,v_h,R_h,dimR,dim,dim);
      }
   }
}

void run_complex_QR_least_squares
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
   copy_complex_matrices(v,w,dim,dim+1);
   complex<T>* sol_h = new complex<T>[dim];
   complexH<T1>* sol = new complexH<T1>[dim];

   if(mode == 6 || mode == 8)
      GPU_complex_mgs2qrls(v_h,R_h,sol_h,dimR,dim,dim+1,freq,BS);

   if(mode == 7 || mode == 8)
   {
      for(int i=0; i<freq; i++) CPU_complex_mgs2qrls(v,R,sol,dim,dim+1);
      if(mode == 8)
      {
         //print_complex matrices(v,v_h,dim,dim+1); // is column dim+1 zero ?
         complex_checkCPUnormal(v,dim,dim);
         complex_checkGPUnormal(v_h,dim,dim);
         complex_checkCPUorthogonal(v,dim,dim);
         complex_checkGPUorthogonal(v_h,dim,dim);
         complex_checkCPUdecomposition(w,v,R,dim,dim+1);
         complex_checkGPUdecomposition(w,v_h,R_h,dimR,dim,dim+1);
         //cout << "the matrix R : " << endl; print_matrix(R,dim+1,dim+1);
         //cout << "the matrix R_h : " << endl; print_complex_vector(R_h,dimR);
         complex_checkCPUsolution(w,sol,dim,dim+1);
         //cout << "the CPU solution sol : " << endl;
         //print_complex_vector(sol,dim);
         complex_checkGPUsolution(w,sol_h,dim,dim+1);
         //cout << "the GPU solution sol_h : " << endl;
         //print_complex_vector(sol_h,dim);
      }
   }
}

void run_real_QR_least_squares
 ( T1** v, T* v_h, int BS, int dim, int freq, int mode )
{
   int dimR = (dim+1)*(dim+2)/2;      // one extra column for R
   T* R_h = new T[dimR];
   T1** R = new T1*[dim+1];
   for(int i=0; i<dim+1; i++)
   {
      R[i] = new T1[dim+1];
      for(int j=0; j<dim+1; j++) R[i][j] = 0.0;
   }
   T1** w = new T1*[dim+1];
   for(int i=0; i<dim+1; i++) w[i] = new T1[dim];
   copy_real_matrices(v,w,dim,dim+1);
   T* sol_h = new T[dim];
   T1* sol = new T1[dim];

   if(mode == 18 || mode == 20)
      GPU_real_mgs2qrls(v_h,R_h,sol_h,dimR,dim,dim+1,freq,BS);

   if(mode == 19 || mode == 20)
   {
      for(int i=0; i<freq; i++) CPU_real_mgs2qrls(v,R,sol,dim,dim+1);
      if(mode == 20)
      {
         //print_real_matrices(v,v_h,dim,dim+1); // is column dim+1 zero ?
         real_checkCPUnormal(v,dim,dim);
         real_checkGPUnormal(v_h,dim,dim);
         real_checkCPUorthogonal(v,dim,dim);
         real_checkGPUorthogonal(v_h,dim,dim);
         real_checkCPUdecomposition(w,v,R,dim,dim+1);
         real_checkGPUdecomposition(w,v_h,R_h,dimR,dim,dim+1);
         //cout << "the matrix R : " << endl; print_real matrix(R,dim+1,dim+1);
         //cout << "the matrix R_h : " << endl; print_real_vector(R_h,dimR);
         real_checkCPUsolution(w,sol,dim,dim+1);
         //cout << "the CPU solution sol : " << endl;
         //print_real_vector(sol,dim);
         real_checkGPUsolution(w,sol_h,dim,dim+1);
         //cout << "the GPU solution sol_h : " << endl;
         //print_real_vector(sol_h,dim);
      }
   }
}

void run_complex_Newton_on_chandra
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
   complex<T>* R_h = new complex<T>[dimR];

   for(int i=0; i<dim+1; i++)
   {
      R[i] = new complexH<T1>[dim+1];
      for(int j=0; j<dim+1; j++) R[i][j].init(0.0,0.0);
   }

   for(int i=0; i<dim; i++) x[i].init(1.0,0.0);

   for(int i=0; i<freq; i++)
   {
      complex_chandra_evaluate_and_differentiate(c,dim,x,v);
      cout << "      f[0] = " << v[dim][0];
      if(mode == 9 || mode == 11)
      {
         int index=0;
         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) comp1_qd2gqd(&v[k][L],&v_h[index++]);
         for(int L=0; L<dim; L++) comp1_qd2gqd(&v[dim][L],&v_h[index++]);
         GPU_complex_mgs2qrls(v_h,R_h,sol_h,dimR,dim,dim+1,1,BS);
         for(int k=0; k<dim; k++) comp1_gqd2qd(&sol_h[k],&sol[k]);
      }
      if(mode == 10 || mode == 11) CPU_complex_mgs2qrls(v,R,sol,dim,dim+1);
      for(int k=0; k<dim; k++) x[k] = x[k] + sol[k];
      cout << "delta x[0] = " << sol[0];
      cout << "      x[0] = " << x[0];
   }
   if(mode == 11)
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] = " << x[i];
}

void run_real_Newton_on_chandra
 ( T1** v, T* v_h, int BS, int dim, int freq, int mode )
{
   T1 c_a = 33.0;
   T1 c_b = 64.0;
   T1 c = c_a/c_b; // avoid representation errors, unlike 0.51234
   T1* x = new T1[dim];
   T1* sol = new T1[dim];
   T* sol_h = new T[dim];
   T1** R = new T1*[dim+1];
   int dimR = (dim+1)*(dim+2)/2;
   T* R_h = new T[dimR];

   for(int i=0; i<dim+1; i++)
   {
      R[i] = new T1[dim+1];
      for(int j=0; j<dim+1; j++) R[i][j] = 0.0;
   }

   for(int i=0; i<dim; i++) x[i] = 1.0;

   for(int i=0; i<freq; i++)
   {
      real_chandra_evaluate_and_differentiate(c,dim,x,v);
      cout << "      f[0] = " << v[dim][0] << endl;
      if(mode == 21 || mode == 23)
      {
         int index=0;
         for(int k=0; k<dim; k++)
            for(int L=0; L<dim; L++) qd2gqd(&v[k][L],&v_h[index++]);
         for(int L=0; L<dim; L++) qd2gqd(&v[dim][L],&v_h[index++]);
         GPU_real_mgs2qrls(v_h,R_h,sol_h,dimR,dim,dim+1,1,BS);
         for(int k=0; k<dim; k++) gqd2qd(&sol_h[k],&sol[k]);
      }
      if(mode == 22 || mode == 23) CPU_real_mgs2qrls(v,R,sol,dim,dim+1);
      for(int k=0; k<dim; k++) x[k] = x[k] + sol[k];
      cout << "delta x[0] = " << sol[0] << endl;
      cout << "      x[0] = " << x[0] << endl;
   }
   if(mode == 23)
      for(int i=0; i<dim; i++)
         cout << scientific << setprecision(64)
              << "x[" << i << "] = " << x[i] << endl;
}

void run_GPU_newton_on_chandra
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

   for(int i=0; i<dim; i++) x[i].init(1.0, 0.0);

   for(int L=0; L<dim; L++) comp1_qd2gqd(&x[L],&sol_h[L]);

   GPU_complex_Newton_chandra(v_h,R_h,sol_h,dimR,dim,dim+1,freq,BS);

   for(int k=0; k<dim; k++) comp1_gqd2qd(&sol_h[k],&sol[k]);

   if(mode == 25)
      for(int i=0; i<dim; i++)
         cout << scientific << setprecision(64)
              << "x[" << i << "] = " << sol[i];
}
