#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <assert.h>
#include <cmath>
#include "vector_functions.h"
#include "gqd_type.h"
#include "gqd_qd_utilT.h"
#include "mgs_kernelsT.h"
#include "separate.h"
#include "DefineType.h"

using namespace std;

int main ( int argc, char *argv[] )
{
   // Initialization of the execution parameters

   int BS,dim,r,mode;
   if(parse_arguments(argc,argv,&BS,&dim,&r,&mode) == 1) return 1;

   int k = dim+1;
   /* k is the number of columns of the matrix for which QR is performed 
      the last column of the matrix is the right side vector 
      of the linear system */ 

   /* Generating linear system with random coefficients, 
      the same system is generated for both (GPU and CPU) executions */
  
   int timevalue;
   if(mode == 2)
      timevalue = time(NULL); // no fixed seed to verify correctness
   else
      timevalue = 1287178355; // fixed seed for timings
   srand(timevalue);

   complexH<T1>** v = new complexH<T1>*[k];
   for(int i=0; i<k; i++) v[i] = new complexH<T1>[dim];

   complex<T>* v_h = new complex<T>[k*dim];        
   // complexH<T1>** v_copy = new complexH<T1>*[k];
   for(int i=0; i<k; i++)
   {
      // v_copy[i] = new complexH<T1>[dim];
      random_point(dim,v_h+(dim*i),v[i]);
      // random_point(dim,v_h+(dim*i)*sizeof(complex<T>),v[i]);
   }

   // Allocating memory

   int dimR = k*(k+1)/2;
   complexH<T1>* Rt = new complexH<T1>[dimR];
   complex<T>* Rt_h = new complex<T>[dimR];
   complex<T>* sol_h = new complex<T>[dim];

   // print(v,k,dim,"v");
   // for(int i=0; i<k; i++)
   //   for(int j=0; j<dim; j++) v_copy[i][j] = v[i][j]; 
   // matrixTr(v,v_copy,dim,dim); // storing transpose of matrix v in v_copy

   complexH<T1>** R = new complexH<T1>*[k];               // allocating R
   for (int i=0;i<k;i++) R[i] = new complexH<T1>[k];

   complex<T>* uNorm_h;
   uNorm_h = new complex<T>[dim];

   // GPU QR Linear System Solving

   if(mode==0 || mode==2) GPU_GS(v_h,Rt_h,sol_h,dim,dimR,k,r,BS);

   // CPU QR Linear System Solving

   if(mode==1 || mode==2) for(int i=0; i<r; i++) CPU_GS(v, R,dim,k);
   // CPU_GS_wcopy(v,v_copy,R,dim,k);

   complexH<T1>* z = new complexH<T1>[k-1];
   for(int i=0; i<k-1; i++) z[i] = R[i][k-1];

   complexH<T1>* x = new complexH<T1>[dim];

   // for Table III in the paper we do not do the back substitution
   // if(mode==1 || mode==2) for (int i=0;i<r;i++) BackSubsSec(dim,R,z,x);

   if(mode == 2) // GPU vs CPU correctness verification
      cout << "printing GPU vs CPU linear solver error : " << endl
           << GPUvsCPUVectorError(sol_h,x,dim) <<  endl;

   return 0;
}
