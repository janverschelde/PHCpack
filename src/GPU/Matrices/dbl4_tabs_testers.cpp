// The file dbl4_tabs_testers.cpp defines the functions specified in
// the file dbl4_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "quad_double_functions.h"
#include "random4_matrices.h"
#include "dbl4_factorizations.h"
#include "dbl4_factorizations.h"
#include "dbl4_tabs_host.h"
#include "dbl4_tabs_kernels.h"
#include "dbl_test_utilities.h"
#include "dbl4_test_utilities.h"
#include "write_dbl4_bstimeflops.h"
#include "dbl_data_files.h"
#include "dbl_tabs_testers.h"

using namespace std;

void test_real4_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Ahihi = new double*[dim];
   double **Alohi = new double*[dim];
   double **Ahilo = new double*[dim];
   double **Alolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Ahihi[i] = new double[dim];
      Alohi[i] = new double[dim];
      Ahilo[i] = new double[dim];
      Alolo[i] = new double[dim];
   }
   // random_dbl4_upper_matrix(dim,dim,Ahi,Alo);
   dbl4_random_upper_factor(dim,Ahihi,Alohi,Ahilo,Alolo);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihi[i][j] << "  " << Alohi[i][j] << endl
                 << "          "
                 << Ahilo[i][j] << "  " << Alolo[i][j] << endl;
   }
   double *solhihi = new double[dim];
   double *sollohi = new double[dim];
   double *solhilo = new double[dim];
   double *sollolo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solhihi[i] = 1.0; sollohi[i] = 0.0;
      solhilo[i] = 0.0; sollolo[i] = 0.0;
   }
   double *rhshihi = new double[dim];
   double *rhslohi = new double[dim];
   double *rhshilo = new double[dim];
   double *rhslolo = new double[dim];
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<dim; i++)
   {
      rhshihi[i] = 0.0; rhslohi[i] = 0.0;
      rhshilo[i] = 0.0; rhslolo[i] = 0.0;

      for(int j=0; j<dim; j++)  // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         qdf_mul(Ahihi[i][j],Alohi[i][j],Ahilo[i][j],Alolo[i][j],
               solhihi[j], sollohi[j], solhilo[j], sollolo[j],
              &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdf_inc(&rhshihi[i],&rhslohi[i],&rhshilo[i],&rhslolo[i],
                  acchihi,    acclohi,    acchilo,    acclolo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : "
              << rhshihi[i] << "  " << rhslohi[i] << endl
              << "       "
              << rhshilo[i] << "  " << rhslolo[i] << endl;
   }
   double **invAhihi_h = new double*[dim];
   double **invAlohi_h = new double*[dim];
   double **invAhilo_h = new double*[dim];
   double **invAlolo_h = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      invAhihi_h[i] = new double[dim];
      invAlohi_h[i] = new double[dim];
      invAhilo_h[i] = new double[dim];
      invAlolo_h[i] = new double[dim];
   }
   double timelapsed_h,timelapsed_d,elapsedms;

   cout << "-> CPU computes the inverse ..." << endl;

   CPU_dbl4_upper_inverse
      (dim,Ahihi,     Alohi,     Ahilo,     Alolo,
        invAhihi_h,invAlohi_h,invAhilo_h,invAlolo_h,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The CPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_h[" << i << "][" << j << "] : "
                 << invAhihi_h[i][j] << "  " << invAlohi_h[i][j] << endl
                 << "               "
                 << invAhilo_h[i][j] << "  " << invAlolo_h[i][j] << endl;
   }
   double **invAhihi_d = new double*[dim];
   double **invAlohi_d = new double*[dim];
   double **invAhilo_d = new double*[dim];
   double **invAlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      invAhihi_d[i] = new double[dim];
      invAlohi_d[i] = new double[dim];
      invAhilo_d[i] = new double[dim];
      invAlolo_d[i] = new double[dim];
   }
   cout << "-> GPU computes the inverse ..." << endl;

   GPU_dbl4_upper_inverse
      (dim,Ahihi,     Alohi,     Ahilo,     Alolo,
        invAhihi_d,invAlohi_d,invAhilo_d,invAlolo_d,
       &elapsedms,&timelapsed_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_d[" << i << "][" << j << "] : "
                 << invAhihi_d[i][j] << "  " << invAlohi_d[i][j] << endl
                 << "               "
                 << invAhilo_d[i][j] << "  " << invAlolo_d[i][j] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on inverse : "
        << dbl4_Matrix_Difference_Sum
              (dim,invAhihi_h,invAlohi_h,invAhilo_h,invAlolo_h,
                   invAhihi_d,invAlohi_d,invAhilo_d,invAlolo_d)
        << endl;

   double *xhihi = new double[dim];
   double *xlohi = new double[dim];
   double *xhilo = new double[dim];
   double *xlolo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      xhihi[i] = 0.0; xlohi[i] = 0.0;
      xhilo[i] = 0.0; xlolo[i] = 0.0;

      for(int j=0; j<dim; j++)   // x[i] = x[i] + invA_h[i][j]*rhs[j];
      {
         qdf_mul(invAhihi_h[i][j],invAlohi_h[i][j],
                 invAhilo_h[i][j],invAlolo_h[i][j],
                 rhshihi[j],rhslohi[j],rhshilo[j],rhslolo[j],
                &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdf_inc(&xhihi[i],&xlohi[i],&xhilo[i],&xlolo[i],
                acchihi,  acclohi,  acchilo,  acclolo);
      }
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The solution computed with the CPU inverse :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhihi[i] << "  " << xlohi[i] << endl
              << "       "
              << xhilo[i] << "  " << xlolo[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on solution : "
        << dbl4_Difference_Sum(dim, solhihi,sollohi,solhilo,sollolo,
                                      xhihi,  xlohi,  xhilo,  xlolo) << endl;
   cout << "Condition number : "
        << dbl4_condition(dim,Ahihi,     Alohi,      Ahilo,     Alolo,
                           invAhihi_h,invAlohi_h, invAhilo_h,invAlolo_h)
        << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
   cout << "                     Time spent by the kernel : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
}

void test_cmplx4_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Arehihi = new double*[dim];
   double **Arelohi = new double*[dim];
   double **Arehilo = new double*[dim];
   double **Arelolo = new double*[dim];
   double **Aimhihi = new double*[dim];
   double **Aimlohi = new double*[dim];
   double **Aimhilo = new double*[dim];
   double **Aimlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Arehihi[i] = new double[dim];
      Arelohi[i] = new double[dim];
      Arehilo[i] = new double[dim];
      Arelolo[i] = new double[dim];
      Aimhihi[i] = new double[dim];
      Aimlohi[i] = new double[dim];
      Aimhilo[i] = new double[dim];
      Aimlolo[i] = new double[dim];
   }
   // random_cmplx4_upper_matrix(dim,dim,Arehi,Arelo,Aimhi,Aimlo);
   cmplx4_random_upper_factor
      (dim,Arehihi,Arelohi,Arehilo,Arelolo,Aimhihi,Aimlohi,Aimhilo,Aimlolo);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihi[i][j] << "  " << Arelohi[i][j] << endl
                 << "            "
                 << Arehilo[i][j] << "  " << Arelolo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihi[i][j] << "  " << Aimlohi[i][j] << endl
                 << "            "
                 << Aimhilo[i][j] << "  " << Aimlolo[i][j] << endl;
         }
   }
   double *solrehihi = new double[dim];
   double *solrelohi = new double[dim];
   double *solrehilo = new double[dim];
   double *solrelolo = new double[dim];
   double *solimhihi = new double[dim];
   double *solimlohi = new double[dim];
   double *solimhilo = new double[dim];
   double *solimlolo = new double[dim];

   for(int i=0; i<dim; i++) 
   {
      solrehihi[i] = 1.0; solrelohi[i] = 0.0;
      solrehilo[i] = 0.0; solrelolo[i] = 0.0;
      solimhihi[i] = 0.0; solimlohi[i] = 0.0;
      solimhilo[i] = 0.0; solimlolo[i] = 0.0;
   }
   double *rhsrehihi = new double[dim];
   double *rhsrelohi = new double[dim];
   double *rhsrehilo = new double[dim];
   double *rhsrelolo = new double[dim];
   double *rhsimhihi = new double[dim];
   double *rhsimlohi = new double[dim];
   double *rhsimhilo = new double[dim];
   double *rhsimlolo = new double[dim];
   double acc1hihi,acc1lohi,acc1hilo,acc1lolo;
   double acc2hihi,acc2lohi,acc2hilo,acc2lolo;
   double acc3hihi,acc3lohi,acc3hilo,acc3lolo;
   double acc4hihi,acc4lohi,acc4hilo,acc4lolo;

   for(int i=0; i<dim; i++)
   {
      rhsrehihi[i] = 0.0; rhsrelohi[i] = 0.0;
      rhsrehilo[i] = 0.0; rhsrelolo[i] = 0.0;
      rhsimhihi[i] = 0.0; rhsimlohi[i] = 0.0;
      rhsimhilo[i] = 0.0; rhsimlolo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
               solrehihi[j], solrelohi[j], solrehilo[j], solrelolo[j],
               &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
         qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
               solimhihi[j], solimlohi[j], solimhilo[j], solimlolo[j],
               &acc2hihi,    &acc2lohi,    &acc2hilo,   &acc2lolo);
         qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
               solrehihi[j], solrelohi[j], solrehilo[j], solrelolo[j],
               &acc3hihi,    &acc3lohi,    &acc3hilo,    &acc3lolo);
         qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
               solimhihi[j], solimlohi[j], solimhilo[j], solimlolo[j],
               &acc4hihi,    &acc4lohi,    &acc4hilo,    &acc4lolo);
         qdf_inc(&rhsrehihi[i],&rhsrelohi[i],&rhsrehilo[i],&rhsrelolo[i],
                   acc1hihi,     acc1lohi,     acc1hilo,     acc1lolo);
         qdf_dec(&rhsrehihi[i],&rhsrelohi[i],&rhsrehilo[i],&rhsrelolo[i],
                   acc2hihi,     acc2lohi,     acc2hilo,     acc2lolo);
         qdf_inc(&rhsimhihi[i],&rhsimlohi[i],&rhsimhilo[i],&rhsimlolo[i],
                   acc3hihi,     acc3lohi,     acc3hilo,     acc3lolo);
         qdf_inc(&rhsimhihi[i],&rhsimlohi[i],&rhsimhilo[i],&rhsimlolo[i],
                   acc4hihi,     acc4lohi,     acc4hilo,     acc4lolo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehihi[i] << "  " << rhsrelohi[i] << endl
              << "         "
              << rhsrehilo[i] << "  " << rhsrelolo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhihi[i] << "  " << rhsimlohi[i] << endl
              << "         "
              << rhsimhilo[i] << "  " << rhsimlolo[i] << endl;
      }
   }
   double **invArehihi_h = new double*[dim];
   double **invArelohi_h = new double*[dim];
   double **invArehilo_h = new double*[dim];
   double **invArelolo_h = new double*[dim];
   double **invAimhihi_h = new double*[dim];
   double **invAimlohi_h = new double*[dim];
   double **invAimhilo_h = new double*[dim];
   double **invAimlolo_h = new double*[dim];
 
   for(int i=0; i<dim; i++)
   {
      invArehihi_h[i] = new double[dim];
      invArelohi_h[i] = new double[dim];
      invArehilo_h[i] = new double[dim];
      invArelolo_h[i] = new double[dim];
      invAimhihi_h[i] = new double[dim];
      invAimlohi_h[i] = new double[dim];
      invAimhilo_h[i] = new double[dim];
      invAimlolo_h[i] = new double[dim];
   }
   double timelapsed_h,timelapsed_d,elapsedms;

   cout << "-> CPU computes the inverse ..." << endl;

   CPU_cmplx4_upper_inverse
      (dim,   Arehihi,     Arelohi,     Arehilo,     Arelolo, 
              Aimhihi,     Aimlohi,     Aimhilo,     Aimlolo,
           invArehihi_h,invArelohi_h,invArehilo_h,invArelolo_h,
           invAimhihi_h,invAimlohi_h,invAimhilo_h,invAimlolo_h,
       &timelapsed_h);

   if(verbose > 0)
   {
      cout << "The CPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "invA_h[" << i << "][" << j << "]re : "
                 << invArehihi_h[i][j] << "  " << invArelohi_h[i][j] << endl
                 << "                 "
                 << invArehilo_h[i][j] << "  " << invArelolo_h[i][j] << endl;
            cout << "invA_h[" << i << "][" << j << "]im : "
                 << invAimhihi_h[i][j] << "  " << invAimlohi_h[i][j] << endl
                 << "                 "
                 << invAimhilo_h[i][j] << "  " << invAimlolo_h[i][j] << endl;
         }
   }
   double **invArehihi_d = new double*[dim];
   double **invArelohi_d = new double*[dim];
   double **invArehilo_d = new double*[dim];
   double **invArelolo_d = new double*[dim];
   double **invAimhihi_d = new double*[dim];
   double **invAimlohi_d = new double*[dim];
   double **invAimhilo_d = new double*[dim];
   double **invAimlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      invArehihi_d[i] = new double[dim];
      invArelohi_d[i] = new double[dim];
      invArehilo_d[i] = new double[dim];
      invArelolo_d[i] = new double[dim];
      invAimhihi_d[i] = new double[dim];
      invAimlohi_d[i] = new double[dim];
      invAimhilo_d[i] = new double[dim];
      invAimlolo_d[i] = new double[dim];
   }

   cout << "-> GPU computes the inverse ..." << endl;

   GPU_cmplx4_upper_inverse
      (dim,   Arehihi,     Arelohi,     Arehilo,     Arelolo,     
              Aimhihi,     Aimlohi,     Aimhilo,     Aimlolo,
           invArehihi_d,invArelohi_d,invArehilo_d,invArelolo_d,
           invAimhihi_d,invAimlohi_d,invAimhilo_d,invAimlolo_d,
       &elapsedms,&timelapsed_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "invA_d[" << i << "][" << j << "]re : "
                 << invArehihi_d[i][j] << "  " << invArelohi_d[i][j] << endl
                 << "                 "
                 << invArehilo_d[i][j] << "  " << invArelolo_d[i][j] << endl;
            cout << "invA_d[" << i << "][" << j << "]im : "
                 << invAimhihi_d[i][j] << "  " << invAimlohi_d[i][j] << endl
                 << "                 "
                 << invAimhilo_d[i][j] << "  " << invAimlolo_d[i][j] << endl;
         }
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on inverse : "
        << cmplx4_Matrix_Difference_Sum
              (dim,invArehihi_h,invArelohi_h,invArehilo_h,invArelolo_h,
                   invAimhihi_h,invAimlohi_h,invAimhilo_h,invAimlolo_h,
                   invArehihi_d,invArelohi_d,invArehilo_d,invArelolo_d,
                   invAimhihi_d,invAimlohi_d,invAimhilo_d,invAimlolo_d)
        << endl;

   double *xrehihi = new double[dim];
   double *xrelohi = new double[dim];
   double *xrehilo = new double[dim];
   double *xrelolo = new double[dim];
   double *ximhihi = new double[dim];
   double *ximlohi = new double[dim];
   double *ximhilo = new double[dim];
   double *ximlolo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      xrehihi[i] = 0.0; xrelohi[i] = 0.0;
      xrehilo[i] = 0.0; xrelolo[i] = 0.0;
      ximhihi[i] = 0.0; ximlohi[i] = 0.0;
      ximhilo[i] = 0.0; ximlolo[i] = 0.0;

      for(int j=0; j<dim; j++) // x[i] = x[i] + invA_h[i][j]*rhs[j];
      {
         qdf_mul(invArehihi_h[i][j],invArelohi_h[i][j],
                 invArehilo_h[i][j],invArelolo_h[i][j],
                 rhsrehihi[j],rhsrelohi[j],rhsrehilo[j],rhsrelolo[j],
                 &acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo);
         qdf_mul(invAimhihi_h[i][j],invAimlohi_h[i][j],
                 invAimhilo_h[i][j],invAimlolo_h[i][j],
                 rhsimhihi[j],rhsimlohi[j],rhsimhilo[j],rhsimlolo[j],
                 &acc2hihi,&acc2lohi,&acc2hilo,&acc2lolo);
         qdf_mul(invAimhihi_h[i][j],invAimlohi_h[i][j],
                 invAimhilo_h[i][j],invAimlolo_h[i][j],
                 rhsrehihi[j],rhsrelohi[j],rhsrehilo[j],rhsrelolo[j],
                 &acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo);
         qdf_mul(invArehihi_h[i][j],invArelohi_h[i][j],
                 invArehilo_h[i][j],invArelolo_h[i][j],
                 rhsimhihi[j],rhsimlohi[j],rhsimhilo[j],rhsimlolo[j],
                 &acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo);
         qdf_inc(&xrehihi[i],&xrelohi[i],&xrehilo[i],&xrelolo[i],
                 acc1hihi,acc1lohi,acc1hilo,acc1lolo);
         qdf_dec(&xrehihi[i],&xrelohi[i],&xrehilo[i],&xrelolo[i],
                 acc2hihi,acc2lohi,acc2hilo,acc2lolo);
         qdf_inc(&ximhihi[i],&ximlohi[i],&ximhilo[i],&ximlolo[i],
                 acc3hihi,acc3lohi,acc3hilo,acc3lolo);
         qdf_inc(&ximhihi[i],&ximlohi[i],&ximhilo[i],&ximlolo[i],
                 acc4hihi,acc4lohi,acc4hilo,acc4lolo);
      }
   }
   if(verbose > 0)
   {
      cout << "The solution computed with the CPU inverse :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehihi[i] << "  " << xrelohi[i] << endl
              << "         "
              << xrehilo[i] << "  " << xrelolo[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhihi[i] << "  " << ximlohi[i] << endl
              << "         "
              << ximhilo[i] << "  " << ximlolo[i] << endl;
      }
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on solution : "
        << cmplx4_Difference_Sum
              (dim,solrehihi,solrelohi,solrehilo,solrelolo,
                   solimhihi,solimlohi,solimhilo,solimlolo,
                     xrehihi,  xrelohi,  xrehilo,  xrelolo,
                     ximhihi,  ximlohi,  ximhilo,  ximlolo) << endl;
   cout << "Condition number : "
        << cmplx4_condition(dim, 
                 Arehihi,     Arelohi,     Arehilo,     Arelolo, 
                 Aimhihi,     Aimlohi,     Aimhilo,     Aimlolo,
              invArehihi_h,invArelohi_h,invArehilo_h,invArelolo_h,
              invAimhihi_h,invAimlohi_h,invAimhilo_h,invAimlolo_h)
        << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
   cout << "                     Time spent by the kernel : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
}

void test_real4_upper_tiling ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   const int dim = sizetile*numtiles;

   cout << "Generate a random matrix (1 = yes, 0 = read matrix) : ";
   int rndmat; cin >> rndmat;

   if(rndmat == 1)
      cout << "-> generating a random upper triangular matrix of dimension "
           << dim << " ..." << endl;
   else
      cout << "-> reading a random upper triangular matrix of dimension "
           << dim << " ..." << endl;
      
   double **Ahihi = new double*[dim];
   double **Alohi = new double*[dim];
   double **Ahilo = new double*[dim];
   double **Alolo = new double*[dim];

   for(int i=0; i<dim; i++) 
   {
      Ahihi[i] = new double[dim];
      Alohi[i] = new double[dim];
      Ahilo[i] = new double[dim];
      Alolo[i] = new double[dim];
      for(int j=0; j<dim; j++)
      {
         Alohi[i][j] = 0.0;
         Ahilo[i][j] = 0.0;
         Alolo[i][j] = 0.0;
      }
   }
   // random_dbl_upper_matrix(dim,dim,A);
   // dbl4_random_upper_factor(dim,Ahihi,Alohi,Ahilo,Alolo);
   if(rndmat == 1)
      dbl_random_upper_factor(dim,Ahihi);
   else
   {
      cout << "Give the name of a file : ";
      string filename; cin >> filename;
      cout << "-> reading " << dim*dim
           << " numbers from " << filename << " ..." << endl;
      dbl_read_matrix(filename,dim,Ahihi);
   }
   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihi[i][j] << "  " << Alohi[i][j] << endl
                 << "          "
                 << Ahilo[i][j] << "  " << Alolo[i][j] << endl;
   }
   double *solhihi = new double[dim];
   double *sollohi = new double[dim];
   double *solhilo = new double[dim];
   double *sollolo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solhihi[i] = 1.0; sollohi[i] = 0.0;
      solhilo[i] = 0.0; sollolo[i] = 0.0;
   }
   double *xhihi = new double[dim];
   double *xlohi = new double[dim];
   double *xhilo = new double[dim];
   double *xlolo = new double[dim];
   double *rhshihi = new double[dim];
   double *rhslohi = new double[dim];
   double *rhshilo = new double[dim];
   double *rhslolo = new double[dim];
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<dim; i++)
   {
      rhshihi[i] = 0.0; rhslohi[i] = 0.0;
      rhshilo[i] = 0.0; rhslolo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         qdf_mul(Ahihi[i][j],Alohi[i][j],Ahilo[i][j],Alolo[i][j],
               solhihi[j], sollohi[j], solhilo[j], sollolo[j],
              &acchihi,   &acclohi,   &acchilo,  &acclolo);
         qdf_inc(&rhshihi[i],&rhslohi[i],&rhshilo[i],&rhslolo[i],
                  acchihi,    acclohi,    acchilo,   acclolo);
      }
   }
   double *xhihi_d = new double[dim];
   double *xlohi_d = new double[dim];
   double *xhilo_d = new double[dim];
   double *xlolo_d = new double[dim];
   double *rhshihi_d = new double[dim];
   double *rhslohi_d = new double[dim];
   double *rhshilo_d = new double[dim];
   double *rhslolo_d = new double[dim];
   double **Ahihi_d = new double*[dim];
   double **Alohi_d = new double*[dim];
   double **Ahilo_d = new double*[dim];
   double **Alolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      rhshihi_d[i] = rhshihi[i];
      rhslohi_d[i] = rhslohi[i];
      rhshilo_d[i] = rhshilo[i];
      rhslolo_d[i] = rhslolo[i];

      Ahihi_d[i] = new double[dim];
      Alohi_d[i] = new double[dim];
      Ahilo_d[i] = new double[dim];
      Alolo_d[i] = new double[dim];

      for(int j=0; j<dim; j++)
      {
         Ahihi_d[i][j] = Ahihi[i][j];
         Alohi_d[i][j] = Alohi[i][j];
         Ahilo_d[i][j] = Ahilo[i][j];
         Alolo_d[i][j] = Alolo[i][j];
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : "
              << rhshihi[i] << "  " << rhslohi[i] << endl
              << "       "
              << rhshilo[i] << "  " << rhslolo[i] << endl;
   }
   double timelapsed_h,timelapsed_d,elapsedms;
   double invlapsed,mullapsed,sublapsed;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;
   double addover = 0.0;
   double mulover = 0.0;
   double divover = 0.0;

   cout << "-> CPU solves an upper triangular system ..." << endl;

   CPU_dbl4_upper_tiled_solver
      (dim,sizetile,numtiles,Ahihi,Alohi,Ahilo,Alolo,
       rhshihi,rhslohi,rhshilo,rhslolo,xhihi,xlohi,xhilo,xlolo,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihi[i][j] << "  " << Alohi[i][j] << endl
                 << "          "
                 << Ahilo[i][j] << "  " << Alolo[i][j] << endl;
   }

   cout << "-> GPU solves an upper triangular system ..." << endl;

   GPU_dbl4_upper_tiled_solver
      (dim,sizetile,numtiles,
         Ahihi_d,  Alohi_d,  Ahilo_d,  Alolo_d,
       rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
         xhihi_d,  xlohi_d,  xhilo_d,  xlolo_d,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&timelapsed_d,
       &addcnt,&addover,&mulcnt,&mulover,&divcnt,&divover);

   if(verbose > 0)
   {
      cout << "The matrix returned by the device :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihi_d[i][j] << "  " << Alohi_d[i][j] << endl
                 << "          "
                 << Ahilo_d[i][j] << "  " << Alolo_d[i][j] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on diagonal tiles : "
        << dbl4_Diagonal_Difference_Sum
              (numtiles,sizetile,
               Ahihi,Alohi,Ahilo,Alolo,Ahihi_d,Alohi_d,Ahilo_d,Alolo_d)
        << endl;

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "CPU solution computed with tiling :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhihi[i] << "  " << xlohi[i] << endl
              << "       "
              << xhilo[i] << "  " << xlolo[i] << endl;

      cout << "GPU solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhihi_d[i] << "  " << xlohi_d[i] << endl
              << "       "
              << xhilo_d[i] << "  " << xlolo_d[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of CPU errors on solution : "
        << dbl4_Difference_Sum(dim,solhihi,sollohi,solhilo,sollolo,
                                     xhihi,  xlohi,  xhilo,  xlolo)
        << endl;
   cout << "   Sum of GPU errors on solution : "
        << dbl4_Difference_Sum(dim,solhihi,sollohi,solhilo,sollolo,
                                     xhihi_d,xlohi_d,xhilo_d,xlolo_d)
        << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;

   write_dbl4_bstimeflops
     (sizetile,numtiles,0,invlapsed,mullapsed,sublapsed,elapsedms,
      timelapsed_d,addcnt,addover,mulcnt,mulover,divcnt,divover);

   for(int i=0; i<dim; i++)
   {
      free(Ahihi[i]); free(Alohi[i]); free(Ahilo[i]); free(Alolo[i]);
      free(Ahihi_d[i]); free(Alohi_d[i]); free(Ahilo_d[i]); free(Alolo_d[i]);
   }
   free(Ahihi); free(Alohi); free(Ahilo); free(Alolo);
   free(Ahihi_d); free(Alohi_d); free(Ahilo_d); free(Alolo_d);
   free(solhihi); free(sollohi); free(solhilo); free(sollolo);
   free(rhshihi); free(rhslohi); free(rhshilo); free(rhslolo);
   free(xhihi); free(xlohi); free(xhilo); free(xlolo); 
   free(rhshihi_d); free(rhslohi_d); free(rhshilo_d); free(rhslolo_d);
   free(xhihi_d); free(xlohi_d); free(xhilo_d); free(xlolo_d);
}

void test_cmplx4_upper_tiling ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   const int dim = sizetile*numtiles;

   cout << "Generate a random matrix (1 = yes, 0 = read matrix) : ";
   int rndmat; cin >> rndmat;

   if(rndmat == 1)
      cout << "-> generating a random upper triangular matrix of dimension "
           << dim << " ..." << endl;
   else
      cout << "-> reading a random upper triangular matrix of dimension "
           << dim << " ..." << endl;

   double **Arehihi = new double*[dim];
   double **Arelohi = new double*[dim];
   double **Arehilo = new double*[dim];
   double **Arelolo = new double*[dim];
   double **Aimhihi = new double*[dim];
   double **Aimlohi = new double*[dim];
   double **Aimhilo = new double*[dim];
   double **Aimlolo = new double*[dim];

   for(int i=0; i<dim; i++) 
   {
      Arehihi[i] = new double[dim];
      Arelohi[i] = new double[dim];
      Arehilo[i] = new double[dim];
      Arelolo[i] = new double[dim];
      Aimhihi[i] = new double[dim];
      Aimlohi[i] = new double[dim];
      Aimhilo[i] = new double[dim];
      Aimlolo[i] = new double[dim];
      for(int j=0; j<dim; j++)
      {
         Arelohi[i][j] = 0.0;
         Arehilo[i][j] = 0.0;
         Arelolo[i][j] = 0.0;
         Aimlohi[i][j] = 0.0;
         Aimhilo[i][j] = 0.0;
         Aimlolo[i][j] = 0.0;
      }
   }
   // random_dbl_upper_matrix(dim,dim,Arehi,Arelo,Aimhi,Aimlo);
   // cmplx4_random_upper_factor
   //    (dim,Arehihi,Arelohi,Arehilo,Arelolo,Aimhihi,Aimlohi,Aimhilo,Aimlolo);
   // cmplx_random_upper_factor(dim,Arehihi,Aimhihi);

   if(rndmat == 1)
      cmplx_random_upper_factor(dim,Arehihi,Aimhihi);
   else
   {
      cout << "Give the name of a file : ";
      string filename; cin >> filename;
      cout << "-> reading " << dim*dim
           << " numbers from " << filename << " ..." << endl;
      cmplx_read_matrix(filename,dim,Arehihi,Aimhihi);
   }
   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihi[i][j] << "  " << Arelohi[i][j] << endl
                 << "            "
                 << Arehilo[i][j] << "  " << Arelolo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihi[i][j] << "  " << Aimlohi[i][j] << endl
                 << "            "
                 << Aimhilo[i][j] << "  " << Aimlolo[i][j] << endl;
         }
   }
   double *solrehihi = new double[dim];
   double *solrelohi = new double[dim];
   double *solrehilo = new double[dim];
   double *solrelolo = new double[dim];
   double *solimhihi = new double[dim];
   double *solimlohi = new double[dim];
   double *solimhilo = new double[dim];
   double *solimlolo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solrehihi[i] = 1.0; solrelohi[i] = 0.0;
      solrehilo[i] = 0.0; solrelolo[i] = 0.0;
      solimhihi[i] = 0.0; solimlohi[i] = 0.0;
      solimhilo[i] = 0.0; solimlolo[i] = 0.0;
   }
   double *xrehihi = new double[dim];
   double *xrelohi = new double[dim];
   double *xrehilo = new double[dim];
   double *xrelolo = new double[dim];
   double *ximhihi = new double[dim];
   double *ximlohi = new double[dim];
   double *ximhilo = new double[dim];
   double *ximlolo = new double[dim];
   double *rhsrehihi = new double[dim];
   double *rhsrelohi = new double[dim];
   double *rhsrehilo = new double[dim];
   double *rhsrelolo = new double[dim];
   double *rhsimhihi = new double[dim];
   double *rhsimlohi = new double[dim];
   double *rhsimhilo = new double[dim];
   double *rhsimlolo = new double[dim];
   double acc1hihi,acc1lohi,acc1hilo,acc1lolo;
   double acc2hihi,acc2lohi,acc2hilo,acc2lolo;
   double acc3hihi,acc3lohi,acc3hilo,acc3lolo;
   double acc4hihi,acc4lohi,acc4hilo,acc4lolo;

   for(int i=0; i<dim; i++)
   {
      rhsrehihi[i] = 0.0; rhsrelohi[i] = 0.0;
      rhsrehilo[i] = 0.0; rhsrelolo[i] = 0.0;
      rhsimhihi[i] = 0.0; rhsimlohi[i] = 0.0;
      rhsimhilo[i] = 0.0; rhsimlolo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
               solrehihi[j] ,solrelohi[j], solrehilo[j], solrelolo[j],
               &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
         qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
               solimhihi[j], solimlohi[j], solimhilo[j], solimlolo[j],
               &acc2hihi,    &acc2lohi,    &acc2hilo,   &acc2lolo);
         qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
               solrehihi[j], solrelohi[j], solrehilo[j], solrelolo[j],
               &acc3hihi,    &acc3lohi,    &acc3hilo,    &acc3lolo);
         qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
               solimhihi[j], solimlohi[j], solimhilo[j], solimlolo[j],
               &acc4hihi,    &acc4lohi,    &acc4hilo,    &acc4lolo);
         qdf_dec(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
                  acc2hihi, acc2lohi, acc2hilo, acc2lolo);
         qdf_inc(&rhsrehihi[i],&rhsrelohi[i],&rhsrehilo[i],&rhsrelolo[i],
                   acc1hihi,     acc1lohi,     acc1hilo,     acc1lolo);
         qdf_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
                  acc4hihi, acc4lohi, acc4hilo, acc4lolo);
         qdf_inc(&rhsimhihi[i],&rhsimlohi[i],&rhsimhilo[i],&rhsimlolo[i],
                   acc3hihi,     acc3lohi,     acc3hilo,     acc3lolo);
      }
   }
   double *xrehihi_d = new double[dim];
   double *xrelohi_d = new double[dim];
   double *xrehilo_d = new double[dim];
   double *xrelolo_d = new double[dim];
   double *ximhihi_d = new double[dim];
   double *ximlohi_d = new double[dim];
   double *ximhilo_d = new double[dim];
   double *ximlolo_d = new double[dim];
   double *rhsrehihi_d = new double[dim];
   double *rhsrelohi_d = new double[dim];
   double *rhsrehilo_d = new double[dim];
   double *rhsrelolo_d = new double[dim];
   double *rhsimhihi_d = new double[dim];
   double *rhsimlohi_d = new double[dim];
   double *rhsimhilo_d = new double[dim];
   double *rhsimlolo_d = new double[dim];
   double **Arehihi_d = new double*[dim];
   double **Arelohi_d = new double*[dim];
   double **Arehilo_d = new double*[dim];
   double **Arelolo_d = new double*[dim];
   double **Aimhihi_d = new double*[dim];
   double **Aimlohi_d = new double*[dim];
   double **Aimhilo_d = new double*[dim];
   double **Aimlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      rhsrehihi_d[i] = rhsrehihi[i];
      rhsrelohi_d[i] = rhsrelohi[i];
      rhsrehilo_d[i] = rhsrehilo[i];
      rhsrelolo_d[i] = rhsrelolo[i];
      rhsimhihi_d[i] = rhsimhihi[i];
      rhsimlohi_d[i] = rhsimlohi[i];
      rhsimhilo_d[i] = rhsimhilo[i];
      rhsimlolo_d[i] = rhsimlolo[i];
      Arehihi_d[i] = new double[dim];
      Arelohi_d[i] = new double[dim];
      Arehilo_d[i] = new double[dim];
      Arelolo_d[i] = new double[dim];
      Aimhihi_d[i] = new double[dim];
      Aimlohi_d[i] = new double[dim];
      Aimhilo_d[i] = new double[dim];
      Aimlolo_d[i] = new double[dim];

      for(int j=0; j<dim; j++)
      {
         Arehihi_d[i][j] = Arehihi[i][j];
         Arelohi_d[i][j] = Arelohi[i][j];
         Arehilo_d[i][j] = Arehilo[i][j];
         Arelolo_d[i][j] = Arelolo[i][j];
         Aimhihi_d[i][j] = Aimhihi[i][j];
         Aimlohi_d[i][j] = Aimlohi[i][j];
         Aimhilo_d[i][j] = Aimhilo[i][j];
         Aimlolo_d[i][j] = Aimlolo[i][j];
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehihi[i] << "  " << rhsrelohi[i] << endl
              << "         "
              << rhsrehilo[i] << "  " << rhsrelolo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhihi[i] << "  " << rhsimlohi[i] << endl
              << "         "
              << rhsimhilo[i] << "  " << rhsimlolo[i] << endl;
      }
   }
   double timelapsed_h,timelapsed_d,elapsedms;
   double invlapsed,mullapsed,sublapsed;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;
   double addover = 0.0;
   double mulover = 0.0;
   double divover = 0.0;

   cout << "-> CPU solves an upper triangular system ..." << endl;

   CPU_cmplx4_upper_tiled_solver
      (dim,sizetile,numtiles,
       Arehihi,Arelohi,Arehilo,Arelolo,Aimhihi,Aimlohi,Aimhilo,Aimlolo,
       rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
       rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
       xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
       &timelapsed_h);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihi[i][j] << "  " << Arelohi[i][j] << endl
                 << "            "
                 << Arehilo[i][j] << "  " << Arelolo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihi[i][j] << "  " << Aimlohi[i][j] << endl
                 << "            "
                 << Aimhilo[i][j] << "  " << Aimlolo[i][j] << endl;
         }
   }
   cout << "-> GPU solves an upper triangular system ..." << endl;

   GPU_cmplx4_upper_tiled_solver
      (dim,sizetile,numtiles,
         Arehihi_d,  Arelohi_d,  Arehilo_d,  Arelolo_d,
         Aimhihi_d,  Aimlohi_d,  Aimhilo_d,  Aimlolo_d,
       rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
       rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,
         xrehihi_d,  xrelohi_d,  xrehilo_d,  xrelolo_d,
         ximhihi_d,  ximlohi_d,  ximhilo_d,  ximlolo_d,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&timelapsed_d,
       &addcnt,&addover,&mulcnt,&mulover,&divcnt,&divover);

   if(verbose > 0)
   {
      cout << "The matrix returned by the device :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihi_d[i][j] << "  " << Arelohi_d[i][j] << endl
                 << "            "
                 << Arehilo_d[i][j] << "  " << Arelolo_d[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihi_d[i][j] << "  " << Aimlohi_d[i][j] << endl
                 << "            "
                 << Aimhilo_d[i][j] << "  " << Aimlolo_d[i][j] << endl;
         }
   }

   cout << scientific << setprecision(2);
   cout << "   Sum of errors on diagonal tiles : "
        << cmplx4_Diagonal_Difference_Sum
             (numtiles,sizetile,
              Arehihi,  Arelohi,  Arehilo,  Arelolo,
              Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,
              Arehihi_d,Arelohi_d,Arehilo_d,Arelolo_d,
              Aimhihi_d,Aimlohi_d,Aimhilo_d,Aimlolo_d)
        << endl;

   if(verbose > 0)
   {
      cout << "CPU solution computed with tiling :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehihi[i] << "  " << xrelohi[i] << endl
              << "         "
              << xrehilo[i] << "  " << xrelolo[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhihi[i] << "  " << ximlohi[i] << endl
              << "         "
              << ximhilo[i] << "  " << ximlolo[i] << endl;
      }
      cout << "GPU solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehihi_d[i] << "  " << xrelohi_d[i] << endl
              << "         "
              << xrehilo_d[i] << "  " << xrelolo_d[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhihi_d[i] << "  " << ximlohi_d[i] << endl
              << "         "
              << ximhilo_d[i] << "  " << ximlolo_d[i] << endl;
      }
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of CPU errors on solution : "
        << cmplx4_Difference_Sum(dim,
              solrehihi,solrelohi,solrehilo,solrelolo,
              solimhihi,solimlohi,solimhilo,solimlolo,
                xrehihi,  xrelohi,  xrehilo,  xrelolo, 
                ximhihi,  ximlohi,  ximhilo,  ximlolo)
        << endl;
   cout << "   Sum of GPU errors on solution : "
        << cmplx4_Difference_Sum(dim,
               solrehihi,solrelohi,solrehilo,solrelolo,
               solimhihi,solimlohi,solimhilo,solimlolo,
                 xrehihi_d,xrelohi_d,xrehilo_d,xrelolo_d,
                 ximhihi_d,ximlohi_d,ximhilo_d,ximlolo_d)
        << endl;
   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;

   write_dbl4_bstimeflops
     (sizetile,numtiles,1,invlapsed,mullapsed,sublapsed,elapsedms,
      timelapsed_d,addcnt,addover,mulcnt,mulover,divcnt,divover);

   for(int i=0; i<dim; i++)
   {
      free(Arehihi[i]); free(Arelohi[i]);
      free(Arehilo[i]); free(Arelolo[i]);
      free(Aimhihi[i]); free(Aimlohi[i]);
      free(Aimhilo[i]); free(Aimlolo[i]);
      free(Arehihi_d[i]); free(Arelohi_d[i]);
      free(Arehilo_d[i]); free(Arelolo_d[i]);
      free(Aimhihi_d[i]); free(Aimlohi_d[i]);
      free(Aimhilo_d[i]); free(Aimlolo_d[i]);
   }
   free(Arehihi); free(Arelohi); free(Arehilo); free(Arelolo);
   free(Aimhihi); free(Aimlohi); free(Aimhilo); free(Aimlolo);
   free(Arehihi_d); free(Arelohi_d); free(Arehilo_d); free(Arelolo_d);
   free(Aimhihi_d); free(Aimlohi_d); free(Aimhilo_d); free(Aimlolo_d);
   free(solrehihi); free(solrelohi); free(solrehilo); free(solrelolo);
   free(solimhihi); free(solimlohi); free(solimhilo); free(solimlolo);
   free(rhsrehihi); free(rhsrelohi); free(rhsrehilo); free(rhsrelolo);
   free(rhsimhihi); free(rhsimlohi); free(rhsimhilo); free(rhsimlolo);
   free(rhsrehihi_d); free(rhsrelohi_d); free(rhsrehilo_d); free(rhsrelolo_d);
   free(rhsimhihi_d); free(rhsimlohi_d); free(rhsimhilo_d); free(rhsimlolo_d);
   free(xrehihi); free(xrelohi); free(xrehilo); free(xrelolo);
   free(ximhihi); free(ximlohi); free(ximhilo); free(ximlolo);
   free(xrehihi_d); free(xrelohi_d); free(xrehilo_d); free(xrelolo_d);
   free(ximhihi_d); free(ximlohi_d); free(ximhilo_d); free(ximlolo_d);
}
