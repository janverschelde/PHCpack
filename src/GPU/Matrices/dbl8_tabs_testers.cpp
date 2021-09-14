// The file dbl8_tabs_testers.cpp defines the functions specified in
// the file dbl8_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
// #include <vector_types.h>
#include "octo_double_functions.h"
#include "random8_matrices.h"
#include "dbl8_factorizations.h"
#include "dbl8_tabs_host.h"
// #include "dbl8_tabs_kernels.h"
#include "dbl8_test_utilities.h"

using namespace std;

void test_real8_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Ahihihi = new double*[dim];
   double **Alohihi = new double*[dim];
   double **Ahilohi = new double*[dim];
   double **Alolohi = new double*[dim];
   double **Ahihilo = new double*[dim];
   double **Alohilo = new double*[dim];
   double **Ahilolo = new double*[dim];
   double **Alololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Ahihihi[i] = new double[dim];
      Alohihi[i] = new double[dim];
      Ahilohi[i] = new double[dim];
      Alolohi[i] = new double[dim];
      Ahihilo[i] = new double[dim];
      Alohilo[i] = new double[dim];
      Ahilolo[i] = new double[dim];
      Alololo[i] = new double[dim];
   }
 /*
   random_dbl8_upper_matrix
      (dim,dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
               Ahihilo,Alohilo,Ahilolo,Alololo);
  */
   dbl8_random_upper_factor
      (dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
           Ahihilo,Alohilo,Ahilolo,Alololo);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihihi[i][j] << "  " << Alohihi[i][j] << endl
                 << "          "
                 << Ahilohi[i][j] << "  " << Alolohi[i][j] << endl
                 << "          "
                 << Ahihilo[i][j] << "  " << Alohilo[i][j] << endl
                 << "          "
                 << Ahilolo[i][j] << "  " << Alololo[i][j] << endl;
   }
   double *solhihihi = new double[dim];
   double *sollohihi = new double[dim];
   double *solhilohi = new double[dim];
   double *sollolohi = new double[dim];
   double *solhihilo = new double[dim];
   double *sollohilo = new double[dim];
   double *solhilolo = new double[dim];
   double *sollololo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solhihihi[i] = 1.0; sollohihi[i] = 0.0;
      solhilohi[i] = 0.0; sollolohi[i] = 0.0;
      solhihilo[i] = 0.0; sollohilo[i] = 0.0;
      solhilolo[i] = 0.0; sollololo[i] = 0.0;
   }
   double *rhshihihi = new double[dim];
   double *rhslohihi = new double[dim];
   double *rhshilohi = new double[dim];
   double *rhslolohi = new double[dim];
   double *rhshihilo = new double[dim];
   double *rhslohilo = new double[dim];
   double *rhshilolo = new double[dim];
   double *rhslololo = new double[dim];
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<dim; i++)
   {
      rhshihihi[i] = 0.0; rhslohihi[i] = 0.0;
      rhshilohi[i] = 0.0; rhslolohi[i] = 0.0;
      rhshihilo[i] = 0.0; rhslohilo[i] = 0.0;
      rhshilolo[i] = 0.0; rhslololo[i] = 0.0;

      for(int j=0; j<dim; j++)  // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         odf_mul(Ahihihi[i][j],Alohihi[i][j],Ahilohi[i][j],Alolohi[i][j],
                 Ahihilo[i][j],Alohilo[i][j],Ahilolo[i][j],Alololo[i][j],
               solhihihi[j], sollohihi[j], solhilohi[j], sollolohi[j],
               solhihilo[j], sollohilo[j], solhilolo[j], sollololo[j],
              &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
              &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
         odf_inc(&rhshihihi[i],&rhslohihi[i],&rhshilohi[i],&rhslolohi[i],
                 &rhshihilo[i],&rhslohilo[i],&rhshilolo[i],&rhslololo[i],
                  acchihihi,    acclohihi,    acchilohi,    acclolohi,
                  acchihilo,    acclohilo,    acchilolo,    acclololo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : "
              << rhshihihi[i] << "  " << rhslohihi[i] << endl
              << "       "
              << rhshilohi[i] << "  " << rhslolohi[i] << endl
              << "       "
              << rhshihilo[i] << "  " << rhslohilo[i] << endl
              << "       "
              << rhshilolo[i] << "  " << rhslololo[i] << endl;
   }
   double **invAhihihi_h = new double*[dim];
   double **invAlohihi_h = new double*[dim];
   double **invAhilohi_h = new double*[dim];
   double **invAlolohi_h = new double*[dim];
   double **invAhihilo_h = new double*[dim];
   double **invAlohilo_h = new double*[dim];
   double **invAhilolo_h = new double*[dim];
   double **invAlololo_h = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      invAhihihi_h[i] = new double[dim];
      invAlohihi_h[i] = new double[dim];
      invAhilohi_h[i] = new double[dim];
      invAlolohi_h[i] = new double[dim];
      invAhihilo_h[i] = new double[dim];
      invAlohilo_h[i] = new double[dim];
      invAhilolo_h[i] = new double[dim];
      invAlololo_h[i] = new double[dim];
   }
   double timelapsed_h,timelapsed_d,elapsedms;

   cout << "-> CPU computes the inverse ..." << endl;

   CPU_dbl8_upper_inverse
      (dim,Ahihihi,     Alohihi,     Ahilohi,     Alolohi,
           Ahihilo,     Alohilo,     Ahilolo,     Alololo,
        invAhihihi_h,invAlohihi_h,invAhilohi_h,invAlolohi_h,
        invAhihilo_h,invAlohilo_h,invAhilolo_h,invAlololo_h,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The CPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_h[" << i << "][" << j << "] : "
                 << invAhihihi_h[i][j] << "  " << invAlohihi_h[i][j] << endl
                 << "               "
                 << invAhilohi_h[i][j] << "  " << invAlolohi_h[i][j] << endl
                 << "               "
                 << invAhihilo_h[i][j] << "  " << invAlohilo_h[i][j] << endl
                 << "               "
                 << invAhilolo_h[i][j] << "  " << invAlololo_h[i][j] << endl;
   }
/*
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
 */
   double *xhihihi = new double[dim];
   double *xlohihi = new double[dim];
   double *xhilohi = new double[dim];
   double *xlolohi = new double[dim];
   double *xhihilo = new double[dim];
   double *xlohilo = new double[dim];
   double *xhilolo = new double[dim];
   double *xlololo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      xhihihi[i] = 0.0; xlohihi[i] = 0.0;
      xhilohi[i] = 0.0; xlolohi[i] = 0.0;
      xhihilo[i] = 0.0; xlohilo[i] = 0.0;
      xhilolo[i] = 0.0; xlololo[i] = 0.0;

      for(int j=0; j<dim; j++)   // x[i] = x[i] + invA_h[i][j]*rhs[j];
      {
         odf_mul(invAhihihi_h[i][j],invAlohihi_h[i][j],
                 invAhilohi_h[i][j],invAlolohi_h[i][j],
                 invAhihilo_h[i][j],invAlohilo_h[i][j],
                 invAhilolo_h[i][j],invAlololo_h[i][j],
                 rhshihihi[j],rhslohihi[j],rhshilohi[j],rhslolohi[j],
                 rhshihilo[j],rhslohilo[j],rhshilolo[j],rhslololo[j],
                &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odf_inc(&xhihihi[i],&xlohihi[i],&xhilohi[i],&xlolohi[i],
                 &xhihilo[i],&xlohilo[i],&xhilolo[i],&xlololo[i],
                acchihihi,  acclohihi,  acchilohi,  acclolohi,
                acchihilo,  acclohilo,  acchilolo,  acclololo);
      }
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The solution computed with the CPU inverse :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhihihi[i] << "  " << xlohihi[i] << endl
              << "       "
              << xhilohi[i] << "  " << xlolohi[i] << endl
              << "       "
              << xhihilo[i] << "  " << xlohilo[i] << endl
              << "       "
              << xhilolo[i] << "  " << xlololo[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on solution : "
        << dbl8_Difference_Sum
              (dim, solhihihi,sollohihi,solhilohi,sollolohi,
                    solhihilo,sollohilo,solhilolo,sollololo,
                      xhihihi,  xlohihi,  xhilohi,  xlolohi,
                      xhihilo,  xlohilo,  xhilolo,  xlololo) << endl;
   cout << "Condition number : "
        << dbl8_condition
              (dim,Ahihihi,     Alohihi,      Ahilohi,     Alolohi,
                   Ahihilo,     Alohilo,      Ahilolo,     Alololo,
                invAhihihi_h,invAlohihi_h, invAhilohi_h,invAlolohi_h,
                invAhihilo_h,invAlohilo_h, invAhilolo_h,invAlololo_h)
        << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
/*
   cout << "                     Time spent by the kernel : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
 */
}

void test_cmplx8_upper_inverse ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;


   double **Arehihihi = new double*[dim];
   double **Arelohihi = new double*[dim];
   double **Arehilohi = new double*[dim];
   double **Arelolohi = new double*[dim];
   double **Arehihilo = new double*[dim];
   double **Arelohilo = new double*[dim];
   double **Arehilolo = new double*[dim];
   double **Arelololo = new double*[dim];
   double **Aimhihihi = new double*[dim];
   double **Aimlohihi = new double*[dim];
   double **Aimhilohi = new double*[dim];
   double **Aimlolohi = new double*[dim];
   double **Aimhihilo = new double*[dim];
   double **Aimlohilo = new double*[dim];
   double **Aimhilolo = new double*[dim];
   double **Aimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Arehihihi[i] = new double[dim];
      Arelohihi[i] = new double[dim];
      Arehilohi[i] = new double[dim];
      Arelolohi[i] = new double[dim];
      Arehihilo[i] = new double[dim];
      Arelohilo[i] = new double[dim];
      Arehilolo[i] = new double[dim];
      Arelololo[i] = new double[dim];
      Aimhihihi[i] = new double[dim];
      Aimlohihi[i] = new double[dim];
      Aimhilohi[i] = new double[dim];
      Aimlolohi[i] = new double[dim];
      Aimhihilo[i] = new double[dim];
      Aimlohilo[i] = new double[dim];
      Aimhilolo[i] = new double[dim];
      Aimlololo[i] = new double[dim];
   }
  /*
   random_cmplx8_upper_matrix
      (dim,dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
               Arehihilo,Arelohilo,Arehilolo,Arelololo,
               Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
               Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);
   */
   cmplx8_random_upper_factor
      (dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
           Arehihilo,Arelohilo,Arehilolo,Arelololo,
           Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
           Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihihi[i][j] << "  " << Arelohihi[i][j] << endl
                 << "            "
                 << Arehilohi[i][j] << "  " << Arelolohi[i][j] << endl
                 << "            "
                 << Arehihilo[i][j] << "  " << Arelohilo[i][j] << endl
                 << "            "
                 << Arehilolo[i][j] << "  " << Arelololo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihihi[i][j] << "  " << Aimlohihi[i][j] << endl
                 << "            "
                 << Aimhilohi[i][j] << "  " << Aimlolohi[i][j] << endl
                 << "            "
                 << Aimhihilo[i][j] << "  " << Aimlohilo[i][j] << endl
                 << "            "
                 << Aimhilolo[i][j] << "  " << Aimlololo[i][j] << endl;
         }
   }
   double *solrehihihi = new double[dim];
   double *solrelohihi = new double[dim];
   double *solrehilohi = new double[dim];
   double *solrelolohi = new double[dim];
   double *solrehihilo = new double[dim];
   double *solrelohilo = new double[dim];
   double *solrehilolo = new double[dim];
   double *solrelololo = new double[dim];
   double *solimhihihi = new double[dim];
   double *solimlohihi = new double[dim];
   double *solimhilohi = new double[dim];
   double *solimlolohi = new double[dim];
   double *solimhihilo = new double[dim];
   double *solimlohilo = new double[dim];
   double *solimhilolo = new double[dim];
   double *solimlololo = new double[dim];

   for(int i=0; i<dim; i++) 
   {
      solrehihihi[i] = 1.0; solrelohihi[i] = 0.0;
      solrehilohi[i] = 0.0; solrelolohi[i] = 0.0;
      solrehihilo[i] = 0.0; solrelohilo[i] = 0.0;
      solrehilolo[i] = 0.0; solrelololo[i] = 0.0;
      solimhihihi[i] = 0.0; solimlohihi[i] = 0.0;
      solimhilohi[i] = 0.0; solimlolohi[i] = 0.0;
      solimhihilo[i] = 0.0; solimlohilo[i] = 0.0;
      solimhilolo[i] = 0.0; solimlololo[i] = 0.0;
   }
   double *rhsrehihihi = new double[dim];
   double *rhsrelohihi = new double[dim];
   double *rhsrehilohi = new double[dim];
   double *rhsrelolohi = new double[dim];
   double *rhsrehihilo = new double[dim];
   double *rhsrelohilo = new double[dim];
   double *rhsrehilolo = new double[dim];
   double *rhsrelololo = new double[dim];
   double *rhsimhihihi = new double[dim];
   double *rhsimlohihi = new double[dim];
   double *rhsimhilohi = new double[dim];
   double *rhsimlolohi = new double[dim];
   double *rhsimhihilo = new double[dim];
   double *rhsimlohilo = new double[dim];
   double *rhsimhilolo = new double[dim];
   double *rhsimlololo = new double[dim];
   double acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi;
   double acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo;
   double acc2hihihi,acc2lohihi,acc2hilohi,acc2lolohi;
   double acc2hihilo,acc2lohilo,acc2hilolo,acc2lololo;
   double acc3hihihi,acc3lohihi,acc3hilohi,acc3lolohi;
   double acc3hihilo,acc3lohilo,acc3hilolo,acc3lololo;
   double acc4hihihi,acc4lohihi,acc4hilohi,acc4lolohi;
   double acc4hihilo,acc4lohilo,acc4hilolo,acc4lololo;

   for(int i=0; i<dim; i++)
   {
      rhsrehihihi[i] = 0.0; rhsrelohihi[i] = 0.0;
      rhsrehilohi[i] = 0.0; rhsrelolohi[i] = 0.0;
      rhsrehihilo[i] = 0.0; rhsrelohilo[i] = 0.0;
      rhsrehilolo[i] = 0.0; rhsrelololo[i] = 0.0;
      rhsimhihihi[i] = 0.0; rhsimlohihi[i] = 0.0;
      rhsimhilohi[i] = 0.0; rhsimlolohi[i] = 0.0;
      rhsimhihilo[i] = 0.0; rhsimlohilo[i] = 0.0;
      rhsimhilolo[i] = 0.0; rhsimlololo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         odf_mul(Arehihihi[i][j],Arelohihi[i][j],
                 Arehilohi[i][j],Arelolohi[i][j],
                 Arehihilo[i][j],Arelohilo[i][j],
                 Arehilolo[i][j],Arelololo[i][j],
               solrehihihi[j], solrelohihi[j], solrehilohi[j], solrelolohi[j],
               solrehihilo[j], solrelohilo[j], solrehilolo[j], solrelololo[j],
               &acc1hihihi,    &acc1lohihi,    &acc1hilohi,    &acc1lolohi,
               &acc1hihilo,    &acc1lohilo,    &acc1hilolo,    &acc1lololo);
         odf_mul(Aimhihihi[i][j],Aimlohihi[i][j],
                 Aimhilohi[i][j],Aimlolohi[i][j],
                 Aimhihilo[i][j],Aimlohilo[i][j],
                 Aimhilolo[i][j],Aimlololo[i][j],
               solimhihihi[j], solimlohihi[j], solimhilohi[j], solimlolohi[j],
               solimhihilo[j], solimlohilo[j], solimhilolo[j], solimlololo[j],
               &acc2hihihi,    &acc2lohihi,    &acc2hilohi,   &acc2lolohi,
               &acc2hihilo,    &acc2lohilo,    &acc2hilolo,   &acc2lololo);
         odf_mul(Aimhihihi[i][j],Aimlohihi[i][j],
                 Aimhilohi[i][j],Aimlolohi[i][j],
                 Aimhihilo[i][j],Aimlohilo[i][j],
                 Aimhilolo[i][j],Aimlololo[i][j],
               solrehihihi[j], solrelohihi[j], solrehilohi[j], solrelolohi[j],
               solrehihilo[j], solrelohilo[j], solrehilolo[j], solrelololo[j],
               &acc3hihihi,    &acc3lohihi,    &acc3hilohi,    &acc3lolohi,
               &acc3hihilo,    &acc3lohilo,    &acc3hilolo,    &acc3lololo);
         odf_mul(Arehihihi[i][j],Arelohihi[i][j],
                 Arehilohi[i][j],Arelolohi[i][j],
                 Arehihilo[i][j],Arelohilo[i][j],
                 Arehilolo[i][j],Arelololo[i][j],
               solimhihihi[j], solimlohihi[j], solimhilohi[j], solimlolohi[j],
               solimhihilo[j], solimlohilo[j], solimhilolo[j], solimlololo[j],
               &acc4hihihi,    &acc4lohihi,    &acc4hilohi,    &acc4lolohi,
               &acc4hihilo,    &acc4lohilo,    &acc4hilolo,    &acc4lololo);
         odf_inc(&rhsrehihihi[i],&rhsrelohihi[i],
                 &rhsrehilohi[i],&rhsrelolohi[i],
                 &rhsrehihilo[i],&rhsrelohilo[i],
                 &rhsrehilolo[i],&rhsrelololo[i],
                   acc1hihihi,     acc1lohihi,     acc1hilohi,     acc1lolohi,
                   acc1hihilo,     acc1lohilo,     acc1hilolo,     acc1lololo);
         odf_dec(&rhsrehihihi[i],&rhsrelohihi[i],
                 &rhsrehilohi[i],&rhsrelolohi[i],
                 &rhsrehihilo[i],&rhsrelohilo[i],
                 &rhsrehilolo[i],&rhsrelololo[i],
                   acc2hihihi,     acc2lohihi,     acc2hilohi,     acc2lolohi,
                   acc2hihilo,     acc2lohilo,     acc2hilolo,     acc2lololo);
         odf_inc(&rhsimhihihi[i],&rhsimlohihi[i],
                 &rhsimhilohi[i],&rhsimlolohi[i],
                 &rhsimhihilo[i],&rhsimlohilo[i],
                 &rhsimhilolo[i],&rhsimlololo[i],
                   acc3hihihi,     acc3lohihi,     acc3hilohi,     acc3lolohi,
                   acc3hihilo,     acc3lohilo,     acc3hilolo,     acc3lololo);
         odf_inc(&rhsimhihihi[i],&rhsimlohihi[i],
                 &rhsimhilohi[i],&rhsimlolohi[i],
                 &rhsimhihilo[i],&rhsimlohilo[i],
                 &rhsimhilolo[i],&rhsimlololo[i],
                   acc4hihihi,     acc4lohihi,     acc4hilohi,     acc4lolohi,
                   acc4hihilo,     acc4lohilo,     acc4hilolo,     acc4lololo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehihihi[i] << "  " << rhsrelohihi[i] << endl
              << "         "
              << rhsrehilohi[i] << "  " << rhsrelolohi[i] << endl
              << "         "
              << rhsrehihilo[i] << "  " << rhsrelohilo[i] << endl
              << "         "
              << rhsrehilolo[i] << "  " << rhsrelololo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhihihi[i] << "  " << rhsimlohihi[i] << endl
              << "         "
              << rhsimhilohi[i] << "  " << rhsimlolohi[i] << endl
              << "         "
              << rhsimhihilo[i] << "  " << rhsimlohilo[i] << endl
              << "         "
              << rhsimhilolo[i] << "  " << rhsimlololo[i] << endl;
      }
   }
   double **invArehihihi_h = new double*[dim];
   double **invArelohihi_h = new double*[dim];
   double **invArehilohi_h = new double*[dim];
   double **invArelolohi_h = new double*[dim];
   double **invArehihilo_h = new double*[dim];
   double **invArelohilo_h = new double*[dim];
   double **invArehilolo_h = new double*[dim];
   double **invArelololo_h = new double*[dim];
   double **invAimhihihi_h = new double*[dim];
   double **invAimlohihi_h = new double*[dim];
   double **invAimhilohi_h = new double*[dim];
   double **invAimlolohi_h = new double*[dim];
   double **invAimhihilo_h = new double*[dim];
   double **invAimlohilo_h = new double*[dim];
   double **invAimhilolo_h = new double*[dim];
   double **invAimlololo_h = new double*[dim];
 
   for(int i=0; i<dim; i++)
   {
      invArehihihi_h[i] = new double[dim];
      invArelohihi_h[i] = new double[dim];
      invArehilohi_h[i] = new double[dim];
      invArelolohi_h[i] = new double[dim];
      invArehihilo_h[i] = new double[dim];
      invArelohilo_h[i] = new double[dim];
      invArehilolo_h[i] = new double[dim];
      invArelololo_h[i] = new double[dim];
      invAimhihihi_h[i] = new double[dim];
      invAimlohihi_h[i] = new double[dim];
      invAimhilohi_h[i] = new double[dim];
      invAimlolohi_h[i] = new double[dim];
      invAimhihilo_h[i] = new double[dim];
      invAimlohilo_h[i] = new double[dim];
      invAimhilolo_h[i] = new double[dim];
      invAimlololo_h[i] = new double[dim];
   }
   double timelapsed_h,timelapsed_d,elapsedms;

   cout << "-> CPU computes the inverse ..." << endl;

   CPU_cmplx8_upper_inverse
      (dim,   Arehihihi,     Arelohihi,     Arehilohi,     Arelolohi,
              Arehihilo,     Arelohilo,     Arehilolo,     Arelololo, 
              Aimhihihi,     Aimlohihi,     Aimhilohi,     Aimlolohi,
              Aimhihilo,     Aimlohilo,     Aimhilolo,     Aimlololo,
           invArehihihi_h,invArelohihi_h,invArehilohi_h,invArelolohi_h,
           invArehihilo_h,invArelohilo_h,invArehilolo_h,invArelololo_h,
           invAimhihihi_h,invAimlohihi_h,invAimhilohi_h,invAimlolohi_h,
           invAimhihilo_h,invAimlohilo_h,invAimhilolo_h,invAimlololo_h,
       &timelapsed_h);

   if(verbose > 0)
   {
      cout << "The CPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "invA_h[" << i << "][" << j << "]re : "
                 << invArehihihi_h[i][j] << "  "
                 << invArelohihi_h[i][j] << endl
                 << "                 "
                 << invArehilohi_h[i][j] << "  "
                 << invArelolohi_h[i][j] << endl
                 << "                 "
                 << invArehihilo_h[i][j] << "  "
                 << invArelohilo_h[i][j] << endl
                 << "                 "
                 << invArehilolo_h[i][j] << "  "
                 << invArelololo_h[i][j] << endl;
            cout << "invA_h[" << i << "][" << j << "]im : "
                 << invAimhihihi_h[i][j] << "  "
                 << invAimlohihi_h[i][j] << endl
                 << "                 "
                 << invAimhilohi_h[i][j] << "  "
                 << invAimlolohi_h[i][j] << endl
                 << "                 "
                 << invAimhihilo_h[i][j] << "  "
                 << invAimlohilo_h[i][j] << endl
                 << "                 "
                 << invAimhilolo_h[i][j] << "  "
                 << invAimlololo_h[i][j] << endl;
         }
   }
/*
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
 */
   double *xrehihihi = new double[dim];
   double *xrelohihi = new double[dim];
   double *xrehilohi = new double[dim];
   double *xrelolohi = new double[dim];
   double *xrehihilo = new double[dim];
   double *xrelohilo = new double[dim];
   double *xrehilolo = new double[dim];
   double *xrelololo = new double[dim];
   double *ximhihihi = new double[dim];
   double *ximlohihi = new double[dim];
   double *ximhilohi = new double[dim];
   double *ximlolohi = new double[dim];
   double *ximhihilo = new double[dim];
   double *ximlohilo = new double[dim];
   double *ximhilolo = new double[dim];
   double *ximlololo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      xrehihihi[i] = 0.0; xrelohihi[i] = 0.0;
      xrehilohi[i] = 0.0; xrelolohi[i] = 0.0;
      xrehihilo[i] = 0.0; xrelohilo[i] = 0.0;
      xrehilolo[i] = 0.0; xrelololo[i] = 0.0;
      ximhihihi[i] = 0.0; ximlohihi[i] = 0.0;
      ximhilohi[i] = 0.0; ximlolohi[i] = 0.0;
      ximhihilo[i] = 0.0; ximlohilo[i] = 0.0;
      ximhilolo[i] = 0.0; ximlololo[i] = 0.0;

      for(int j=0; j<dim; j++) // x[i] = x[i] + invA_h[i][j]*rhs[j];
      {
         odf_mul(invArehihihi_h[i][j],invArelohihi_h[i][j],
                 invArehilohi_h[i][j],invArelolohi_h[i][j],
                 invArehihilo_h[i][j],invArelohilo_h[i][j],
                 invArehilolo_h[i][j],invArelololo_h[i][j],
                 rhsrehihihi[j],rhsrelohihi[j],rhsrehilohi[j],rhsrelolohi[j],
                 rhsrehihilo[j],rhsrelohilo[j],rhsrehilolo[j],rhsrelololo[j],
                 &acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                 &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo);
         odf_mul(invAimhihihi_h[i][j],invAimlohihi_h[i][j],
                 invAimhilohi_h[i][j],invAimlolohi_h[i][j],
                 invAimhihilo_h[i][j],invAimlohilo_h[i][j],
                 invAimhilolo_h[i][j],invAimlololo_h[i][j],
                 rhsimhihihi[j],rhsimlohihi[j],rhsimhilohi[j],rhsimlolohi[j],
                 rhsimhihilo[j],rhsimlohilo[j],rhsimhilolo[j],rhsimlololo[j],
                 &acc2hihihi,&acc2lohihi,&acc2hilohi,&acc2lolohi,
                 &acc2hihilo,&acc2lohilo,&acc2hilolo,&acc2lololo);
         odf_mul(invAimhihihi_h[i][j],invAimlohihi_h[i][j],
                 invAimhilohi_h[i][j],invAimlolohi_h[i][j],
                 invAimhihilo_h[i][j],invAimlohilo_h[i][j],
                 invAimhilolo_h[i][j],invAimlololo_h[i][j],
                 rhsrehihihi[j],rhsrelohihi[j],rhsrehilohi[j],rhsrelolohi[j],
                 rhsrehihilo[j],rhsrelohilo[j],rhsrehilolo[j],rhsrelololo[j],
                 &acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                 &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo);
         odf_mul(invArehihihi_h[i][j],invArelohihi_h[i][j],
                 invArehilolo_h[i][j],invArelololo_h[i][j],
                 invArehihihi_h[i][j],invArelohihi_h[i][j],
                 invArehilolo_h[i][j],invArelololo_h[i][j],
                 rhsimhihihi[j],rhsimlohihi[j],rhsimhilohi[j],rhsimlolohi[j],
                 rhsimhihilo[j],rhsimlohilo[j],rhsimhilolo[j],rhsimlololo[j],
                 &acc4hihihi,&acc4lohihi,&acc4hilohi,&acc4lolohi,
                 &acc4hihilo,&acc4lohilo,&acc4hilolo,&acc4lololo);
         odf_inc(&xrehihihi[i],&xrelohihi[i],&xrehilohi[i],&xrelolohi[i],
                 &xrehihilo[i],&xrelohilo[i],&xrehilolo[i],&xrelololo[i],
                 acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi,
                 acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo);
         odf_dec(&xrehihihi[i],&xrelohihi[i],&xrehilohi[i],&xrelolohi[i],
                 &xrehihilo[i],&xrelohilo[i],&xrehilolo[i],&xrelololo[i],
                 acc2hihihi,acc2lohihi,acc2hilohi,acc2lolohi,
                 acc2hihilo,acc2lohilo,acc2hilolo,acc2lololo);
         odf_inc(&ximhihihi[i],&ximlohihi[i],&ximhilohi[i],&ximlolohi[i],
                 &ximhihilo[i],&ximlohilo[i],&ximhilolo[i],&ximlololo[i],
                 acc3hihihi,acc3lohihi,acc3hilohi,acc3lolohi,
                 acc3hihilo,acc3lohilo,acc3hilolo,acc3lololo);
         odf_inc(&ximhihihi[i],&ximlohihi[i],&ximhilohi[i],&ximlolohi[i],
                 &ximhihilo[i],&ximlohilo[i],&ximhilolo[i],&ximlololo[i],
                 acc4hihihi,acc4lohihi,acc4hilohi,acc4lolohi,
                 acc4hihilo,acc4lohilo,acc4hilolo,acc4lololo);
      }
   }
   if(verbose > 0)
   {
      cout << "The solution computed with the CPU inverse :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehihihi[i] << "  " << xrelohihi[i] << endl
              << "         "
              << xrehilohi[i] << "  " << xrelolohi[i] << endl
              << "         "
              << xrehihilo[i] << "  " << xrelohilo[i] << endl
              << "         "
              << xrehilolo[i] << "  " << xrelololo[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhihihi[i] << "  " << ximlohihi[i] << endl
              << "         "
              << ximhilohi[i] << "  " << ximlolohi[i] << endl
              << "         "
              << ximhihilo[i] << "  " << ximlohilo[i] << endl
              << "         "
              << ximhilolo[i] << "  " << ximlololo[i] << endl;
      }
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on solution : "
        << cmplx8_Difference_Sum
              (dim,solrehihihi,solrelohihi,solrehilohi,solrelolohi,
                   solrehihilo,solrelohilo,solrehilolo,solrelololo,
                   solimhihihi,solimlohihi,solimhilohi,solimlolohi,
                   solimhihilo,solimlohilo,solimhilolo,solimlololo,
                     xrehihihi,  xrelohihi,  xrehilohi,  xrelolohi,
                     xrehihilo,  xrelohilo,  xrehilolo,  xrelololo,
                     ximhihihi,  ximlohihi, ximhilohi,  ximlolohi,
                     ximhihilo,  ximlohilo, ximhilolo,  ximlololo) << endl;
   cout << "Condition number : "
        << cmplx8_condition(dim, 
                 Arehihihi,     Arelohihi,     Arehilohi,     Arelolohi, 
                 Arehihilo,     Arelohilo,     Arehilolo,     Arelololo, 
                 Aimhihihi,     Aimlohihi,     Aimhilohi,     Aimlolohi,
                 Aimhihilo,     Aimlohilo,     Aimhilolo,     Aimlololo,
              invArehihihi_h,invArelohihi_h,invArehilohi_h,invArelolohi_h,
              invArehihilo_h,invArelohilo_h,invArehilolo_h,invArelololo_h,
              invAimhihihi_h,invAimlohihi_h,invAimhilohi_h,invAimlolohi_h,
              invAimhihilo_h,invAimlohilo_h,invAimhilolo_h,invAimlololo_h)
        << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
   cout << "                     Time spent by the kernel : ";
 /*
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;

 */
}

void test_real8_upper_tiling ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   const int dim = sizetile*numtiles;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Ahihihi = new double*[dim];
   double **Alohihi = new double*[dim];
   double **Ahilohi = new double*[dim];
   double **Alolohi = new double*[dim];
   double **Ahihilo = new double*[dim];
   double **Alohilo = new double*[dim];
   double **Ahilolo = new double*[dim];
   double **Alololo = new double*[dim];

   for(int i=0; i<dim; i++) 
   {
      Ahihihi[i] = new double[dim];
      Alohihi[i] = new double[dim];
      Ahilohi[i] = new double[dim];
      Alolohi[i] = new double[dim];
      Ahihilo[i] = new double[dim];
      Alohilo[i] = new double[dim];
      Ahilolo[i] = new double[dim];
      Alololo[i] = new double[dim];
   }
   random_dbl8_upper_matrix
      (dim,dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
               Ahihilo,Alohilo,Ahilolo,Alololo);
   // dbl4_random_upper_factor(dim,Ahihi,Alohi,Ahilo,Alolo);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihihi[i][j] << "  " << Alohihi[i][j] << endl
                 << "          "
                 << Ahilohi[i][j] << "  " << Alolohi[i][j] << endl
                 << "          "
                 << Ahihilo[i][j] << "  " << Alohilo[i][j] << endl
                 << "          "
                 << Ahilolo[i][j] << "  " << Alololo[i][j] << endl;
   }
   double *solhihihi = new double[dim];
   double *sollohihi = new double[dim];
   double *solhilohi = new double[dim];
   double *sollolohi = new double[dim];
   double *solhihilo = new double[dim];
   double *sollohilo = new double[dim];
   double *solhilolo = new double[dim];
   double *sollololo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solhihihi[i] = 1.0; sollohihi[i] = 0.0;
      solhilohi[i] = 0.0; sollolohi[i] = 0.0;
      solhihilo[i] = 0.0; sollohilo[i] = 0.0;
      solhilolo[i] = 0.0; sollololo[i] = 0.0;
   }
   double *xhihihi = new double[dim];
   double *xlohihi = new double[dim];
   double *xhilohi = new double[dim];
   double *xlolohi = new double[dim];
   double *xhihilo = new double[dim];
   double *xlohilo = new double[dim];
   double *xhilolo = new double[dim];
   double *xlololo = new double[dim];
   double *rhshihihi = new double[dim];
   double *rhslohihi = new double[dim];
   double *rhshilohi = new double[dim];
   double *rhslolohi = new double[dim];
   double *rhshihilo = new double[dim];
   double *rhslohilo = new double[dim];
   double *rhshilolo = new double[dim];
   double *rhslololo = new double[dim];
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<dim; i++)
   {
      rhshihihi[i] = 0.0; rhslohihi[i] = 0.0;
      rhshilohi[i] = 0.0; rhslolohi[i] = 0.0;
      rhshihilo[i] = 0.0; rhslohilo[i] = 0.0;
      rhshilolo[i] = 0.0; rhslololo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         odf_mul(Ahihihi[i][j],Alohihi[i][j],Ahilohi[i][j],Alolohi[i][j],
                 Ahihilo[i][j],Alohilo[i][j],Ahilolo[i][j],Alololo[i][j],
               solhihihi[j], sollohihi[j], solhilohi[j], sollolohi[j],
               solhihilo[j], sollohilo[j], solhilolo[j], sollololo[j],
              &acchihihi,   &acclohihi,   &acchilohi,  &acclolohi,
              &acchihilo,   &acclohilo,   &acchilolo,  &acclololo);
         odf_inc(&rhshihihi[i],&rhslohihi[i],&rhshilohi[i],&rhslolohi[i],
                 &rhshihilo[i],&rhslohilo[i],&rhshilolo[i],&rhslololo[i],
                  acchihihi,    acclohihi,    acchilohi,   acclolohi,
                  acchihilo,    acclohilo,    acchilolo,   acclololo);
      }
   }
   double timelapsed_h;

   cout << "-> CPU solves an upper triangular system ..." << endl;

   CPU_dbl8_upper_tiled_solver
      (dim,sizetile,numtiles,
         Ahihihi,  Alohihi,  Ahilohi,  Alolohi,
         Ahihilo,  Alohilo,  Ahilolo,  Alololo,
       rhshihihi,rhslohihi,rhshilohi,rhslolohi,
       rhshihilo,rhslohilo,rhshilolo,rhslololo,
         xhihihi,  xlohihi,  xhilohi,  xlolohi,
         xhihilo,  xlohilo,  xhilolo,  xlololo,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihihi[i][j] << "  " << Alohihi[i][j] << endl
                 << "          "
                 << Ahilohi[i][j] << "  " << Alolohi[i][j] << endl
                 << "          "
                 << Ahihilo[i][j] << "  " << Alohilo[i][j] << endl
                 << "          "
                 << Ahilolo[i][j] << "  " << Alololo[i][j] << endl;
   }
/*
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
   double timelapsed_d,elapsedms;
   double invlapsed,mullapsed,sublapsed;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;

   cout << "-> GPU solves an upper triangular system ..." << endl;

   GPU_dbl4_upper_tiled_solver
      (dim,sizetile,numtiles,
         Ahihi_d,  Alohi_d,  Ahilo_d,  Alolo_d,
       rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
         xhihi_d,  xlohi_d,  xhilo_d,  xlolo_d,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&timelapsed_d,
       &addcnt,&mulcnt,&divcnt);

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
*/
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "CPU solution computed with tiling :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhihihi[i] << "  " << xlohihi[i] << endl
              << "       "
              << xhilohi[i] << "  " << xlolohi[i] << endl
              << "       "
              << xhihilo[i] << "  " << xlohilo[i] << endl
              << "       "
              << xhilolo[i] << "  " << xlololo[i] << endl;
    /*
      cout << "GPU solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhihi_d[i] << "  " << xlohi_d[i] << endl
              << "       "
              << xhilo_d[i] << "  " << xlolo_d[i] << endl;
     */
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of CPU errors on solution : "
        << dbl8_Difference_Sum
              (dim,solhihihi,sollohihi,solhilohi,sollolohi,
                   solhihilo,sollohilo,solhilolo,sollololo,
                     xhihihi,  xlohihi,  xhilohi,  xlolohi,
                     xhihilo,  xlohilo,  xhilolo,  xlololo)
        << endl;
 /*
   cout << "   Sum of GPU errors on solution : "
        << dbl4_Difference_Sum(dim,solhihi,sollohi,solhilo,sollolo,
                                     xhihi_d,xlohi_d,xhilo_d,xlolo_d)
        << endl;
  */
   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
/*
   cout << "          Time spent to invert diagonal tiles : ";
   cout << invlapsed << " milliseconds." << endl;
   cout << "   Time spent to multiply with inverted tiles : ";
   cout << mullapsed << " milliseconds." << endl;
   cout << "             Time spent for back substitution : ";
   cout << sublapsed << " milliseconds." << endl;
   cout << "                    Time spent by all kernels : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
   cout << endl;
   cout << "             Number of additions/subtractions : "
        << addcnt << " x 89 " << endl;
   cout << "                    Number of multiplications : "
        << mulcnt << " x 336 " << endl;
   cout << "                          Number of divisions : "
        << divcnt << " x 893 " << endl;
   long long int flopcnt = 89*addcnt + 336*mulcnt + 893*divcnt;
   cout << "    Total number of floating-point operations : "
        << flopcnt << endl;
   cout << endl;
   double kernflops = 1000.0*((double) flopcnt)/elapsedms;
   double wallflops = ((double) flopcnt)/timelapsed_d;
   const int gigacnt = pow(2.0,30);
   cout << "Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
   cout << " Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;
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
 */
}

void test_cmplx8_upper_tiling ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   const int dim = sizetile*numtiles;

   cout << "-> generating a random upper triangular matrix of dimension "
        << dim << " ..." << endl;

   double **Arehihihi = new double*[dim];
   double **Arelohihi = new double*[dim];
   double **Arehilohi = new double*[dim];
   double **Arelolohi = new double*[dim];
   double **Arehihilo = new double*[dim];
   double **Arelohilo = new double*[dim];
   double **Arehilolo = new double*[dim];
   double **Arelololo = new double*[dim];
   double **Aimhihihi = new double*[dim];
   double **Aimlohihi = new double*[dim];
   double **Aimhilohi = new double*[dim];
   double **Aimlolohi = new double*[dim];
   double **Aimhihilo = new double*[dim];
   double **Aimlohilo = new double*[dim];
   double **Aimhilolo = new double*[dim];
   double **Aimlololo = new double*[dim];

   for(int i=0; i<dim; i++) 
   {
      Arehihihi[i] = new double[dim];
      Arelohihi[i] = new double[dim];
      Arehilohi[i] = new double[dim];
      Arelolohi[i] = new double[dim];
      Arehihilo[i] = new double[dim];
      Arelohilo[i] = new double[dim];
      Arehilolo[i] = new double[dim];
      Arelololo[i] = new double[dim];
      Aimhihihi[i] = new double[dim];
      Aimlohihi[i] = new double[dim];
      Aimhilohi[i] = new double[dim];
      Aimlolohi[i] = new double[dim];
      Aimhihilo[i] = new double[dim];
      Aimlohilo[i] = new double[dim];
      Aimhilolo[i] = new double[dim];
      Aimlololo[i] = new double[dim];
   }

   random_cmplx8_upper_matrix
      (dim,dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
               Arehihilo,Arelohilo,Arehilolo,Arelololo,
               Aimhihihi,Aimlohihi,Aimhilolo,Aimlolohi,
               Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);
   // cmplx4_random_upper_factor
   //   (dim,Arehihi,Arelohi,Arehilo,Arelolo,Aimhihi,Aimlohi,Aimhilo,Aimlolo);

   cout << scientific << setprecision(16);

   if(verbose > 0)
   {
      cout << "A random upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihihi[i][j] << "  " << Arelohihi[i][j] << endl
                 << "            "
                 << Arehilohi[i][j] << "  " << Arelolohi[i][j] << endl
                 << "            "
                 << Arehihilo[i][j] << "  " << Arelohilo[i][j] << endl
                 << "            "
                 << Arehilolo[i][j] << "  " << Arelololo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihihi[i][j] << "  " << Aimlohihi[i][j] << endl
                 << "            "
                 << Aimhilohi[i][j] << "  " << Aimlolohi[i][j] << endl
                 << "            "
                 << Aimhihilo[i][j] << "  " << Aimlohilo[i][j] << endl
                 << "            "
                 << Aimhilolo[i][j] << "  " << Aimlololo[i][j] << endl;
         }
   }
   double *solrehihihi = new double[dim];
   double *solrelohihi = new double[dim];
   double *solrehilohi = new double[dim];
   double *solrelolohi = new double[dim];
   double *solrehihilo = new double[dim];
   double *solrelohilo = new double[dim];
   double *solrehilolo = new double[dim];
   double *solrelololo = new double[dim];
   double *solimhihihi = new double[dim];
   double *solimlohihi = new double[dim];
   double *solimhilohi = new double[dim];
   double *solimlolohi = new double[dim];
   double *solimhihilo = new double[dim];
   double *solimlohilo = new double[dim];
   double *solimhilolo = new double[dim];
   double *solimlololo = new double[dim];

   for(int i=0; i<dim; i++)
   {
      solrehihihi[i] = 1.0; solrelohihi[i] = 0.0;
      solrehilohi[i] = 0.0; solrelolohi[i] = 0.0;
      solrehihilo[i] = 0.0; solrelohilo[i] = 0.0;
      solrehilolo[i] = 0.0; solrelololo[i] = 0.0;
      solimhihihi[i] = 0.0; solimlohihi[i] = 0.0;
      solimhilohi[i] = 0.0; solimlolohi[i] = 0.0;
      solimhihilo[i] = 0.0; solimlohilo[i] = 0.0;
      solimhilolo[i] = 0.0; solimlololo[i] = 0.0;
   }
   double *xrehihihi = new double[dim];
   double *xrelohihi = new double[dim];
   double *xrehilohi = new double[dim];
   double *xrelolohi = new double[dim];
   double *xrehihilo = new double[dim];
   double *xrelohilo = new double[dim];
   double *xrehilolo = new double[dim];
   double *xrelololo = new double[dim];
   double *ximhihihi = new double[dim];
   double *ximlohihi = new double[dim];
   double *ximhilohi = new double[dim];
   double *ximlolohi = new double[dim];
   double *ximhihilo = new double[dim];
   double *ximlohilo = new double[dim];
   double *ximhilolo = new double[dim];
   double *ximlololo = new double[dim];
   double *rhsrehihihi = new double[dim];
   double *rhsrelohihi = new double[dim];
   double *rhsrehilohi = new double[dim];
   double *rhsrelolohi = new double[dim];
   double *rhsrehihilo = new double[dim];
   double *rhsrelohilo = new double[dim];
   double *rhsrehilolo = new double[dim];
   double *rhsrelololo = new double[dim];
   double *rhsimhihihi = new double[dim];
   double *rhsimlohihi = new double[dim];
   double *rhsimhilohi = new double[dim];
   double *rhsimlolohi = new double[dim];
   double *rhsimhihilo = new double[dim];
   double *rhsimlohilo = new double[dim];
   double *rhsimhilolo = new double[dim];
   double *rhsimlololo = new double[dim];
   double acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi;
   double acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo;
   double acc2hihihi,acc2lohihi,acc2hilohi,acc2lolohi;
   double acc2hihilo,acc2lohilo,acc2hilolo,acc2lololo;
   double acc3hihihi,acc3lohihi,acc3hilohi,acc3lolohi;
   double acc3hihilo,acc3lohilo,acc3hilolo,acc3lololo;
   double acc4hihihi,acc4lohihi,acc4hilohi,acc4lolohi;
   double acc4hihilo,acc4lohilo,acc4hilolo,acc4lololo;

   for(int i=0; i<dim; i++)
   {
      rhsrehihihi[i] = 0.0; rhsrelohihi[i] = 0.0;
      rhsrehilohi[i] = 0.0; rhsrelolohi[i] = 0.0;
      rhsrehihilo[i] = 0.0; rhsrelohilo[i] = 0.0;
      rhsrehilolo[i] = 0.0; rhsrelololo[i] = 0.0;
      rhsimhihihi[i] = 0.0; rhsimlohihi[i] = 0.0;
      rhsimhilohi[i] = 0.0; rhsimlolohi[i] = 0.0;
      rhsimhihilo[i] = 0.0; rhsimlohilo[i] = 0.0;
      rhsimhilolo[i] = 0.0; rhsimlololo[i] = 0.0;

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         odf_mul(Arehihihi[i][j],Arelohihi[i][j],
                 Arehilohi[i][j],Arelolohi[i][j],
                 Arehihilo[i][j],Arelohilo[i][j],
                 Arehilolo[i][j],Arelololo[i][j],
               solrehihihi[j] ,solrelohihi[j], solrehilohi[j], solrelolohi[j],
               solrehihilo[j] ,solrelohilo[j], solrehilolo[j], solrelololo[j],
               &acc1hihihi,    &acc1lohihi,    &acc1hilohi,    &acc1lolohi,
               &acc1hihilo,    &acc1lohilo,    &acc1hilolo,    &acc1lololo);
         odf_mul(Aimhihihi[i][j],Aimlohihi[i][j],
                 Aimhilohi[i][j],Aimlolohi[i][j],
                 Aimhihilo[i][j],Aimlohilo[i][j],
                 Aimhilolo[i][j],Aimlololo[i][j],
               solimhihihi[j], solimlohihi[j], solimhilohi[j], solimlolohi[j],
               solimhihilo[j], solimlohilo[j], solimhilolo[j], solimlololo[j],
               &acc2hihihi,    &acc2lohihi,    &acc2hilohi,   &acc2lolohi,
               &acc2hihilo,    &acc2lohilo,    &acc2hilolo,   &acc2lololo);
         odf_mul(Aimhihihi[i][j],Aimlohihi[i][j],
                 Aimhilohi[i][j],Aimlolohi[i][j],
                 Aimhihilo[i][j],Aimlohilo[i][j],
                 Aimhilolo[i][j],Aimlololo[i][j],
               solrehihihi[j], solrelohihi[j], solrehilohi[j], solrelolohi[j],
               solrehihilo[j], solrelohilo[j], solrehilolo[j], solrelololo[j],
               &acc3hihihi,    &acc3lohihi,    &acc3hilohi,    &acc3lolohi,
               &acc3hihilo,    &acc3lohilo,    &acc3hilolo,    &acc3lololo);
         odf_mul(Arehihihi[i][j],Arelohihi[i][j],
                 Arehilohi[i][j],Arelolohi[i][j],
                 Arehihilo[i][j],Arelohilo[i][j],
                 Arehilolo[i][j],Arelololo[i][j],
               solimhihihi[j], solimlohihi[j], solimhilohi[j], solimlolohi[j],
               solimhihilo[j], solimlohilo[j], solimhilolo[j], solimlololo[j],
               &acc4hihihi,    &acc4lohihi,    &acc4hilohi,    &acc4lolohi,
               &acc4hihilo,    &acc4lohilo,    &acc4hilolo,    &acc4lololo);
         odf_dec(&acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                 &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo,
                  acc2hihihi, acc2lohihi, acc2hilohi, acc2lolohi,
                  acc2hihilo, acc2lohilo, acc2hilolo, acc2lololo);
         odf_inc(&rhsrehihihi[i],&rhsrelohihi[i],
                 &rhsrehilohi[i],&rhsrelolohi[i],
                 &rhsrehihilo[i],&rhsrelohilo[i],
                 &rhsrehilolo[i],&rhsrelololo[i],
                   acc1hihihi,     acc1lohihi,    acc1hilohi,    acc1lolohi,
                   acc1hihilo,     acc1lohilo,    acc1hilolo,    acc1lololo);
         odf_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                 &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
                  acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
                  acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
         odf_inc(&rhsimhihihi[i],&rhsimlohihi[i],
                 &rhsimhilohi[i],&rhsimlolohi[i],
                 &rhsimhihilo[i],&rhsimlohilo[i],
                 &rhsimhilolo[i],&rhsimlololo[i],
                   acc3hihihi,     acc3lohihi,    acc3hilohi,    acc3lolohi,
                   acc3hihilo,     acc3lohilo,    acc3hilolo,    acc3lololo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehihihi[i] << "  " << rhsrelohihi[i] << endl
              << "         "
              << rhsrehilohi[i] << "  " << rhsrelolohi[i] << endl
              << "         "
              << rhsrehihilo[i] << "  " << rhsrelohilo[i] << endl
              << "         "
              << rhsrehilolo[i] << "  " << rhsrelololo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhihihi[i] << "  " << rhsimlohihi[i] << endl
              << "         "
              << rhsimhilohi[i] << "  " << rhsimlolohi[i] << endl
              << "         "
              << rhsimhihilo[i] << "  " << rhsimlohilo[i] << endl
              << "         "
              << rhsimhilolo[i] << "  " << rhsimlololo[i] << endl;
      }
   }
   double timelapsed_h;

   cout << "-> CPU solves an upper triangular system ..." << endl;

   CPU_cmplx8_upper_tiled_solver
      (dim,sizetile,numtiles,
       Arehihihi,Arelohihi,Arehilohi,Arelolohi,
       Arehihilo,Arelohilo,Arehilolo,Arelololo,
       Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
       Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo,
       rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
       rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
       rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
       rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
       xrehihihi,xrelohihi,xrehilohi,xrelolohi,
       xrehihilo,xrelohilo,xrehilolo,xrelololo,
       ximhihihi,ximlohihi,ximhilohi,ximlolohi,
       ximhihilo,ximlohilo,ximhilolo,ximlololo,&timelapsed_h);

   if(verbose > 0)
   {
      cout << "The matrix computed by the host :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihihi[i][j] << "  " << Arelohihi[i][j] << endl
                 << "            "
                 << Arehilohi[i][j] << "  " << Arelolohi[i][j] << endl
                 << "            "
                 << Arehihilo[i][j] << "  " << Arelohilo[i][j] << endl
                 << "            "
                 << Arehilolo[i][j] << "  " << Arelololo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihihi[i][j] << "  " << Aimlohihi[i][j] << endl
                 << "            "
                 << Aimhilohi[i][j] << "  " << Aimlolohi[i][j] << endl
                 << "            "
                 << Aimhihilo[i][j] << "  " << Aimlohilo[i][j] << endl
                 << "            "
                 << Aimhilolo[i][j] << "  " << Aimlololo[i][j] << endl;
         }
   }
/*
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
   double timelapsed_d,elapsedms;
   double invlapsed,mullapsed,sublapsed;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;

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
       &addcnt,&mulcnt,&divcnt);

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

 */

   if(verbose > 0)
   {
      cout << "CPU solution computed with tiling :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehihihi[i] << "  " << xrelohihi[i] << endl
              << "         "
              << xrehilohi[i] << "  " << xrelolohi[i] << endl
              << "         "
              << xrehihilo[i] << "  " << xrelohilo[i] << endl
              << "         "
              << xrehilolo[i] << "  " << xrelololo[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhihihi[i] << "  " << ximlohihi[i] << endl
              << "         "
              << ximhilohi[i] << "  " << ximlolohi[i] << endl
              << "         "
              << ximhihilo[i] << "  " << ximlohilo[i] << endl
              << "         "
              << ximhilolo[i] << "  " << ximlololo[i] << endl;
      }
    /*
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
    */
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of CPU errors on solution : "
        << cmplx8_Difference_Sum(dim,
              solrehihihi,solrelohihi,solrehilohi,solrelolohi,
              solrehihilo,solrelohilo,solrehilolo,solrelololo,
              solimhihihi,solimlohihi,solimhilohi,solimlolohi,
              solimhihilo,solimlohilo,solimhilolo,solimlololo,
                xrehihihi,  xrelohihi,  xrehilohi,  xrelolohi, 
                xrehihilo,  xrelohilo,  xrehilolo,  xrelololo, 
                ximhihihi,  ximlohihi,  ximhilohi,  ximlolohi,
                ximhihilo,  ximlohilo,  ximhilolo,  ximlololo)
        << endl;
/*
   cout << "   Sum of GPU errors on solution : "
        << cmplx4_Difference_Sum(dim,
               solrehihi,solrelohi,solrehilo,solrelolo,
               solimhihi,solimlohi,solimhilo,solimlolo,
                 xrehihi_d,xrelohi_d,xrehilo_d,xrelolo_d,
                 ximhihi_d,ximlohi_d,ximhilo_d,ximlolo_d)
        << endl;
 */
   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
/*
   cout << "          Time spent to invert diagonal tiles : ";
   cout << invlapsed << " milliseconds." << endl;
   cout << "   Time spent to multiply with inverted tiles : ";
   cout << mullapsed << " milliseconds." << endl;
   cout << "             Time spent for back substitution : ";
   cout << sublapsed << " milliseconds." << endl;
   cout << "                    Time spent by all kernels : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
   cout << endl;
   cout << "             Number of additions/subtractions : "
        << addcnt << " x 89 " << endl;
   cout << "                    Number of multiplications : "
        << mulcnt << " x 336 " << endl;
   cout << "                          Number of divisions : "
        << divcnt << " x 893 " << endl;
   long long int flopcnt = 89*addcnt + 336*mulcnt + 893*divcnt;
   cout << "    Total number of floating-point operations : "
        << flopcnt << endl;
   cout << endl;
   double kernflops = 1000.0*((double) flopcnt)/elapsedms;
   double wallflops = ((double) flopcnt)/timelapsed_d;
   const int gigacnt = pow(2.0,30);
   cout << "Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
   cout << " Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;

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
 */
}
