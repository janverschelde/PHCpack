// The file dbl8_tabs_testers.cpp defines the functions specified in
// the file dbl8_tabs_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "octo_double_functions.h"
#include "random8_matrices.h"
#include "dbl8_factorizations.h"
#include "dbl8_tabs_host.h"
#include "dbl8_tabs_kernels.h"
#include "dbl_test_utilities.h"
#include "dbl8_test_utilities.h"
#include "write_dbl8_bstimeflops.h"
#include "dbl_data_files.h"
#include "dbl_tabs_testers.h"

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
   double **invAhihihi_d = new double*[dim];
   double **invAlohihi_d = new double*[dim];
   double **invAhilohi_d = new double*[dim];
   double **invAlolohi_d = new double*[dim];
   double **invAhihilo_d = new double*[dim];
   double **invAlohilo_d = new double*[dim];
   double **invAhilolo_d = new double*[dim];
   double **invAlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      invAhihihi_d[i] = new double[dim];
      invAlohihi_d[i] = new double[dim];
      invAhilohi_d[i] = new double[dim];
      invAlolohi_d[i] = new double[dim];
      invAhihilo_d[i] = new double[dim];
      invAlohilo_d[i] = new double[dim];
      invAhilolo_d[i] = new double[dim];
      invAlololo_d[i] = new double[dim];
   }
   cout << "-> GPU computes the inverse ..." << endl;

   GPU_dbl8_upper_inverse
      (dim,Ahihihi,     Alohihi,     Ahilohi,     Alolohi,
           Ahihilo,     Alohilo,     Ahilolo,     Alololo,
        invAhihihi_d,invAlohihi_d,invAhilohi_d,invAlolohi_d,
        invAhihilo_d,invAlohilo_d,invAhilolo_d,invAlololo_d,
       &elapsedms,&timelapsed_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "invA_d[" << i << "][" << j << "] : "
                 << invAhihihi_d[i][j] << "  " << invAlohihi_d[i][j] << endl
                 << "               "
                 << invAhilohi_d[i][j] << "  " << invAlolohi_d[i][j] << endl
                 << "               "
                 << invAhihilo_d[i][j] << "  " << invAlohilo_d[i][j] << endl
                 << "               "
                 << invAhilolo_d[i][j] << "  " << invAlololo_d[i][j] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on inverse : "
        << dbl8_Matrix_Difference_Sum
              (dim,invAhihihi_h,invAlohihi_h,invAhilohi_h,invAlolohi_h,
                   invAhihilo_h,invAlohilo_h,invAhilolo_h,invAlololo_h,
                   invAhihihi_d,invAlohihi_d,invAhilohi_d,invAlolohi_d,
                   invAhihilo_d,invAlohilo_d,invAhilolo_d,invAlololo_d)
        << endl;

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
   cout << "                     Time spent by the kernel : ";
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
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
      for(int j=0; j<dim; j++)
      {
         Arelohihi[i][j] = 0.0;
         Arehilohi[i][j] = 0.0;
         Arelolohi[i][j] = 0.0;
         Arehihilo[i][j] = 0.0;
         Arelohilo[i][j] = 0.0;
         Arehilolo[i][j] = 0.0;
         Arelololo[i][j] = 0.0;
         Aimhihihi[i][j] = 0.0;
         Aimlohihi[i][j] = 0.0;
         Aimhilohi[i][j] = 0.0;
         Aimlolohi[i][j] = 0.0;
         Aimhihilo[i][j] = 0.0;
         Aimlohilo[i][j] = 0.0;
         Aimhilolo[i][j] = 0.0;
         Aimlololo[i][j] = 0.0;
      }
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
   // cmplx_random_upper_factor(dim,Arehihihi,Aimhihihi);
  /*
   dbl8_random_upper_factor
      (dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
           Arehihilo,Arelohilo,Arehilolo,Arelololo);
   */
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
   double **invArehihihi_d = new double*[dim];
   double **invArelohihi_d = new double*[dim];
   double **invArehilohi_d = new double*[dim];
   double **invArelolohi_d = new double*[dim];
   double **invArehihilo_d = new double*[dim];
   double **invArelohilo_d = new double*[dim];
   double **invArehilolo_d = new double*[dim];
   double **invArelololo_d = new double*[dim];
   double **invAimhihihi_d = new double*[dim];
   double **invAimlohihi_d = new double*[dim];
   double **invAimhilohi_d = new double*[dim];
   double **invAimlolohi_d = new double*[dim];
   double **invAimhihilo_d = new double*[dim];
   double **invAimlohilo_d = new double*[dim];
   double **invAimhilolo_d = new double*[dim];
   double **invAimlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      invArehihihi_d[i] = new double[dim];
      invArelohihi_d[i] = new double[dim];
      invArehilohi_d[i] = new double[dim];
      invArelolohi_d[i] = new double[dim];
      invArehihilo_d[i] = new double[dim];
      invArelohilo_d[i] = new double[dim];
      invArehilolo_d[i] = new double[dim];
      invArelololo_d[i] = new double[dim];
      invAimhihihi_d[i] = new double[dim];
      invAimlohihi_d[i] = new double[dim];
      invAimhilohi_d[i] = new double[dim];
      invAimlolohi_d[i] = new double[dim];
      invAimhihilo_d[i] = new double[dim];
      invAimlohilo_d[i] = new double[dim];
      invAimhilolo_d[i] = new double[dim];
      invAimlololo_d[i] = new double[dim];
   }

   cout << "-> GPU computes the inverse ..." << endl;

   GPU_cmplx8_upper_inverse
      (dim,   Arehihihi,     Arelohihi,     Arehilohi,     Arelolohi,     
              Arehihilo,     Arelohilo,     Arehilolo,     Arelololo,     
              Aimhihihi,     Aimlohihi,     Aimhilohi,     Aimlolohi,
              Aimhihilo,     Aimlohilo,     Aimhilolo,     Aimlololo,
           invArehihihi_d,invArelohihi_d,invArehilohi_d,invArelolohi_d,
           invArehihilo_d,invArelohilo_d,invArehilolo_d,invArelololo_d,
           invAimhihihi_d,invAimlohihi_d,invAimhilohi_d,invAimlolohi_d,
           invAimhihilo_d,invAimlohilo_d,invAimhilolo_d,invAimlololo_d,
       &elapsedms,&timelapsed_d);

   if(verbose > 0)
   {
      cout << "The GPU inverse of the upper triangular matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "invA_d[" << i << "][" << j << "]re : "
                 << invArehihihi_d[i][j] << "  "
                 << invArelohihi_d[i][j] << endl
                 << "                 "
                 << invArehilohi_d[i][j] << "  "
                 << invArelolohi_d[i][j] << endl
                 << "                 "
                 << invArehihilo_d[i][j] << "  "
                 << invArelohilo_d[i][j] << endl
                 << "                 "
                 << invArehilolo_d[i][j] << "  "
                 << invArelololo_d[i][j] << endl;
            cout << "invA_d[" << i << "][" << j << "]im : "
                 << invAimhihihi_d[i][j] << "  "
                 << invAimlohihi_d[i][j] << endl
                 << "                 "
                 << invAimhilohi_d[i][j] << "  "
                 << invAimlolohi_d[i][j] << endl
                 << "                 "
                 << invAimhihilo_d[i][j] << "  "
                 << invAimlohilo_d[i][j] << endl
                 << "                 "
                 << invAimhilolo_d[i][j] << "  "
                 << invAimlololo_d[i][j] << endl;
         }
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on inverse : "
        << cmplx8_Matrix_Difference_Sum
              (dim,invArehihihi_h,invArelohihi_h,invArehilohi_h,invArelolohi_h,
                   invArehihilo_h,invArelohilo_h,invArehilolo_h,invArelololo_h,
                   invAimhihihi_h,invAimlohihi_h,invAimhilohi_h,invAimlolohi_h,
                   invAimhihilo_h,invAimlohilo_h,invAimhilolo_h,invAimlololo_h,
                   invArehihihi_d,invArelohihi_d,invArehilohi_d,invArelolohi_d,
                   invArehihilo_d,invArelohilo_d,invArehilolo_d,invArelololo_d,
                   invAimhihihi_d,invAimlohihi_d,invAimhilohi_d,invAimlolohi_d,
                   invAimhihilo_d,invAimlohilo_d,invAimhilolo_d,invAimlololo_d)
        << endl;
 
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
                 invArehilohi_h[i][j],invArelolohi_h[i][j],
                 invArehihilo_h[i][j],invArelohilo_h[i][j],
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
   cout << elapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
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

   cout << "Generate a random matrix (1 = yes, 0 = read matrix) : ";
   int rndmat; cin >> rndmat;

   if(rndmat == 1)
      cout << "-> generating a random upper triangular matrix of dimension "
           << dim << " ..." << endl;
   else
      cout << "-> reading a random upper triangular matrix of dimension "
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
      for(int j=0; j<dim; j++)
      {
         Alohihi[i][j] = 0.0;
         Ahilohi[i][j] = 0.0;
         Alolohi[i][j] = 0.0;
         Ahihilo[i][j] = 0.0;
         Alohilo[i][j] = 0.0;
         Ahilolo[i][j] = 0.0;
         Alololo[i][j] = 0.0;
      }
   
   }
   // random_dbl8_upper_matrix
   //   (dim,dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
   //            Ahihilo,Alohilo,Ahilolo,Alololo);
   if(rndmat == 1)
      dbl_random_upper_factor(dim,Ahihihi);
   else
   {
      cout << "Give the name of a file : ";
      string filename; cin >> filename;
      cout << "-> reading " << dim*dim
           << " numbers from " << filename << " ..." << endl;
      dbl_read_matrix(filename,dim,Ahihihi);
   }
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
   double *xhihihi_d = new double[dim];
   double *xlohihi_d = new double[dim];
   double *xhilohi_d = new double[dim];
   double *xlolohi_d = new double[dim];
   double *xhihilo_d = new double[dim];
   double *xlohilo_d = new double[dim];
   double *xhilolo_d = new double[dim];
   double *xlololo_d = new double[dim];
   double *rhshihihi_d = new double[dim];
   double *rhslohihi_d = new double[dim];
   double *rhshilohi_d = new double[dim];
   double *rhslolohi_d = new double[dim];
   double *rhshihilo_d = new double[dim];
   double *rhslohilo_d = new double[dim];
   double *rhshilolo_d = new double[dim];
   double *rhslololo_d = new double[dim];
   double **Ahihihi_d = new double*[dim];
   double **Alohihi_d = new double*[dim];
   double **Ahilohi_d = new double*[dim];
   double **Alolohi_d = new double*[dim];
   double **Ahihilo_d = new double*[dim];
   double **Alohilo_d = new double*[dim];
   double **Ahilolo_d = new double*[dim];
   double **Alololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      rhshihihi_d[i] = rhshihihi[i];
      rhslohihi_d[i] = rhslohihi[i];
      rhshilohi_d[i] = rhshilohi[i];
      rhslolohi_d[i] = rhslolohi[i];
      rhshihilo_d[i] = rhshihilo[i];
      rhslohilo_d[i] = rhslohilo[i];
      rhshilolo_d[i] = rhshilolo[i];
      rhslololo_d[i] = rhslololo[i];

      Ahihihi_d[i] = new double[dim];
      Alohihi_d[i] = new double[dim];
      Ahilohi_d[i] = new double[dim];
      Alolohi_d[i] = new double[dim];
      Ahihilo_d[i] = new double[dim];
      Alohilo_d[i] = new double[dim];
      Ahilolo_d[i] = new double[dim];
      Alololo_d[i] = new double[dim];

      for(int j=0; j<dim; j++)
      {
         Ahihihi_d[i][j] = Ahihihi[i][j];
         Alohihi_d[i][j] = Alohihi[i][j];
         Ahilohi_d[i][j] = Ahilohi[i][j];
         Alolohi_d[i][j] = Alolohi[i][j];
         Ahihilo_d[i][j] = Ahihilo[i][j];
         Alohilo_d[i][j] = Alohilo[i][j];
         Ahilolo_d[i][j] = Ahilolo[i][j];
         Alololo_d[i][j] = Alololo[i][j];
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
   double timelapsed_d,elapsedms;
   double invlapsed,mullapsed,sublapsed;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;

   cout << "-> GPU solves an upper triangular system ..." << endl;

   GPU_dbl8_upper_tiled_solver
      (dim,sizetile,numtiles,
         Ahihihi_d,  Alohihi_d,  Ahilohi_d,  Alolohi_d,
         Ahihilo_d,  Alohilo_d,  Ahilolo_d,  Alololo_d,
       rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
       rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
         xhihihi_d,  xlohihi_d,  xhilohi_d,  xlolohi_d,
         xhihilo_d,  xlohilo_d,  xhilolo_d,  xlololo_d,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&timelapsed_d,
       &addcnt,&mulcnt,&divcnt);

   if(verbose > 0)
   {
      cout << "The matrix returned by the device :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihihi_d[i][j] << "  " << Alohihi_d[i][j] << endl
                 << "          "
                 << Ahilohi_d[i][j] << "  " << Alolohi_d[i][j] << endl
                 << "          "
                 << Ahihilo_d[i][j] << "  " << Alohilo_d[i][j] << endl
                 << "          "
                 << Ahilolo_d[i][j] << "  " << Alololo_d[i][j] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on diagonal tiles : "
        << dbl8_Diagonal_Difference_Sum
              (numtiles,sizetile,
               Ahihihi,  Alohihi,  Ahilohi,  Alolohi,
               Ahihilo,  Alohilo,  Ahilolo,  Alololo,
               Ahihihi_d,Alohihi_d,Ahilohi_d,Alolohi_d,
               Ahihilo_d,Alohilo_d,Ahilolo_d,Alololo_d)
        << endl;

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
      cout << "GPU solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : "
              << xhihihi_d[i] << "  " << xlohihi_d[i] << endl
              << "       "
              << xhilohi_d[i] << "  " << xlolohi_d[i] << endl
              << "       "
              << xhihilo_d[i] << "  " << xlohilo_d[i] << endl
              << "       "
              << xhilolo_d[i] << "  " << xlololo_d[i] << endl;
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of CPU errors on solution : "
        << dbl8_Difference_Sum
              (dim,solhihihi,sollohihi,solhilohi,sollolohi,
                   solhihilo,sollohilo,solhilolo,sollololo,
                     xhihihi,  xlohihi,  xhilohi,  xlolohi,
                     xhihilo,  xlohilo,  xhilolo,  xlololo)
        << endl;
   cout << "   Sum of GPU errors on solution : "
        << dbl8_Difference_Sum
              (dim,solhihihi,sollohihi,solhilohi,sollolohi,
                   solhihilo,sollohilo,solhilolo,sollololo,
                     xhihihi_d,xlohihi_d,xhilohi_d,xlolohi_d,
                     xhihilo_d,xlohilo_d,xhilolo_d,xlololo_d)
        << endl;
   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;

   write_dbl8_bstimeflops
     (sizetile,numtiles,0,invlapsed,mullapsed,sublapsed,elapsedms,
      timelapsed_d,addcnt,mulcnt,divcnt);

   for(int i=0; i<dim; i++)
   {
      free(Ahihihi[i]); free(Alohihi[i]); free(Ahilohi[i]); free(Alolohi[i]);
      free(Ahihilo[i]); free(Alohilo[i]); free(Ahilolo[i]); free(Alololo[i]);
      free(Ahihihi_d[i]); free(Alohihi_d[i]);
      free(Ahilohi_d[i]); free(Alolohi_d[i]);
      free(Ahihilo_d[i]); free(Alohilo_d[i]);
      free(Ahilolo_d[i]); free(Alololo_d[i]);
   }
   free(Ahihihi); free(Alohihi); free(Ahilohi); free(Alolohi);
   free(Ahihilo); free(Alohilo); free(Ahilolo); free(Alololo);
   free(Ahihihi_d); free(Alohihi_d); free(Ahilohi_d); free(Alolohi_d);
   free(Ahihilo_d); free(Alohilo_d); free(Ahilolo_d); free(Alololo_d);
   free(solhihihi); free(sollohihi); free(solhilohi); free(sollolohi);
   free(solhihilo); free(sollohilo); free(solhilolo); free(sollololo);
   free(rhshihihi); free(rhslohihi); free(rhshilohi); free(rhslolohi);
   free(rhshihilo); free(rhslohilo); free(rhshilolo); free(rhslololo);
   free(xhihihi); free(xlohihi); free(xhilohi); free(xlolohi); 
   free(xhihilo); free(xlohilo); free(xhilolo); free(xlololo); 
   free(rhshihihi_d); free(rhslohihi_d); free(rhshilohi_d); free(rhslolohi_d);
   free(rhshihilo_d); free(rhslohilo_d); free(rhshilolo_d); free(rhslololo_d);
   free(xhihihi_d); free(xlohihi_d); free(xhilohi_d); free(xlolohi_d);
   free(xhihilo_d); free(xlohilo_d); free(xhilolo_d); free(xlololo_d);
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

   cout << "Generate a random matrix (1 = yes, 0 = read matrix) : ";
   int rndmat; cin >> rndmat;

   if(rndmat == 1)
      cout << "-> generating a random upper triangular matrix of dimension "
           << dim << " ..." << endl;
   else
      cout << "-> reading a random upper triangular matrix of dimension "
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
      for(int j=0; j<dim; j++)
      {
         Arelohihi[i][j] = 0.0;
         Arehilohi[i][j] = 0.0;
         Arelolohi[i][j] = 0.0;
         Arehihilo[i][j] = 0.0;
         Arelohilo[i][j] = 0.0;
         Arehilolo[i][j] = 0.0;
         Arelololo[i][j] = 0.0;
         Aimlohihi[i][j] = 0.0;
         Aimhilohi[i][j] = 0.0;
         Aimlolohi[i][j] = 0.0;
         Aimhihilo[i][j] = 0.0;
         Aimlohilo[i][j] = 0.0;
         Aimhilolo[i][j] = 0.0;
         Aimlololo[i][j] = 0.0;
      }
   }
  /*
   random_cmplx8_upper_matrix
      (dim,dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
               Arehihilo,Arelohilo,Arehilolo,Arelololo,
               Aimhihihi,Aimlohihi,Aimhilolo,Aimlolohi,
               Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);
   */
   if(rndmat == 1)
      cmplx_random_upper_factor(dim,Arehihihi,Aimhihihi);
   else
   {
      cout << "Give the name of a file : ";
      string filename; cin >> filename;
      cout << "-> reading " << dim*dim
           << " numbers from " << filename << " ..." << endl;
      cmplx_read_matrix(filename,dim,Arehihihi,Aimhihihi);
   }
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
   double *xrehihihi_d = new double[dim];
   double *xrelohihi_d = new double[dim];
   double *xrehilohi_d = new double[dim];
   double *xrelolohi_d = new double[dim];
   double *xrehihilo_d = new double[dim];
   double *xrelohilo_d = new double[dim];
   double *xrehilolo_d = new double[dim];
   double *xrelololo_d = new double[dim];
   double *ximhihihi_d = new double[dim];
   double *ximlohihi_d = new double[dim];
   double *ximhilohi_d = new double[dim];
   double *ximlolohi_d = new double[dim];
   double *ximhihilo_d = new double[dim];
   double *ximlohilo_d = new double[dim];
   double *ximhilolo_d = new double[dim];
   double *ximlololo_d = new double[dim];
   double *rhsrehihihi_d = new double[dim];
   double *rhsrelohihi_d = new double[dim];
   double *rhsrehilohi_d = new double[dim];
   double *rhsrelolohi_d = new double[dim];
   double *rhsrehihilo_d = new double[dim];
   double *rhsrelohilo_d = new double[dim];
   double *rhsrehilolo_d = new double[dim];
   double *rhsrelololo_d = new double[dim];
   double *rhsimhihihi_d = new double[dim];
   double *rhsimlohihi_d = new double[dim];
   double *rhsimhilohi_d = new double[dim];
   double *rhsimlolohi_d = new double[dim];
   double *rhsimhihilo_d = new double[dim];
   double *rhsimlohilo_d = new double[dim];
   double *rhsimhilolo_d = new double[dim];
   double *rhsimlololo_d = new double[dim];
   double **Arehihihi_d = new double*[dim];
   double **Arelohihi_d = new double*[dim];
   double **Arehilohi_d = new double*[dim];
   double **Arelolohi_d = new double*[dim];
   double **Arehihilo_d = new double*[dim];
   double **Arelohilo_d = new double*[dim];
   double **Arehilolo_d = new double*[dim];
   double **Arelololo_d = new double*[dim];
   double **Aimhihihi_d = new double*[dim];
   double **Aimlohihi_d = new double*[dim];
   double **Aimhilohi_d = new double*[dim];
   double **Aimlolohi_d = new double*[dim];
   double **Aimhihilo_d = new double*[dim];
   double **Aimlohilo_d = new double*[dim];
   double **Aimhilolo_d = new double*[dim];
   double **Aimlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      rhsrehihihi_d[i] = rhsrehihihi[i];
      rhsrelohihi_d[i] = rhsrelohihi[i];
      rhsrehilohi_d[i] = rhsrehilohi[i];
      rhsrelolohi_d[i] = rhsrelolohi[i];
      rhsrehihilo_d[i] = rhsrehihilo[i];
      rhsrelohilo_d[i] = rhsrelohilo[i];
      rhsrehilolo_d[i] = rhsrehilolo[i];
      rhsrelololo_d[i] = rhsrelololo[i];
      rhsimhihihi_d[i] = rhsimhihihi[i];
      rhsimlohihi_d[i] = rhsimlohihi[i];
      rhsimhilohi_d[i] = rhsimhilohi[i];
      rhsimlolohi_d[i] = rhsimlolohi[i];
      rhsimhihilo_d[i] = rhsimhihilo[i];
      rhsimlohilo_d[i] = rhsimlohilo[i];
      rhsimhilolo_d[i] = rhsimhilolo[i];
      rhsimlololo_d[i] = rhsimlololo[i];
      Arehihihi_d[i] = new double[dim];
      Arelohihi_d[i] = new double[dim];
      Arehilohi_d[i] = new double[dim];
      Arelolohi_d[i] = new double[dim];
      Arehihilo_d[i] = new double[dim];
      Arelohilo_d[i] = new double[dim];
      Arehilolo_d[i] = new double[dim];
      Arelololo_d[i] = new double[dim];
      Aimhihihi_d[i] = new double[dim];
      Aimlohihi_d[i] = new double[dim];
      Aimhilohi_d[i] = new double[dim];
      Aimlolohi_d[i] = new double[dim];
      Aimhihilo_d[i] = new double[dim];
      Aimlohilo_d[i] = new double[dim];
      Aimhilolo_d[i] = new double[dim];
      Aimlololo_d[i] = new double[dim];

      for(int j=0; j<dim; j++)
      {
         Arehihihi_d[i][j] = Arehihihi[i][j];
         Arelohihi_d[i][j] = Arelohihi[i][j];
         Arehilohi_d[i][j] = Arehilohi[i][j];
         Arelolohi_d[i][j] = Arelolohi[i][j];
         Arehihilo_d[i][j] = Arehihilo[i][j];
         Arelohilo_d[i][j] = Arelohilo[i][j];
         Arehilolo_d[i][j] = Arehilolo[i][j];
         Arelololo_d[i][j] = Arelololo[i][j];
         Aimhihihi_d[i][j] = Aimhihihi[i][j];
         Aimlohihi_d[i][j] = Aimlohihi[i][j];
         Aimhilohi_d[i][j] = Aimhilohi[i][j];
         Aimlolohi_d[i][j] = Aimlolohi[i][j];
         Aimhihilo_d[i][j] = Aimhihilo[i][j];
         Aimlohilo_d[i][j] = Aimlohilo[i][j];
         Aimhilolo_d[i][j] = Aimhilolo[i][j];
         Aimlololo_d[i][j] = Aimlololo[i][j];
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
   double timelapsed_d,elapsedms;
   double invlapsed,mullapsed,sublapsed;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;

   cout << "-> GPU solves an upper triangular system ..." << endl;

   GPU_cmplx8_upper_tiled_solver
      (dim,sizetile,numtiles,
         Arehihihi_d,  Arelohihi_d,  Arehilohi_d,  Arelolohi_d,
         Arehihilo_d,  Arelohilo_d,  Arehilolo_d,  Arelololo_d,
         Aimhihihi_d,  Aimlohihi_d,  Aimhilohi_d,  Aimlolohi_d,
         Aimhihilo_d,  Aimlohilo_d,  Aimhilolo_d,  Aimlololo_d,
       rhsrehihihi_d,rhsrelohihi_d,rhsrehilohi_d,rhsrelolohi_d,
       rhsrehihilo_d,rhsrelohilo_d,rhsrehilolo_d,rhsrelololo_d,
       rhsimhihihi_d,rhsimlohihi_d,rhsimhilohi_d,rhsimlolohi_d,
       rhsimhihilo_d,rhsimlohilo_d,rhsimhilolo_d,rhsimlololo_d,
         xrehihihi_d,  xrelohihi_d,  xrehilohi_d,  xrelolohi_d,
         xrehihilo_d,  xrelohilo_d,  xrehilolo_d,  xrelololo_d,
         ximhihihi_d,  ximlohihi_d,  ximhilohi_d,  ximlolohi_d,
         ximhihilo_d,  ximlohilo_d,  ximhilolo_d,  ximlololo_d,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&timelapsed_d,
       &addcnt,&mulcnt,&divcnt);

   if(verbose > 0)
   {
      cout << "The matrix returned by the device :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihihi_d[i][j] << "  " << Arelohihi_d[i][j] << endl
                 << "            "
                 << Arehilohi_d[i][j] << "  " << Arelolohi_d[i][j] << endl
                 << "            "
                 << Arehihilo_d[i][j] << "  " << Arelohilo_d[i][j] << endl
                 << "            "
                 << Arehilolo_d[i][j] << "  " << Arelolohi_d[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihihi_d[i][j] << "  " << Aimlohihi_d[i][j] << endl
                 << "            "
                 << Aimhilohi_d[i][j] << "  " << Aimlolohi_d[i][j] << endl
                 << "            "
                 << Aimhihilo_d[i][j] << "  " << Aimlohilo_d[i][j] << endl
                 << "            "
                 << Aimhilolo_d[i][j] << "  " << Aimlololo_d[i][j] << endl;
         }
   }
   cout << scientific << setprecision(2);
   cout << "   Sum of errors on diagonal tiles : "
        << cmplx8_Diagonal_Difference_Sum
             (numtiles,sizetile,
              Arehihihi,  Arelohihi,  Arehilohi,  Arelolohi,
              Arehihilo,  Arelohilo,  Arehilolo,  Arelololo,
              Aimhihihi,  Aimlohihi,  Aimhilohi,  Aimlolohi,
              Aimhihilo,  Aimlohilo,  Aimhilolo,  Aimlololo,
              Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
              Arehihilo_d,Arelohilo_d,Arehilolo_d,Arelololo_d,
              Aimhihihi_d,Aimlohihi_d,Aimhilohi_d,Aimlolohi_d,
              Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d)
        << endl;

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
      cout << "GPU solution computed with tiling :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "x[" << i << "]re : "
              << xrehihihi_d[i] << "  " << xrelohihi_d[i] << endl
              << "         "
              << xrehilohi_d[i] << "  " << xrelolohi_d[i] << endl
              << "         "
              << xrehihilo_d[i] << "  " << xrelohilo_d[i] << endl
              << "         "
              << xrehilolo_d[i] << "  " << xrelololo_d[i] << endl;
         cout << "x[" << i << "]im : "
              << ximhihihi_d[i] << "  " << ximlohihi_d[i] << endl
              << "         "
              << ximhilohi_d[i] << "  " << ximlolohi_d[i] << endl
              << "         "
              << ximhihilo_d[i] << "  " << ximlohilo_d[i] << endl
              << "         "
              << ximhilolo_d[i] << "  " << ximlololo_d[i] << endl;
      }
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
   cout << "   Sum of GPU errors on solution : "
        << cmplx8_Difference_Sum(dim,
               solrehihihi,solrelohihi,solrehilohi,solrelolohi,
               solrehihilo,solrelohilo,solrehilolo,solrelololo,
               solimhihihi,solimlohihi,solimhilohi,solimlolohi,
               solimhihilo,solimlohilo,solimhilolo,solimlololo,
                 xrehihihi_d,xrelohihi_d,xrehilohi_d,xrelolohi_d,
                 xrehihilo_d,xrelohilo_d,xrehilolo_d,xrelololo_d,
                 ximhihihi_d,ximlohihi_d,ximhilohi_d,ximlolohi_d,
                 ximhihilo_d,ximlohilo_d,ximhilolo_d,ximlololo_d)
        << endl;
   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;

   write_dbl8_bstimeflops
     (sizetile,numtiles,1,invlapsed,mullapsed,sublapsed,elapsedms,
      timelapsed_d,addcnt,mulcnt,divcnt);

   for(int i=0; i<dim; i++)
   {
      free(Arehihihi[i]); free(Arelohihi[i]);
      free(Arehilohi[i]); free(Arelolohi[i]);
      free(Arehihilo[i]); free(Arelohilo[i]);
      free(Arehilolo[i]); free(Arelololo[i]);
      free(Aimhihihi[i]); free(Aimlohihi[i]);
      free(Aimhilohi[i]); free(Aimlolohi[i]);
      free(Aimhihilo[i]); free(Aimlohilo[i]);
      free(Aimhilolo[i]); free(Aimlololo[i]);
      free(Arehihihi_d[i]); free(Arelohihi_d[i]);
      free(Arehilohi_d[i]); free(Arelolohi_d[i]);
      free(Arehihilo_d[i]); free(Arelohilo_d[i]);
      free(Arehilolo_d[i]); free(Arelololo_d[i]);
      free(Aimhihihi_d[i]); free(Aimlohihi_d[i]);
      free(Aimhilohi_d[i]); free(Aimlolohi_d[i]);
      free(Aimhihilo_d[i]); free(Aimlohilo_d[i]);
      free(Aimhilolo_d[i]); free(Aimlololo_d[i]);
   }
   free(Arehihihi); free(Arelohihi); free(Arehilohi); free(Arelolohi);
   free(Arehihilo); free(Arelohilo); free(Arehilolo); free(Arelololo);
   free(Aimhihihi); free(Aimlohihi); free(Aimhilohi); free(Aimlolohi);
   free(Aimhihilo); free(Aimlohilo); free(Aimhilolo); free(Aimlololo);
   free(Arehihihi_d); free(Arelohihi_d); free(Arehilohi_d); free(Arelolohi_d);
   free(Arehihilo_d); free(Arelohilo_d); free(Arehilolo_d); free(Arelololo_d);
   free(Aimhihihi_d); free(Aimlohihi_d); free(Aimhilohi_d); free(Aimlolohi_d);
   free(Aimhihilo_d); free(Aimlohilo_d); free(Aimhilolo_d); free(Aimlololo_d);
   free(solrehihihi); free(solrelohihi); free(solrehilohi); free(solrelolohi);
   free(solrehihilo); free(solrelohilo); free(solrehilolo); free(solrelololo);
   free(solimhihihi); free(solimlohihi); free(solimhilohi); free(solimlolohi);
   free(solimhihilo); free(solimlohilo); free(solimhilolo); free(solimlololo);
   free(rhsrehihihi); free(rhsrelohihi); free(rhsrehilohi); free(rhsrelolohi);
   free(rhsrehihilo); free(rhsrelohilo); free(rhsrehilolo); free(rhsrelololo);
   free(rhsimhihihi); free(rhsimlohihi); free(rhsimhilohi); free(rhsimlolohi);
   free(rhsimhihilo); free(rhsimlohilo); free(rhsimhilolo); free(rhsimlololo);
   free(rhsrehihihi_d); free(rhsrelohihi_d);
   free(rhsrehilohi_d); free(rhsrelolohi_d);
   free(rhsrehihilo_d); free(rhsrelohilo_d);
   free(rhsrehilolo_d); free(rhsrelololo_d);
   free(rhsimhihihi_d); free(rhsimlohihi_d);
   free(rhsimhilohi_d); free(rhsimlolohi_d);
   free(rhsimhihilo_d); free(rhsimlohilo_d);
   free(rhsimhilolo_d); free(rhsimlololo_d);
   free(xrehihihi); free(xrelohihi); free(xrehilohi); free(xrelolohi);
   free(xrehihilo); free(xrelohilo); free(xrehilolo); free(xrelololo);
   free(ximhihihi); free(ximlohihi); free(ximhilohi); free(ximlolohi);
   free(ximhihilo); free(ximlohilo); free(ximhilolo); free(ximlololo);
   free(xrehihihi_d); free(xrelohihi_d); free(xrehilohi_d); free(xrelolohi_d);
   free(xrehihilo_d); free(xrelohilo_d); free(xrehilolo_d); free(xrelololo_d);
   free(ximhihihi_d); free(ximlohihi_d); free(ximhilohi_d); free(ximlolohi_d);
   free(ximhihilo_d); free(ximlohilo_d); free(ximhilolo_d); free(ximlololo_d);
}
