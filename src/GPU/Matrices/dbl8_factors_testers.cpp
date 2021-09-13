/* The file dbl8_factors_testers.cpp define the functions specified in
   the file dbl8_factors_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "octo_double_functions.h"
#include "random8_matrices.h"
#include "dbl8_factorizations.h"

using namespace std;

void test_factors_real8_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

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
   random_dbl8_matrix
      (dim,dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
               Ahihilo,Alohilo,Ahilolo,Alololo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
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

      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
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
   double *xhihihi = new double[dim];
   double *xlohihi = new double[dim];
   double *xhilohi = new double[dim];
   double *xlolohi = new double[dim];
   double *xhihilo = new double[dim];
   double *xlohilo = new double[dim];
   double *xhilolo = new double[dim];
   double *xlololo = new double[dim];
   int *pivots = new int[dim];

   CPU_dbl8_factors_lusolve
      (dim,Ahihihi,  Alohihi,  Ahilohi,  Alolohi,
           Ahihilo,  Alohilo,  Ahilolo,  Alololo,pivots,
         rhshihihi,rhslohihi,rhshilohi,rhslolohi,
         rhshihilo,rhslohilo,rhshilolo,rhslololo,
           xhihihi,  xlohihi,  xhilohi,  xlolohi,
           xhihilo,  xlohilo,  xhilolo,  xlololo);

   if(verbose > 0)
   {
      cout << "The computed solution :" << endl;
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
   double error = 0.0;
   for(int i=0; i<dim; i++)
      error = error + abs(xhihihi[i] - 1.0) + abs(xlohihi[i])
                    + abs(xhilohi[i]) + abs(xlolohi[i])
                    + abs(xhihilo[i]) + abs(xlohilo[i])
                    + abs(xhilolo[i]) + abs(xlololo[i]);

   cout << scientific << setprecision(2);
   cout << "Sum of errors : " << error << endl;

   for(int i=0; i<dim; i++)
   {
      free(Ahihihi[i]); free(Alohihi[i]);
      free(Ahilohi[i]); free(Alolohi[i]);
      free(Ahihilo[i]); free(Alohilo[i]);
      free(Ahilolo[i]); free(Alololo[i]);
   }
   free(Ahihihi); free(Alohihi); free(Ahilohi); free(Alolohi);
   free(Ahihilo); free(Alohilo); free(Ahilolo); free(Alololo);
   free(solhihihi); free(sollohihi); free(solhilohi); free(sollolohi);
   free(solhihilo); free(sollohilo); free(solhilolo); free(sollololo);
   free(rhshihihi); free(rhslohihi); free(rhshilohi); free(rhslolohi);
   free(rhshihilo); free(rhslohilo); free(rhshilolo); free(rhslololo);
   free(xhihihi); free(xlohihi); free(xhilohi); free(xlolohi);
   free(xhihilo); free(xlohilo); free(xhilolo); free(xlololo);
}

void test_factors_cmplx8_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

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
   random_cmplx8_matrix
      (dim,dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
               Arehihilo,Arelohilo,Arehilolo,Arelololo,
               Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
               Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
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
   double *solrehilohi = new double[dim];
   double *solrelohihi = new double[dim];
   double *solrelolohi = new double[dim];
   double *solrehihilo = new double[dim];
   double *solrehilolo = new double[dim];
   double *solrelohilo = new double[dim];
   double *solrelololo = new double[dim];
   double *solimhihihi = new double[dim];
   double *solimhilohi = new double[dim];
   double *solimlohihi = new double[dim];
   double *solimlolohi = new double[dim];
   double *solimhihilo = new double[dim];
   double *solimhilolo = new double[dim];
   double *solimlohilo = new double[dim];
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
   double *rhsrehilohi = new double[dim];
   double *rhsrelohihi = new double[dim];
   double *rhsrelolohi = new double[dim];
   double *rhsrehihilo = new double[dim];
   double *rhsrehilolo = new double[dim];
   double *rhsrelohilo = new double[dim];
   double *rhsrelololo = new double[dim];
   double *rhsimhihihi = new double[dim];
   double *rhsimhilohi = new double[dim];
   double *rhsimlohihi = new double[dim];
   double *rhsimlolohi = new double[dim];
   double *rhsimhihilo = new double[dim];
   double *rhsimhilolo = new double[dim];
   double *rhsimlohilo = new double[dim];
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
               &acc2hihihi,    &acc2lohihi,    &acc2hilohi,    &acc2lolohi,
               &acc2hihilo,    &acc2lohilo,    &acc2hilolo,    &acc2lololo);
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
                   acc1hihihi,     acc1lohihi,   acc1hilohi,   acc1lolohi,
                   acc1hihilo,     acc1lohilo,   acc1hilolo,   acc1lololo);
         odf_dec(&rhsrehihihi[i],&rhsrelohihi[i],
                 &rhsrehilohi[i],&rhsrelolohi[i],
                 &rhsrehihilo[i],&rhsrelohilo[i],
                 &rhsrehilolo[i],&rhsrelololo[i],
                   acc2hihihi,     acc2lohihi,   acc2hilohi,   acc2lolohi,
                   acc2hihilo,     acc2lohilo,   acc2hilolo,   acc2lololo);
         odf_inc(&rhsimhihihi[i],&rhsimlohihi[i],
                 &rhsimhilohi[i],&rhsimlolohi[i],
                 &rhsimhihilo[i],&rhsimlohilo[i],
                 &rhsimhilolo[i],&rhsimlololo[i],
                   acc3hihihi,     acc3lohihi,   acc3hilohi,   acc3lolohi,
                   acc3hihilo,     acc3lohilo,   acc3hilolo,   acc3lololo);
         odf_inc(&rhsimhihihi[i],&rhsimlohihi[i],
                 &rhsimhilohi[i],&rhsimlolohi[i],
                 &rhsimhihilo[i],&rhsimlohilo[i],
                 &rhsimhilolo[i],&rhsimlololo[i],
                   acc4hihihi,     acc4lohihi,   acc4hilohi,   acc4lolohi,
                   acc4hihilo,     acc4lohilo,   acc4hilolo,   acc4lololo);
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
   int *pivots = new int[dim];

   CPU_cmplx8_factors_lusolve
      (dim,Arehihihi,  Arelohihi,  Arehilohi,  Arelolohi,
           Arehihilo,  Arelohilo,  Arehilolo,  Arelololo,
           Aimhihihi,  Aimlohihi,  Aimhilohi,  Aimlolohi,
           Aimhihilo,  Aimlohilo,  Aimhilolo,  Aimlololo,pivots,
         rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
         rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
         rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
         rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
           xrehihihi,  xrelohihi,  xrehilohi,  xrelolohi,
           xrehihilo,  xrelohilo,  xrehilolo,  xrelololo,
           ximhihihi,  ximlohihi,  ximhilohi,  ximlolohi,
           ximhihilo,  ximlohilo,  ximhilolo,  ximlololo);

   if(verbose > 0)
   {
      cout << "The computed solution :" << endl;
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
   double error = 0.0;
   for(int i=0; i<dim; i++)
      error = error + abs(xrehihihi[i] - 1.0) + abs(xrelohihi[i])
                    + abs(xrehilohi[i]) + abs(xrelolohi[i])
                    + abs(xrehihilo[i]) + abs(xrelohilo[i])
                    + abs(xrehilolo[i]) + abs(xrelololo[i])
                    + abs(ximhihihi[i]) + abs(ximlohihi[i])
                    + abs(ximhilohi[i]) + abs(ximlolohi[i])
                    + abs(ximhihilo[i]) + abs(ximlohilo[i])
                    + abs(ximhilolo[i]) + abs(ximlololo[i]);

   cout << scientific << setprecision(2);
   cout << "Sum of errors : " << error << endl;

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
   }
   free(Arehihihi); free(Arelohihi); free(Arehilohi); free(Arelolohi);
   free(Arehihilo); free(Arelohilo); free(Arehilolo); free(Arelololo);
   free(Aimhihihi); free(Aimlohihi); free(Aimhilohi); free(Aimlolohi);
   free(Aimhihilo); free(Aimlohilo); free(Aimhilolo); free(Aimlololo);
   free(solrehihihi); free(solrelohihi); free(solrehilohi); free(solrelolohi);
   free(solrehihilo); free(solrelohilo); free(solrehilolo); free(solrelololo);
   free(solimhihihi); free(solimlohihi); free(solimhilohi); free(solimlolohi);
   free(solimhihilo); free(solimlohilo); free(solimhilolo); free(solimlololo);
   free(rhsrehihihi); free(rhsrelohihi); free(rhsrehilohi); free(rhsrelolohi);
   free(rhsrehihilo); free(rhsrelohilo); free(rhsrehilolo); free(rhsrelololo);
   free(rhsimhihihi); free(rhsimlohihi); free(rhsimhilohi); free(rhsimlolohi);
   free(rhsimhihilo); free(rhsimlohilo); free(rhsimhilolo); free(rhsimlololo);
   free(xrehihihi); free(xrelohihi); free(xrehilohi); free(xrelolohi);
   free(xrehihilo); free(xrelohilo); free(xrehilolo); free(xrelololo);
   free(ximhihihi); free(ximlohihi); free(ximhilohi); free(ximlolohi);
   free(ximhihilo); free(ximlohilo); free(ximhilolo); free(ximlololo);
}

int test_real8_qr_factors
 ( int nrows, int ncols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double tol, int verbose )
{
   double **QThihihi = new double*[nrows];
   double **QTlohihi = new double*[nrows];
   double **QThilohi = new double*[nrows];
   double **QTlolohi = new double*[nrows];
   double **QThihilo = new double*[nrows];
   double **QTlohilo = new double*[nrows];
   double **QThilolo = new double*[nrows];
   double **QTlololo = new double*[nrows];
   double **QTQhihihi = new double*[nrows];
   double **QTQlohihi = new double*[nrows];
   double **QTQhilohi = new double*[nrows];
   double **QTQlolohi = new double*[nrows];
   double **QTQhihilo = new double*[nrows];
   double **QTQlohilo = new double*[nrows];
   double **QTQhilolo = new double*[nrows];
   double **QTQlololo = new double*[nrows];
   double **QTAhihihi = new double*[nrows];
   double **QTAlohihi = new double*[nrows];
   double **QTAhilohi = new double*[nrows];
   double **QTAlolohi = new double*[nrows];
   double **QTAhihilo = new double*[nrows];
   double **QTAlohilo = new double*[nrows];
   double **QTAhilolo = new double*[nrows];
   double **QTAlololo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      QThihihi[i] = new double[nrows];
      QTlohihi[i] = new double[nrows];
      QThilohi[i] = new double[nrows];
      QTlolohi[i] = new double[nrows];
      QThihilo[i] = new double[nrows];
      QTlohilo[i] = new double[nrows];
      QThilolo[i] = new double[nrows];
      QTlololo[i] = new double[nrows];
      QTQhihihi[i] = new double[nrows];
      QTQlohihi[i] = new double[nrows];
      QTQhilohi[i] = new double[nrows];
      QTQlolohi[i] = new double[nrows];
      QTQhihilo[i] = new double[nrows];
      QTQlohilo[i] = new double[nrows];
      QTQhilolo[i] = new double[nrows];
      QTQlololo[i] = new double[nrows];
      QTAhihihi[i] = new double[ncols];
      QTAlohihi[i] = new double[ncols];
      QTAhilohi[i] = new double[ncols];
      QTAlolohi[i] = new double[ncols];
      QTAhihilo[i] = new double[ncols];
      QTAlohilo[i] = new double[ncols];
      QTAhilolo[i] = new double[ncols];
      QTAlololo[i] = new double[ncols];
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix Q :" << endl;
   }
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q[" << i << "][" << j << "] : "
                 << Qhihihi[i][j] << "  " << Qlohihi[i][j] << endl
                 << "          "
                 << Qhilohi[i][j] << "  " << Qlolohi[i][j] << endl
                 << "          "
                 << Qhihilo[i][j] << "  " << Qlohilo[i][j] << endl
                 << "          "
                 << Qhilolo[i][j] << "  " << Qlololo[i][j] << endl;
         QThihihi[j][i] = Qhihihi[i][j];
         QTlohihi[j][i] = Qlohihi[i][j];
         QThilohi[j][i] = Qhilohi[i][j];
         QTlolohi[j][i] = Qlolohi[i][j];
         QThihilo[j][i] = Qhihilo[i][j];
         QTlohilo[j][i] = Qlohilo[i][j];
         QThilolo[j][i] = Qhilolo[i][j];
         QTlololo[j][i] = Qlololo[i][j];
      }

   CPU_dbl8_factors_matmatmul
      (nrows,nrows,nrows,QThihihi, QTlohihi, QThilohi, QTlolohi,
                         QThihilo, QTlohilo, QThilolo, QTlololo,
                          Qhihihi,  Qlohihi,  Qhilohi,  Qlolohi,
                          Qhihilo,  Qlohilo,  Qhilolo,  Qlololo,
                        QTQhihihi,QTQlohihi,QTQhilohi,QTQlolohi,
                        QTQhihilo,QTQlohilo,QTQhilolo,QTQlololo);

   double errorQ = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q^T*Q[" << i << "][" << j << "] : "
                 << QTQhihihi[i][j] << "  " << QTQlohihi[i][j] << endl
                 << "              "
                 << QTQhilohi[i][j] << "  " << QTQlolohi[i][j] << endl
                 << "              "
                 << QTQhihilo[i][j] << "  " << QTQlohilo[i][j] << endl
                 << "              "
                 << QTQhilolo[i][j] << "  " << QTQlololo[i][j] << endl;

         if(i == j)
            errorQ = errorQ + fabs(QTQhihihi[i][j] - 1.0)
                            + fabs(QTQlohihi[i][j])
                            + fabs(QTQhilohi[i][j]) + fabs(QTQlolohi[i][j])
                            + fabs(QTQhihilo[i][j]) + fabs(QTQlohilo[i][j])
                            + fabs(QTQhilolo[i][j]) + fabs(QTQlololo[i][j]);
         else
            errorQ = errorQ + fabs(QTQhihihi[i][j]) + fabs(QTQlohihi[i][j])
                            + fabs(QTQhilohi[i][j]) + fabs(QTQlolohi[i][j])
                            + fabs(QTQhihilo[i][j]) + fabs(QTQlohilo[i][j])
                            + fabs(QTQhilolo[i][j]) + fabs(QTQlololo[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix R :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhihihi[i][j] << "  " << Rlohihi[i][j] << endl
                 << "          "
                 << Rhilohi[i][j] << "  " << Rlolohi[i][j] << endl
                 << "          "
                 << Rhihilo[i][j] << "  " << Rlohilo[i][j] << endl
                 << "          "
                 << Rhilolo[i][j] << "  " << Rlololo[i][j] << endl;
   }
   CPU_dbl8_factors_matmatmul
      (nrows,nrows,ncols,QThihihi, QTlohihi, QThilohi, QTlolohi,
                         QThihilo, QTlohilo, QThilolo, QTlololo,
                          Ahihihi,  Alohihi,  Ahilohi,  Alolohi,
                          Ahihilo,  Alohilo,  Ahilolo,  Alololo,
                        QTAhihihi,QTAlohihi,QTAhilohi,QTAlolohi,
                        QTAhihilo,QTAlohilo,QTAhilolo,QTAlololo);

   double errorR = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(verbose > 0)
            cout << "Q^T*A[" << i << "][" << j << "] : "
                 << QTAhihihi[i][j] << "  " << QTAlohihi[i][j] << endl
                 << "              "
                 << QTAhilohi[i][j] << "  " << QTAlolohi[i][j] << endl
                 << "              "
                 << QTAhihilo[i][j] << "  " << QTAlohilo[i][j] << endl
                 << "              "
                 << QTAhilolo[i][j] << "  " << QTAlololo[i][j] << endl;

         errorR = errorR + fabs(Rhihihi[i][j] - QTAhihihi[i][j])
                         + fabs(Rhilohi[i][j] - QTAhilohi[i][j])
                         + fabs(Rlohihi[i][j] - QTAlohihi[i][j])
                         + fabs(Rlolohi[i][j] - QTAlolohi[i][j])
                         + fabs(Rhihilo[i][j] - QTAhihilo[i][j])
                         + fabs(Rhilolo[i][j] - QTAhilolo[i][j])
                         + fabs(Rlohilo[i][j] - QTAlohilo[i][j])
                         + fabs(Rlololo[i][j] - QTAlololo[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   for(int i=0; i<nrows; i++)
   {
      free(QThihihi[i]); free(QTQhihihi[i]); free(QTAhihihi[i]);
      free(QTlohihi[i]); free(QTQlohihi[i]); free(QTAlohihi[i]);
      free(QThilohi[i]); free(QTQhilohi[i]); free(QTAhilohi[i]);
      free(QTlolohi[i]); free(QTQlolohi[i]); free(QTAlolohi[i]);
      free(QThihilo[i]); free(QTQhihilo[i]); free(QTAhihilo[i]);
      free(QTlohilo[i]); free(QTQlohilo[i]); free(QTAlohilo[i]);
      free(QThilolo[i]); free(QTQhilolo[i]); free(QTAhilolo[i]);
      free(QTlololo[i]); free(QTQlololo[i]); free(QTAlololo[i]);
   }
   free(QThihihi); free(QTQhihihi); free(QTAhihihi);
   free(QTlohihi); free(QTQlohihi); free(QTAlohihi);
   free(QThilohi); free(QTQhilohi); free(QTAhilohi);
   free(QTlolohi); free(QTQlolohi); free(QTAlolohi);
   free(QThihilo); free(QTQhihilo); free(QTAhihilo);
   free(QTlohilo); free(QTQlohilo); free(QTAlohilo);
   free(QThilolo); free(QTQhilolo); free(QTAhilolo);
   free(QTlololo); free(QTQlololo); free(QTAlololo);

   return int(errorQ + errorR > tol);
}

int test_real8_qr_factors_probe
 ( int nrows, int ncols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double tol, int nbprobes, int verbose )
{
   int rowidx,colidx;
   double Qsumhihihi,Qsumlohihi,Qsumhilohi,Qsumlolohi;
   double Qsumhihilo,Qsumlohilo,Qsumhilolo,Qsumlololo;
   double Rsumhihihi,Rsumlohihi,Rsumhilohi,Rsumlolohi;
   double Rsumhihilo,Rsumlohilo,Rsumhilolo,Rsumlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   double errorQ = 0.0;
   double errorR = 0.0;

   for(int p=0; p<nbprobes; p++)
   {
      rowidx = rand() % nrows;
      colidx = rand() % ncols;

      if(verbose > 0)
      {
         cout << "Probing row index : " << rowidx
              << ", column index : " << colidx << "." << endl;
      }
      Qsumhihihi = 0.0; Qsumlohihi = 0.0;
      Qsumhilohi = 0.0; Qsumlolohi = 0.0;
      Qsumhihilo = 0.0; Qsumlohilo = 0.0;
      Qsumhilolo = 0.0; Qsumlololo = 0.0;
      Rsumhihihi = 0.0; Rsumlohihi = 0.0;
      Rsumhilohi = 0.0; Rsumlolohi = 0.0;
      Rsumhihilo = 0.0; Rsumlohilo = 0.0;
      Rsumhilolo = 0.0; Rsumlololo = 0.0;

      for(int i=0; i<nrows; i++)
      {
         odf_mul(Qhihihi[i][rowidx],Qlohihi[i][rowidx],
                 Qhilohi[i][rowidx],Qlolohi[i][rowidx],
                 Qhihilo[i][rowidx],Qlohilo[i][rowidx],
                 Qhilolo[i][rowidx],Qlololo[i][rowidx],
                 Qhihihi[i][colidx],Qlohihi[i][colidx],
                 Qhilohi[i][colidx],Qlolohi[i][colidx],
                 Qhihilo[i][colidx],Qlohilo[i][colidx],
                 Qhilolo[i][colidx],Qlololo[i][colidx],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&Qsumhihihi,&Qsumlohihi,&Qsumhilohi,&Qsumlolohi,
                 &Qsumhihilo,&Qsumlohilo,&Qsumhilolo,&Qsumlololo,
                   acchihihi,  acclohihi,  acchilohi,  acclolohi,
                   acchihilo,  acclohilo,  acchilolo,  acclololo);
         odf_mul(Qhihihi[i][rowidx],Qlohihi[i][rowidx],
                 Qhilohi[i][rowidx],Qlolohi[i][rowidx],
                 Qhihilo[i][rowidx],Qlohilo[i][rowidx],
                 Qhilolo[i][rowidx],Qlololo[i][rowidx],
                 Ahihihi[i][colidx],Alohihi[i][colidx],
                 Ahilohi[i][colidx],Alolohi[i][colidx],
                 Ahihilo[i][colidx],Alohilo[i][colidx],
                 Ahilolo[i][colidx],Alololo[i][colidx],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&Rsumhihihi,&Rsumlohihi,&Rsumhilohi,&Rsumlolohi,
                 &Rsumhihilo,&Rsumlohilo,&Rsumhilolo,&Rsumlololo,
                   acchihihi,  acclohihi,  acchilohi,  acclolohi,
                   acchihilo,  acclohilo,  acchilolo,  acclololo);
      }
      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "] : "
              << Qsumhihihi << "  " << Qsumlohihi << endl
              << "              "
              << Qsumhilohi << "  " << Qsumlolohi << endl
              << "              "
              << Qsumhihilo << "  " << Qsumlohilo << endl
              << "              "
              << Qsumhilolo << "  " << Qsumlololo << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "] : "
              << Rsumhihihi << "  " << Rsumlohihi << endl
              << "              "
              << Rsumhilohi << "  " << Rsumlolohi << endl
              << "              "
              << Rsumhihilo << "  " << Rsumlohilo << endl
              << "              "
              << Rsumhilolo << "  " << Rsumlololo << endl;
         cout << "    R[" << rowidx << "][" << colidx << "] : "
              << Rhihihi[rowidx][colidx] << "  "
              << Rlohihi[rowidx][colidx] << endl
              << "              "
              << Rhilohi[rowidx][colidx] << "  "
              << Rlolohi[rowidx][colidx] << endl
              << "              "
              << Rhihilo[rowidx][colidx] << "  "
              << Rlohilo[rowidx][colidx] << endl
              << "              "
              << Rhilolo[rowidx][colidx] << "  "
              << Rlololo[rowidx][colidx] << endl;
      }
      if(rowidx == colidx)
         errorQ = errorQ + fabs(Qsumhihihi - 1.0) + fabs(Qsumlohihi)
                         + fabs(Qsumhilohi) + fabs(Qsumlolohi)
                         + fabs(Qsumhihilo) + fabs(Qsumlohilo)
                         + fabs(Qsumhilolo) + fabs(Qsumlololo);
      else
         errorQ = errorQ + fabs(Qsumhihihi) + fabs(Qsumlohihi)
                         + fabs(Qsumhilohi) + fabs(Qsumlolohi)
                         + fabs(Qsumhihilo) + fabs(Qsumlohilo)
                         + fabs(Qsumhilolo) + fabs(Qsumlololo);

      errorR = errorR + fabs(Rsumhihihi - Rhihihi[rowidx][colidx])
                      + fabs(Rsumlohihi - Rlohihi[rowidx][colidx])
                      + fabs(Rsumhilohi - Rhilohi[rowidx][colidx])
                      + fabs(Rsumlolohi - Rlolohi[rowidx][colidx])
                      + fabs(Rsumhihilo - Rhihilo[rowidx][colidx])
                      + fabs(Rsumlohilo - Rlohilo[rowidx][colidx])
                      + fabs(Rsumhilolo - Rhilolo[rowidx][colidx])
                      + fabs(Rsumlololo - Rlololo[rowidx][colidx]);
   }
   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   return int(errorQ + errorR > tol);
}

int test_cmplx8_qr_factors
 ( int nrows, int ncols,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double tol, int verbose )
{
   double **QHrehihihi = new double*[nrows];
   double **QHrelohihi = new double*[nrows];
   double **QHrehilohi = new double*[nrows];
   double **QHrelolohi = new double*[nrows];
   double **QHrehihilo = new double*[nrows];
   double **QHrelohilo = new double*[nrows];
   double **QHrehilolo = new double*[nrows];
   double **QHrelololo = new double*[nrows];
   double **QHimhihihi = new double*[nrows];
   double **QHimlohihi = new double*[nrows];
   double **QHimhilohi = new double*[nrows];
   double **QHimlolohi = new double*[nrows];
   double **QHimhihilo = new double*[nrows];
   double **QHimlohilo = new double*[nrows];
   double **QHimhilolo = new double*[nrows];
   double **QHimlololo = new double*[nrows];
   double **QHQrehihihi = new double*[nrows];
   double **QHQrelohihi = new double*[nrows];
   double **QHQrehilohi = new double*[nrows];
   double **QHQrelolohi = new double*[nrows];
   double **QHQrehihilo = new double*[nrows];
   double **QHQrelohilo = new double*[nrows];
   double **QHQrehilolo = new double*[nrows];
   double **QHQrelololo = new double*[nrows];
   double **QHQimhihihi = new double*[nrows];
   double **QHQimlohihi = new double*[nrows];
   double **QHQimhilohi = new double*[nrows];
   double **QHQimlolohi = new double*[nrows];
   double **QHQimhihilo = new double*[nrows];
   double **QHQimlohilo = new double*[nrows];
   double **QHQimhilolo = new double*[nrows];
   double **QHQimlololo = new double*[nrows];
   double **QHArehihihi = new double*[nrows];
   double **QHArelohihi = new double*[nrows];
   double **QHArehilohi = new double*[nrows];
   double **QHArelolohi = new double*[nrows];
   double **QHArehihilo = new double*[nrows];
   double **QHArelohilo = new double*[nrows];
   double **QHArehilolo = new double*[nrows];
   double **QHArelololo = new double*[nrows];
   double **QHAimhihihi = new double*[nrows];
   double **QHAimlohihi = new double*[nrows];
   double **QHAimhilohi = new double*[nrows];
   double **QHAimlolohi = new double*[nrows];
   double **QHAimhihilo = new double*[nrows];
   double **QHAimlohilo = new double*[nrows];
   double **QHAimhilolo = new double*[nrows];
   double **QHAimlololo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      QHrehihihi[i] = new double[nrows];
      QHrelohihi[i] = new double[nrows];
      QHrehilohi[i] = new double[nrows];
      QHrelolohi[i] = new double[nrows];
      QHrehihilo[i] = new double[nrows];
      QHrelohilo[i] = new double[nrows];
      QHrehilolo[i] = new double[nrows];
      QHrelololo[i] = new double[nrows];
      QHimhihihi[i] = new double[nrows];
      QHimlohihi[i] = new double[nrows];
      QHimhilohi[i] = new double[nrows];
      QHimlolohi[i] = new double[nrows];
      QHimhihilo[i] = new double[nrows];
      QHimlohilo[i] = new double[nrows];
      QHimhilolo[i] = new double[nrows];
      QHimlololo[i] = new double[nrows];
      QHQrehihihi[i] = new double[nrows];
      QHQrelohihi[i] = new double[nrows];
      QHQrehilohi[i] = new double[nrows];
      QHQrelolohi[i] = new double[nrows];
      QHQrehihilo[i] = new double[nrows];
      QHQrelohilo[i] = new double[nrows];
      QHQrehilolo[i] = new double[nrows];
      QHQrelololo[i] = new double[nrows];
      QHQimhihihi[i] = new double[nrows];
      QHQimlohihi[i] = new double[nrows];
      QHQimhilohi[i] = new double[nrows];
      QHQimlolohi[i] = new double[nrows];
      QHQimhihilo[i] = new double[nrows];
      QHQimlohilo[i] = new double[nrows];
      QHQimhilolo[i] = new double[nrows];
      QHQimlololo[i] = new double[nrows];
      QHArehihihi[i] = new double[ncols];
      QHArelohihi[i] = new double[ncols];
      QHArehilohi[i] = new double[ncols];
      QHArelolohi[i] = new double[ncols];
      QHArehihilo[i] = new double[ncols];
      QHArelohilo[i] = new double[ncols];
      QHArehilolo[i] = new double[ncols];
      QHArelololo[i] = new double[ncols];
      QHAimhihihi[i] = new double[ncols];
      QHAimlohihi[i] = new double[ncols];
      QHAimhilohi[i] = new double[ncols];
      QHAimlolohi[i] = new double[ncols];
      QHAimhihilo[i] = new double[ncols];
      QHAimlohilo[i] = new double[ncols];
      QHAimhilolo[i] = new double[ncols];
      QHAimlololo[i] = new double[ncols];
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix Q :" << endl;
   }
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
         {
            cout << "Q[" << i << "][" << j << "]re : "
                 << Qrehihihi[i][j] << "  " << Qrelohihi[i][j] << endl
                 << "            "
                 << Qrehilohi[i][j] << "  " << Qrelolohi[i][j] << endl
                 << "            "
                 << Qrehihilo[i][j] << "  " << Qrelohilo[i][j] << endl
                 << "            "
                 << Qrehilolo[i][j] << "  " << Qrelololo[i][j] << endl;
            cout << "Q[" << i << "][" << j << "]im : "
                 << Qimhihihi[i][j] << "  " << Qimlohihi[i][j] << endl
                 << "            "
                 << Qimhilohi[i][j] << "  " << Qimlolohi[i][j] << endl
                 << "            "
                 << Qimhihilo[i][j] << "  " << Qimlohilo[i][j] << endl
                 << "            "
                 << Qimhilolo[i][j] << "  " << Qimlololo[i][j] << endl;
         }
         QHrehihihi[j][i] = Qrehihihi[i][j];
         QHrelohihi[j][i] = Qrelohihi[i][j];
         QHrehilohi[j][i] = Qrehilohi[i][j];
         QHrelolohi[j][i] = Qrelolohi[i][j];
         QHrehihilo[j][i] = Qrehihilo[i][j];
         QHrelohilo[j][i] = Qrelohilo[i][j];
         QHrehilolo[j][i] = Qrehilolo[i][j];
         QHrelololo[j][i] = Qrelololo[i][j];
         QHimhihihi[j][i] = Qimhihihi[i][j];
         QHimlohihi[j][i] = Qimlohihi[i][j];
         QHimhilohi[j][i] = Qimhilohi[i][j];
         QHimlolohi[j][i] = Qimlolohi[i][j];
         QHimhihilo[j][i] = Qimhihilo[i][j];
         QHimlohilo[j][i] = Qimlohilo[i][j];
         QHimhilolo[j][i] = Qimhilolo[i][j];
         QHimlololo[j][i] = Qimlololo[i][j];

         odf_minus(&QHimhihihi[j][i],&QHimlohihi[j][i],
                   &QHimhilohi[j][i],&QHimlolohi[j][i],
                   &QHimhihilo[j][i],&QHimlohilo[j][i],
                   &QHimhilolo[j][i],&QHimlololo[j][i]); // Hermitian transpose
      }

   CPU_cmplx8_factors_matmatmul
      (nrows,nrows,nrows,QHrehihihi, QHrelohihi,  QHrehilohi, QHrelolohi,
                         QHrehihilo, QHrelohilo,  QHrehilolo, QHrelololo,
                         QHimhihihi, QHimlohihi,  QHimhilohi, QHimlolohi,
                         QHimhihilo, QHimlohilo,  QHimhilolo, QHimlololo,
                          Qrehihihi,  Qrelohihi,  Qrehilohi,  Qrelolohi, 
                          Qrehihilo,  Qrelohilo,  Qrehilolo,  Qrelololo, 
                          Qimhihihi,  Qimlohihi,  Qimhilohi,  Qimlolohi,
                          Qimhihilo,  Qimlohilo,  Qimhilolo,  Qimlololo,
                        QHQrehihihi,QHQrelohihi,QHQrehilohi,QHQrelolohi,
                        QHQrehihilo,QHQrelohilo,QHQrehilolo,QHQrelololo,
                        QHQimhihihi,QHQimlohihi,QHQimhilohi,QHQimlolohi,
                        QHQimhihilo,QHQimlohilo,QHQimhilolo,QHQimlololo);

   double errorQ = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
         {
            cout << "Q^H*Q[" << i << "][" << j << "]re : "
                 << QHQrehihihi[i][j] << "  " << QHQrelohihi[i][j] << endl
                 << "                "
                 << QHQrehilohi[i][j] << "  " << QHQrelolohi[i][j] << endl
                 << "                "
                 << QHQrehihilo[i][j] << "  " << QHQrelohilo[i][j] << endl
                 << "                "
                 << QHQrehilolo[i][j] << "  " << QHQrelololo[i][j] << endl;
            cout << "Q^H*Q[" << i << "][" << j << "]im : "
                 << QHQimhihihi[i][j] << "  " << QHQimlohihi[i][j] << endl
                 << "                "
                 << QHQimhilohi[i][j] << "  " << QHQimlolohi[i][j] << endl
                 << "                "
                 << QHQimhihilo[i][j] << "  " << QHQimlohilo[i][j] << endl
                 << "                "
                 << QHQimhilolo[i][j] << "  " << QHQimlololo[i][j] << endl;
         }
         if(i == j)
            errorQ = errorQ + abs(QHQrehihihi[i][j] - 1.0)
                            + abs(QHQrelohihi[i][j])
                            + abs(QHQrehilohi[i][j]) + abs(QHQrelolohi[i][j])
                            + abs(QHQrehihilo[i][j]) + abs(QHQrelohilo[i][j])
                            + abs(QHQrehilolo[i][j]) + abs(QHQrelololo[i][j])
                            + abs(QHQimhihihi[i][j]) + abs(QHQimlohihi[i][j])
                            + abs(QHQimhilohi[i][j]) + abs(QHQimlolohi[i][j])
                            + abs(QHQimhihilo[i][j]) + abs(QHQimlohilo[i][j])
                            + abs(QHQimhilolo[i][j]) + abs(QHQimlololo[i][j]);
         else
            errorQ = errorQ + abs(QHQrehihihi[i][j]) + abs(QHQrelohihi[i][j])
                            + abs(QHQrehilohi[i][j]) + abs(QHQrelolohi[i][j])
                            + abs(QHQrehihilo[i][j]) + abs(QHQrelohilo[i][j])
                            + abs(QHQrehilolo[i][j]) + abs(QHQrelololo[i][j])
                            + abs(QHQimhihihi[i][j]) + abs(QHQimlohihi[i][j])
                            + abs(QHQimhilohi[i][j]) + abs(QHQimlolohi[i][j])
                            + abs(QHQimhihilo[i][j]) + abs(QHQimlohilo[i][j])
                            + abs(QHQimhilolo[i][j]) + abs(QHQimlololo[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^H*Q - I| : " << errorQ << endl;

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix R :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "R[" << i << "][" << j << "]re : "
                 << Rrehihihi[i][j] << "  " << Rrelohihi[i][j] << endl
                 << "            "
                 << Rrehilohi[i][j] << "  " << Rrelolohi[i][j] << endl
                 << "            "
                 << Rrehihilo[i][j] << "  " << Rrelohilo[i][j] << endl
                 << "            "
                 << Rrehilolo[i][j] << "  " << Rrelololo[i][j] << endl;
            cout << "R[" << i << "][" << j << "]im : "
                 << Rimhihihi[i][j] << "  " << Rimlohihi[i][j] << endl
                 << "            "
                 << Rimhilohi[i][j] << "  " << Rimlolohi[i][j] << endl
                 << "            "
                 << Rimhihilo[i][j] << "  " << Rimlohilo[i][j] << endl
                 << "            "
                 << Rimhilolo[i][j] << "  " << Rimlololo[i][j] << endl;
         }
   }
   CPU_cmplx8_factors_matmatmul
      (nrows,nrows,ncols, QHrehihihi, QHrelohihi, QHrehilohi, QHrelolohi,
                          QHrehihilo, QHrelohilo, QHrehilolo, QHrelololo,
                          QHimhihihi, QHimlohihi, QHimhilohi, QHimlolohi,
                          QHimhihilo, QHimlohilo, QHimhilolo, QHimlololo,
                           Arehihihi,  Arelohihi,  Arehilohi,  Arelolohi, 
                           Arehihilo,  Arelohilo,  Arehilolo,  Arelololo, 
                           Aimhihihi,  Aimlohihi,  Aimhilohi,  Aimlolohi,
                           Aimhihilo,  Aimlohilo,  Aimhilolo,  Aimlololo,
                         QHArehihihi,QHArelohihi,QHArehilohi,QHArelolohi,
                         QHArehihilo,QHArelohilo,QHArehilolo,QHArelololo,
                         QHAimhihihi,QHAimlohihi,QHAimhilohi,QHAimlolohi,
                         QHAimhihilo,QHAimlohilo,QHAimhilolo,QHAimlololo);

   double errorR = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(verbose > 0)
         {
            cout << "Q^H*A[" << i << "][" << j << "]re : "
                 << QHArehihihi[i][j] << "  " << QHArelohihi[i][j] << endl
                 << "                "
                 << QHArehilohi[i][j] << "  " << QHArelolohi[i][j] << endl
                 << "                "
                 << QHArehihilo[i][j] << "  " << QHArelohilo[i][j] << endl
                 << "                "
                 << QHArehilolo[i][j] << "  " << QHArelololo[i][j] << endl;
            cout << "Q^H*A[" << i << "][" << j << "]im : "
                 << QHAimhihihi[i][j] << "  " << QHAimlohihi[i][j] << endl
                 << "                "
                 << QHAimhilohi[i][j] << "  " << QHAimlolohi[i][j] << endl
                 << "                "
                 << QHAimhihilo[i][j] << "  " << QHAimlohilo[i][j] << endl
                 << "                "
                 << QHAimhilolo[i][j] << "  " << QHAimlololo[i][j] << endl;
         }
         errorR = errorR + abs(Rrehihihi[i][j] - QHArehihihi[i][j])
                         + abs(Rrelohihi[i][j] - QHArelohihi[i][j])
                         + abs(Rrehilohi[i][j] - QHArehilohi[i][j])
                         + abs(Rrelolohi[i][j] - QHArelolohi[i][j])
                         + abs(Rrehihilo[i][j] - QHArehihilo[i][j])
                         + abs(Rrelohilo[i][j] - QHArelohilo[i][j])
                         + abs(Rrehilolo[i][j] - QHArehilolo[i][j])
                         + abs(Rrelololo[i][j] - QHArelololo[i][j])
                         + abs(Rimhihihi[i][j] - QHAimhihihi[i][j])
                         + abs(Rimlohihi[i][j] - QHAimlohihi[i][j])
                         + abs(Rimhilohi[i][j] - QHAimhilohi[i][j])
                         + abs(Rimlolohi[i][j] - QHAimlolohi[i][j])
                         + abs(Rimhihilo[i][j] - QHAimhihilo[i][j])
                         + abs(Rimlohilo[i][j] - QHAimlohilo[i][j])
                         + abs(Rimhilolo[i][j] - QHAimhilolo[i][j])
                         + abs(Rimlololo[i][j] - QHAimlololo[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^H*A - R| : " << errorR << endl;

   for(int i=0; i<nrows; i++)
   {
      free(QHrehihihi[i]); free(QHimhihihi[i]);
      free(QHrelohihi[i]); free(QHimlohihi[i]);
      free(QHrehilohi[i]); free(QHimhilohi[i]);
      free(QHrelolohi[i]); free(QHimlolohi[i]);
      free(QHrehihilo[i]); free(QHimhihilo[i]);
      free(QHrelohilo[i]); free(QHimlohilo[i]);
      free(QHrehilolo[i]); free(QHimhilolo[i]);
      free(QHrelololo[i]); free(QHimlololo[i]);
      free(QHQrehihihi[i]); free(QHQimhihihi[i]);
      free(QHQrelohihi[i]); free(QHQimlohihi[i]);
      free(QHQrehilohi[i]); free(QHQimhilohi[i]);
      free(QHQrelolohi[i]); free(QHQimlolohi[i]);
      free(QHQrehihilo[i]); free(QHQimhihilo[i]);
      free(QHQrelohilo[i]); free(QHQimlohilo[i]);
      free(QHQrehilolo[i]); free(QHQimhilolo[i]);
      free(QHQrelololo[i]); free(QHQimlololo[i]);
      free(QHArehihihi[i]); free(QHAimhihihi[i]);
      free(QHArelohihi[i]); free(QHAimlohihi[i]);
      free(QHArehilohi[i]); free(QHAimhilohi[i]);
      free(QHArelolohi[i]); free(QHAimlolohi[i]);
      free(QHArehihilo[i]); free(QHAimhihilo[i]);
      free(QHArelohilo[i]); free(QHAimlohilo[i]);
      free(QHArehilolo[i]); free(QHAimhilolo[i]);
      free(QHArelololo[i]); free(QHAimlololo[i]);
   }
   free(QHrehihihi); free(QHQrehihihi); free(QHArehihihi);
   free(QHrelohihi); free(QHQrelohihi); free(QHArelohihi);
   free(QHrehilohi); free(QHQrehilohi); free(QHArehilohi);
   free(QHrelolohi); free(QHQrelolohi); free(QHArelolohi);
   free(QHrehihilo); free(QHQrehihilo); free(QHArehihilo);
   free(QHrelohilo); free(QHQrelohilo); free(QHArelohilo);
   free(QHrehilolo); free(QHQrehilolo); free(QHArehilolo);
   free(QHrelololo); free(QHQrelololo); free(QHArelololo);
   free(QHimhihihi); free(QHQimhihihi); free(QHAimhihihi);
   free(QHimlohihi); free(QHQimlohihi); free(QHAimlohihi);
   free(QHimhilohi); free(QHQimhilohi); free(QHAimhilohi);
   free(QHimlolohi); free(QHQimlolohi); free(QHAimlolohi);
   free(QHimhihilo); free(QHQimhihilo); free(QHAimhihilo);
   free(QHimlohilo); free(QHQimlohilo); free(QHAimlohilo);
   free(QHimhilolo); free(QHQimhilolo); free(QHAimhilolo);
   free(QHimlololo); free(QHQimlololo); free(QHAimlololo);

   return int(errorQ + errorR > tol);
}

int test_cmplx8_qr_factors_probe
 ( int nrows, int ncols,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double tol, int nbprobes, int verbose )
{
   int rowidx,colidx;
   double Qsumrehihihi,Qsumrelohihi,Qsumrehilohi,Qsumrelolohi;
   double Qsumrehihilo,Qsumrelohilo,Qsumrehilolo,Qsumrelololo;
   double Rsumrehihihi,Rsumrelohihi,Rsumrehilohi,Rsumrelolohi;
   double Rsumrehihilo,Rsumrelohilo,Rsumrehilolo,Rsumrelololo;
   double Qsumimhihihi,Qsumimlohihi,Qsumimhilohi,Qsumimlolohi;
   double Qsumimhihilo,Qsumimlohilo,Qsumimhilolo,Qsumimlololo;
   double Rsumimhihihi,Rsumimlohihi,Rsumimhilohi,Rsumimlolohi;
   double Rsumimhihilo,Rsumimlohilo,Rsumimhilolo,Rsumimlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   double errorQ = 0.0;
   double errorR = 0.0;

   for(int p=0; p<nbprobes; p++)
   {
      rowidx = rand() % nrows;
      colidx = rand() % ncols;

      if(verbose > 0)
      {
         cout << "Probing row index : " << rowidx
              << ", column index : " << colidx << "." << endl;
      }
      Qsumrehihihi = 0.0; Qsumrelohihi = 0.0;
      Qsumrehilohi = 0.0; Qsumrelolohi = 0.0;
      Qsumrehihilo = 0.0; Qsumrelohilo = 0.0;
      Qsumrehilolo = 0.0; Qsumrelololo = 0.0;
      Qsumimhihihi = 0.0; Qsumimlohihi = 0.0;
      Qsumimhilohi = 0.0; Qsumimlolohi = 0.0;
      Qsumimhihilo = 0.0; Qsumimlohilo = 0.0;
      Qsumimhilolo = 0.0; Qsumimlololo = 0.0;
      Rsumrehihihi = 0.0; Rsumrelohihi = 0.0;
      Rsumrehilohi = 0.0; Rsumrelolohi = 0.0;
      Rsumrehihilo = 0.0; Rsumrelohilo = 0.0;
      Rsumrehilolo = 0.0; Rsumrelololo = 0.0;
      Rsumimhihihi = 0.0; Rsumimlohihi = 0.0;
      Rsumimhilohi = 0.0; Rsumimlolohi = 0.0;
      Rsumimhihilo = 0.0; Rsumimlohilo = 0.0;
      Rsumimhilolo = 0.0; Rsumimlololo = 0.0;

      for(int i=0; i<nrows; i++)
      {
         // multiply Q^H with Q
         odf_mul(Qrehihihi[i][rowidx],Qrelohihi[i][rowidx],
                 Qrehilohi[i][rowidx],Qrelolohi[i][rowidx],
                 Qrehihilo[i][rowidx],Qrelohilo[i][rowidx],
                 Qrehilolo[i][rowidx],Qrelololo[i][rowidx],
                 Qrehihihi[i][colidx],Qrelohihi[i][colidx],
                 Qrehilohi[i][colidx],Qrelolohi[i][colidx],
                 Qrehihilo[i][colidx],Qrelohilo[i][colidx],
                 Qrehilolo[i][colidx],Qrelololo[i][colidx],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&Qsumrehihihi,&Qsumrelohihi,&Qsumrehilohi,&Qsumrelolohi,
                 &Qsumrehihilo,&Qsumrelohilo,&Qsumrehilolo,&Qsumrelololo,
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
         odf_mul(Qimhihihi[i][rowidx],Qimlohihi[i][rowidx],
                 Qimhilohi[i][rowidx],Qimlolohi[i][rowidx],
                 Qimhihilo[i][rowidx],Qimlohilo[i][rowidx],
                 Qimhilolo[i][rowidx],Qimlololo[i][rowidx],
                 Qimhihihi[i][colidx],Qimlohihi[i][colidx],
                 Qimhilohi[i][colidx],Qimlolohi[i][colidx],
                 Qimhihilo[i][colidx],Qimlohilo[i][colidx],
                 Qimhilolo[i][colidx],Qimlololo[i][colidx],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&Qsumrehihihi,&Qsumrelohihi,&Qsumrehilohi,&Qsumrelolohi,
                 &Qsumrehihilo,&Qsumrelohilo,&Qsumrehilolo,&Qsumrelololo,
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
         odf_mul(Qimhihihi[i][rowidx],Qimlohihi[i][rowidx],
                 Qimhilohi[i][rowidx],Qimlolohi[i][rowidx],
                 Qimhihilo[i][rowidx],Qimlohilo[i][rowidx],
                 Qimhilolo[i][rowidx],Qimlololo[i][rowidx],
                 Qrehihihi[i][colidx],Qrelohihi[i][colidx],
                 Qrehilohi[i][colidx],Qrelolohi[i][colidx],
                 Qrehihilo[i][colidx],Qrelohilo[i][colidx],
                 Qrehilolo[i][colidx],Qrelololo[i][colidx],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_dec(&Qsumimhihihi,&Qsumimlohihi,&Qsumimhilohi,&Qsumimlolohi,
                 &Qsumimhihilo,&Qsumimlohilo,&Qsumimhilolo,&Qsumimlololo,
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
         odf_mul(Qrehihihi[i][rowidx],Qrelohihi[i][rowidx],
                 Qrehilohi[i][rowidx],Qrelolohi[i][rowidx],
                 Qrehihilo[i][rowidx],Qrelohilo[i][rowidx],
                 Qrehilolo[i][rowidx],Qrelololo[i][rowidx],
                 Qimhihihi[i][colidx],Qimlohihi[i][colidx],
                 Qimhilohi[i][colidx],Qimlolohi[i][colidx],
                 Qimhihilo[i][colidx],Qimlohilo[i][colidx],
                 Qimhilolo[i][colidx],Qimlololo[i][colidx],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&Qsumimhihihi,&Qsumimlohihi,&Qsumimhilohi,&Qsumimlolohi,
                 &Qsumimhihilo,&Qsumimlohilo,&Qsumimhilolo,&Qsumimlololo,
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
         // multiply Q^H with A
         odf_mul(Qrehihihi[i][rowidx],Qrelohihi[i][rowidx],
                 Qrehilohi[i][rowidx],Qrelolohi[i][rowidx],
                 Qrehihilo[i][rowidx],Qrelohilo[i][rowidx],
                 Qrehilolo[i][rowidx],Qrelololo[i][rowidx],
                 Arehihihi[i][colidx],Arelohihi[i][colidx],
                 Arehilohi[i][colidx],Arelolohi[i][colidx],
                 Arehihilo[i][colidx],Arelohilo[i][colidx],
                 Arehilolo[i][colidx],Arelololo[i][colidx],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&Rsumrehihihi,&Rsumrelohihi,&Rsumrehilohi,&Rsumrelolohi,
                 &Rsumrehihilo,&Rsumrelohilo,&Rsumrehilolo,&Rsumrelololo,
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
         odf_mul(Qimhihihi[i][rowidx],Qimlohihi[i][rowidx],
                 Qimhilohi[i][rowidx],Qimlolohi[i][rowidx],
                 Qimhihilo[i][rowidx],Qimlohilo[i][rowidx],
                 Qimhilolo[i][rowidx],Qimlololo[i][rowidx],
                 Aimhihihi[i][colidx],Aimlohihi[i][colidx],
                 Aimhilohi[i][colidx],Aimlolohi[i][colidx],
                 Aimhihilo[i][colidx],Aimlohilo[i][colidx],
                 Aimhilolo[i][colidx],Aimlololo[i][colidx],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&Rsumrehihihi,&Rsumrelohihi,&Rsumrehilohi,&Rsumrelolohi,
                 &Rsumrehihilo,&Rsumrelohilo,&Rsumrehilolo,&Rsumrelololo,
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
         odf_mul(Qimhihihi[i][rowidx],Qimlohihi[i][rowidx],
                 Qimhilohi[i][rowidx],Qimlolohi[i][rowidx],
                 Qimhihilo[i][rowidx],Qimlohilo[i][rowidx],
                 Qimhilolo[i][rowidx],Qimlololo[i][rowidx],
                 Arehihihi[i][colidx],Arelohihi[i][colidx],
                 Arehilohi[i][colidx],Arelolohi[i][colidx],
                 Arehihilo[i][colidx],Arelohilo[i][colidx],
                 Arehilolo[i][colidx],Arelololo[i][colidx],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_dec(&Rsumimhihihi,&Rsumimlohihi,&Rsumimhilohi,&Rsumimlolohi,
                 &Rsumimhihilo,&Rsumimlohilo,&Rsumimhilolo,&Rsumimlololo,
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
         odf_mul(Qrehihihi[i][rowidx],Qrelohihi[i][rowidx],
                 Qrehilohi[i][rowidx],Qrelolohi[i][rowidx],
                 Qrehihilo[i][rowidx],Qrelohilo[i][rowidx],
                 Qrehilolo[i][rowidx],Qrelololo[i][rowidx],
                 Aimhihihi[i][colidx],Aimlohihi[i][colidx],
                 Aimhilohi[i][colidx],Aimlolohi[i][colidx],
                 Aimhihilo[i][colidx],Aimlohilo[i][colidx],
                 Aimhilolo[i][colidx],Aimlololo[i][colidx],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_inc(&Rsumimhihihi,&Rsumimlohihi,&Rsumimhilohi,&Rsumimlolohi,
                 &Rsumimhihilo,&Rsumimlohilo,&Rsumimhilolo,&Rsumimlololo,
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
      }
      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "]re : "
              << Qsumrehihihi << "  " << Qsumrelohihi << endl
              << "                "
              << Qsumrehilohi << "  " << Qsumrelolohi << endl
              << "                "
              << Qsumrehihilo << "  " << Qsumrelohilo << endl
              << "                "
              << Qsumrehilolo << "  " << Qsumrelololo << endl;
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "]im : "
              << Qsumimhihihi << "  " << Qsumimlohihi << endl
              << "                "
              << Qsumimhilohi << "  " << Qsumimlolohi << endl
              << "                "
              << Qsumimhihilo << "  " << Qsumimlohilo << endl
              << "                "
              << Qsumimhilolo << "  " << Qsumimlololo << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "]re : "
              << Rsumrehihihi << "  " << Rsumrelohihi << endl
              << "                "
              << Rsumrehilohi << "  " << Rsumrelolohi << endl
              << "                "
              << Rsumrehihilo << "  " << Rsumrelohilo << endl
              << "                "
              << Rsumrehilolo << "  " << Rsumrelololo << endl;
         cout << "    R[" << rowidx << "][" << colidx << "]re : "
              << Rrehihihi[rowidx][colidx] << "  "
              << Rrelohihi[rowidx][colidx] << endl
              << "                "
              << Rrehilohi[rowidx][colidx] << "  "
              << Rrelolohi[rowidx][colidx] << endl
              << "                "
              << Rrehihilo[rowidx][colidx] << "  "
              << Rrelohilo[rowidx][colidx] << endl
              << "                "
              << Rrehilolo[rowidx][colidx] << "  "
              << Rrelololo[rowidx][colidx] << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "]im : "
              << Rsumimhihihi << "  " << Rsumimlohihi << endl
              << "                "
              << Rsumimhilohi << "  " << Rsumimlolohi << endl
              << "                "
              << Rsumimhihilo << "  " << Rsumimlohilo << endl
              << "                "
              << Rsumimhilolo << "  " << Rsumimlololo << endl;
         cout << "    R[" << rowidx << "][" << colidx << "]im : "
              << Rimhihihi[rowidx][colidx] << "  "
              << Rimlohihi[rowidx][colidx] << endl
              << "                "
              << Rimhilohi[rowidx][colidx] << "  "
              << Rimlolohi[rowidx][colidx] << endl
              << "                "
              << Rimhihilo[rowidx][colidx] << "  "
              << Rimlohilo[rowidx][colidx] << endl
              << "                "
              << Rimhilolo[rowidx][colidx] << "  "
              << Rimlololo[rowidx][colidx] << endl;
      }
      if(rowidx == colidx)
      {
         errorQ = errorQ + fabs(Qsumrehihihi - 1.0) + fabs(Qsumrelohihi)
                         + fabs(Qsumrehilohi) + fabs(Qsumrelolohi)
                         + fabs(Qsumimhihihi) + fabs(Qsumimlohihi)
                         + fabs(Qsumimhilohi) + fabs(Qsumimlolohi)
                         + fabs(Qsumrehihilo) + fabs(Qsumrelohilo)
                         + fabs(Qsumrehilolo) + fabs(Qsumrelololo)
                         + fabs(Qsumimhihilo) + fabs(Qsumimlohilo)
                         + fabs(Qsumimhilolo) + fabs(Qsumimlololo);
      }
      else
      {
         errorQ = errorQ + fabs(Qsumrehihihi) + fabs(Qsumrelohihi)
                         + fabs(Qsumrehilohi) + fabs(Qsumrelolohi)
                         + fabs(Qsumimhihihi) + fabs(Qsumimlohihi)
                         + fabs(Qsumimhilohi) + fabs(Qsumimlolohi)
                         + fabs(Qsumrehihilo) + fabs(Qsumrelohilo)
                         + fabs(Qsumrehilolo) + fabs(Qsumrelololo)
                         + fabs(Qsumimhihilo) + fabs(Qsumimlohilo)
                         + fabs(Qsumimhilolo) + fabs(Qsumimlololo);
      }
      errorR = errorR + fabs(Rsumrehihihi - Rrehihihi[rowidx][colidx])
                      + fabs(Rsumrehilohi - Rrehilohi[rowidx][colidx])
                      + fabs(Rsumrelohihi - Rrelohihi[rowidx][colidx])
                      + fabs(Rsumrelolohi - Rrelolohi[rowidx][colidx])
                      + fabs(Rsumimhihilo - Rimhihilo[rowidx][colidx])
                      + fabs(Rsumimhilolo - Rimhilolo[rowidx][colidx])
                      + fabs(Rsumimlohilo - Rimlohilo[rowidx][colidx])
                      + fabs(Rsumimlololo - Rimlololo[rowidx][colidx])
                      + fabs(Rsumrehihihi - Rrehihihi[rowidx][colidx])
                      + fabs(Rsumrehilohi - Rrehilohi[rowidx][colidx])
                      + fabs(Rsumrelohihi - Rrelohihi[rowidx][colidx])
                      + fabs(Rsumrelolohi - Rrelolohi[rowidx][colidx])
                      + fabs(Rsumimhihilo - Rimhihilo[rowidx][colidx])
                      + fabs(Rsumimhilolo - Rimhilolo[rowidx][colidx])
                      + fabs(Rsumimlohilo - Rimlohilo[rowidx][colidx])
                      + fabs(Rsumimlololo - Rimlololo[rowidx][colidx]);
   }
   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   return int(errorQ + errorR > tol);
}

void test_factors_real8_houseqr ( void )
{
   cout << "Give the number of rows : ";
   int nrows; cin >> nrows;

   cout << "Give the number of columns : ";
   int ncols; cin >> ncols;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Ahihihi = new double*[nrows];
   double **Alohihi = new double*[nrows];
   double **Ahilohi = new double*[nrows];
   double **Alolohi = new double*[nrows];
   double **Ahihilo = new double*[nrows];
   double **Alohilo = new double*[nrows];
   double **Ahilolo = new double*[nrows];
   double **Alololo = new double*[nrows];
   double **Qhihihi = new double*[nrows];
   double **Qlohihi = new double*[nrows];
   double **Qhilohi = new double*[nrows];
   double **Qlolohi = new double*[nrows];
   double **Qhihilo = new double*[nrows];
   double **Qlohilo = new double*[nrows];
   double **Qhilolo = new double*[nrows];
   double **Qlololo = new double*[nrows];
   double **Rhihihi = new double*[nrows];
   double **Rlohihi = new double*[nrows];
   double **Rhilohi = new double*[nrows];
   double **Rlolohi = new double*[nrows];
   double **Rhihilo = new double*[nrows];
   double **Rlohilo = new double*[nrows];
   double **Rhilolo = new double*[nrows];
   double **Rlololo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahihihi[i] = new double[ncols];
      Alohihi[i] = new double[ncols];
      Ahilohi[i] = new double[ncols];
      Alolohi[i] = new double[ncols];
      Ahihilo[i] = new double[ncols];
      Alohilo[i] = new double[ncols];
      Ahilolo[i] = new double[ncols];
      Alololo[i] = new double[ncols];
      Qhihihi[i] = new double[nrows];
      Qlohihi[i] = new double[nrows];
      Qhilohi[i] = new double[nrows];
      Qlolohi[i] = new double[nrows];
      Qhihilo[i] = new double[nrows];
      Qlohilo[i] = new double[nrows];
      Qhilolo[i] = new double[nrows];
      Qlololo[i] = new double[nrows];
      Rhihihi[i] = new double[ncols];
      Rlohihi[i] = new double[ncols];
      Rhilohi[i] = new double[ncols];
      Rlolohi[i] = new double[ncols];
      Rhihilo[i] = new double[ncols];
      Rlohilo[i] = new double[ncols];
      Rhilolo[i] = new double[ncols];
      Rlololo[i] = new double[ncols];
   }
   random_dbl8_matrix
      (nrows,ncols,Ahihihi,Alohihi,Ahilohi,Alolohi,
                   Ahihilo,Alohilo,Ahilolo,Alololo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihihi[i][j] << "  " << Alohihi[i][j] << endl
                 << "          "
                 << Ahilohi[i][j] << "  " << Alolohi[i][j] << endl
                 << "          "
                 << Ahihilo[i][j] << "  " << Alohilo[i][j] << endl
                 << "          "
                 << Ahilolo[i][j] << "  " << Alololo[i][j] << endl;
   }
   CPU_dbl8_factors_houseqr
      (nrows,ncols,Ahihihi,Alohihi,Ahilohi,Alolohi,
                   Ahihilo,Alohilo,Ahilolo,Alololo,
                   Qhihihi,Qlohihi,Qhilohi,Qlolohi,
                   Qhihilo,Qlohilo,Qhilolo,Qlololo,
                   Rhihihi,Rlohihi,Rhilohi,Rlolohi,
                   Rhihilo,Rlohilo,Rhilolo,Rlololo);

   const double tol = 1.0E-80;
   const int fail = test_real8_qr_factors
      (nrows,ncols,Ahihihi,Alohihi,Ahilohi,Alolohi,
                   Ahihilo,Alohilo,Ahilolo,Alololo,
                   Qhihihi,Qlohihi,Qhilohi,Qlolohi,
                   Qhihilo,Qlohilo,Qhilolo,Qlololo,
                   Rhihihi,Rlohihi,Rhilohi,Rlolohi,
                   Rhihilo,Rlohilo,Rhilolo,Rlololo,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(Ahihihi[i]); free(Qhihihi[i]); free(Rhihihi[i]);
      free(Alohihi[i]); free(Qlohihi[i]); free(Rlohihi[i]);
      free(Ahilohi[i]); free(Qhilohi[i]); free(Rhilohi[i]);
      free(Alolohi[i]); free(Qlolohi[i]); free(Rlolohi[i]);
      free(Ahihilo[i]); free(Qhihilo[i]); free(Rhihilo[i]);
      free(Alohilo[i]); free(Qlohilo[i]); free(Rlohilo[i]);
      free(Ahilolo[i]); free(Qhilolo[i]); free(Rhilolo[i]);
      free(Alololo[i]); free(Qlololo[i]); free(Rlololo[i]);
   }
   free(Ahihihi); free(Qhihihi); free(Rhihihi);
   free(Alohihi); free(Qlohihi); free(Rlohihi);
   free(Ahilohi); free(Qhilohi); free(Rhilohi);
   free(Alolohi); free(Qlolohi); free(Rlolohi);
   free(Ahihilo); free(Qhihilo); free(Rhihilo);
   free(Alohilo); free(Qlohilo); free(Rlohilo);
   free(Ahilolo); free(Qhilolo); free(Rhilolo);
   free(Alololo); free(Qlololo); free(Rlololo);
}

void test_factors_cmplx8_houseqr ( void )
{
   cout << "Give the number of rows : ";
   int nrows; cin >> nrows;

   cout << "Give the number of columns : ";
   int ncols; cin >> ncols;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Arehihihi = new double*[nrows];
   double **Arelohihi = new double*[nrows];
   double **Arehilohi = new double*[nrows];
   double **Arelolohi = new double*[nrows];
   double **Arehihilo = new double*[nrows];
   double **Arelohilo = new double*[nrows];
   double **Arehilolo = new double*[nrows];
   double **Arelololo = new double*[nrows];
   double **Aimhihihi = new double*[nrows];
   double **Aimlohihi = new double*[nrows];
   double **Aimhilohi = new double*[nrows];
   double **Aimlolohi = new double*[nrows];
   double **Aimhihilo = new double*[nrows];
   double **Aimlohilo = new double*[nrows];
   double **Aimhilolo = new double*[nrows];
   double **Aimlololo = new double*[nrows];
   double **Qrehihihi = new double*[nrows];
   double **Qrelohihi = new double*[nrows];
   double **Qrehilohi = new double*[nrows];
   double **Qrelolohi = new double*[nrows];
   double **Qrehihilo = new double*[nrows];
   double **Qrelohilo = new double*[nrows];
   double **Qrehilolo = new double*[nrows];
   double **Qrelololo = new double*[nrows];
   double **Qimhihihi = new double*[nrows];
   double **Qimlohihi = new double*[nrows];
   double **Qimhilohi = new double*[nrows];
   double **Qimlolohi = new double*[nrows];
   double **Qimhihilo = new double*[nrows];
   double **Qimlohilo = new double*[nrows];
   double **Qimhilolo = new double*[nrows];
   double **Qimlololo = new double*[nrows];
   double **Rrehihihi = new double*[nrows];
   double **Rrelohihi = new double*[nrows];
   double **Rrehilohi = new double*[nrows];
   double **Rrelolohi = new double*[nrows];
   double **Rrehihilo = new double*[nrows];
   double **Rrelohilo = new double*[nrows];
   double **Rrehilolo = new double*[nrows];
   double **Rrelololo = new double*[nrows];
   double **Rimhihihi = new double*[nrows];
   double **Rimlohihi = new double*[nrows];
   double **Rimhilohi = new double*[nrows];
   double **Rimlolohi = new double*[nrows];
   double **Rimhihilo = new double*[nrows];
   double **Rimlohilo = new double*[nrows];
   double **Rimhilolo = new double*[nrows];
   double **Rimlololo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Arehihihi[i] = new double[ncols];
      Arelohihi[i] = new double[ncols];
      Arehilohi[i] = new double[ncols];
      Arelolohi[i] = new double[ncols];
      Arehihilo[i] = new double[ncols];
      Arelohilo[i] = new double[ncols];
      Arehilolo[i] = new double[ncols];
      Arelololo[i] = new double[ncols];
      Aimhihihi[i] = new double[ncols];
      Aimlohihi[i] = new double[ncols];
      Aimhilohi[i] = new double[ncols];
      Aimlolohi[i] = new double[ncols];
      Aimhihilo[i] = new double[ncols];
      Aimlohilo[i] = new double[ncols];
      Aimhilolo[i] = new double[ncols];
      Aimlololo[i] = new double[ncols];
      Qrehihihi[i] = new double[nrows];
      Qrelohihi[i] = new double[nrows];
      Qrehilohi[i] = new double[nrows];
      Qrelolohi[i] = new double[nrows];
      Qrehihilo[i] = new double[nrows];
      Qrelohilo[i] = new double[nrows];
      Qrehilolo[i] = new double[nrows];
      Qrelololo[i] = new double[nrows];
      Qimhihihi[i] = new double[nrows];
      Qimlohihi[i] = new double[nrows];
      Qimhilohi[i] = new double[nrows];
      Qimlolohi[i] = new double[nrows];
      Qimhihilo[i] = new double[nrows];
      Qimlohilo[i] = new double[nrows];
      Qimhilolo[i] = new double[nrows];
      Qimlololo[i] = new double[nrows];
      Rrehihihi[i] = new double[ncols];
      Rrelohihi[i] = new double[ncols];
      Rrehilohi[i] = new double[ncols];
      Rrelolohi[i] = new double[ncols];
      Rrehihilo[i] = new double[ncols];
      Rrelohilo[i] = new double[ncols];
      Rrehilolo[i] = new double[ncols];
      Rrelololo[i] = new double[ncols];
      Rimhihihi[i] = new double[ncols];
      Rimlohihi[i] = new double[ncols];
      Rimhilohi[i] = new double[ncols];
      Rimlolohi[i] = new double[ncols];
      Rimhihilo[i] = new double[ncols];
      Rimlohilo[i] = new double[ncols];
      Rimhilolo[i] = new double[ncols];
      Rimlololo[i] = new double[ncols];
   }
   random_cmplx8_matrix
      (nrows,ncols,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
                   Arehihilo,Arelohilo,Arehilolo,Arelololo,
                   Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
                   Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
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
   CPU_cmplx8_factors_houseqr
      (nrows,ncols,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
                   Arehihilo,Arelohilo,Arehilolo,Arelololo,
                   Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
                   Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo,
                   Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
                   Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
                   Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
                   Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
                   Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
                   Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
                   Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
                   Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo);

   const double tol = 1.0e-80;
   const int fail = test_cmplx8_qr_factors
      (nrows,ncols,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
                   Arehihilo,Arelohilo,Arehilolo,Arelololo,
                   Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
                   Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo,
                   Qrehihihi,Qrelohihi,Qrehilohi,Qrelolohi,
                   Qrehihilo,Qrelohilo,Qrehilolo,Qrelololo,
                   Qimhihihi,Qimlohihi,Qimhilohi,Qimlolohi,
                   Qimhihilo,Qimlohilo,Qimhilolo,Qimlololo,
                   Rrehihihi,Rrelohihi,Rrehilohi,Rrelolohi,
                   Rrehihilo,Rrelohilo,Rrehilolo,Rrelololo,
                   Rimhihihi,Rimlohihi,Rimhilohi,Rimlolohi,
                   Rimhihilo,Rimlohilo,Rimhilolo,Rimlololo,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(Arehihihi[i]); free(Qrehihihi[i]); free(Rrehihihi[i]);
      free(Arelohihi[i]); free(Qrelohihi[i]); free(Rrelohihi[i]);
      free(Arehilohi[i]); free(Qrehilohi[i]); free(Rrehilohi[i]);
      free(Arelolohi[i]); free(Qrelolohi[i]); free(Rrelolohi[i]);
      free(Arehihilo[i]); free(Qrehihilo[i]); free(Rrehihilo[i]);
      free(Arelohilo[i]); free(Qrelohilo[i]); free(Rrelohilo[i]);
      free(Arehilolo[i]); free(Qrehilolo[i]); free(Rrehilolo[i]);
      free(Arelololo[i]); free(Qrelololo[i]); free(Rrelololo[i]);
      free(Aimhihihi[i]); free(Qimhihihi[i]); free(Rimhihihi[i]);
      free(Aimlohihi[i]); free(Qimlohihi[i]); free(Rimlohihi[i]);
      free(Aimhilohi[i]); free(Qimhilohi[i]); free(Rimhilohi[i]);
      free(Aimlolohi[i]); free(Qimlolohi[i]); free(Rimlolohi[i]);
      free(Aimhihilo[i]); free(Qimhihilo[i]); free(Rimhihilo[i]);
      free(Aimlohilo[i]); free(Qimlohilo[i]); free(Rimlohilo[i]);
      free(Aimhilolo[i]); free(Qimhilolo[i]); free(Rimhilolo[i]);
      free(Aimlololo[i]); free(Qimlololo[i]); free(Rimlololo[i]);
   }
   free(Arehihihi); free(Qrehihihi); free(Rrehihihi);
   free(Arelohihi); free(Qrelohihi); free(Rrelohihi);
   free(Arehilohi); free(Qrehilohi); free(Rrehilohi);
   free(Arelolohi); free(Qrelolohi); free(Rrelolohi);
   free(Arehihilo); free(Qrehihilo); free(Rrehihilo);
   free(Arelohilo); free(Qrelohilo); free(Rrelohilo);
   free(Arehilolo); free(Qrehilolo); free(Rrehilolo);
   free(Arelololo); free(Qrelololo); free(Rrelololo);
   free(Aimhihihi); free(Qimhihihi); free(Rimhihihi);
   free(Aimlohihi); free(Qimlohihi); free(Rimlohihi);
   free(Aimhilohi); free(Qimhilohi); free(Rimhilohi);
   free(Aimlolohi); free(Qimlolohi); free(Rimlolohi);
   free(Aimhihilo); free(Qimhihilo); free(Rimhihilo);
   free(Aimlohilo); free(Qimlohilo); free(Rimlohilo);
   free(Aimhilolo); free(Qimhilolo); free(Rimhilolo);
   free(Aimlololo); free(Qimlololo); free(Rimlololo);
}
