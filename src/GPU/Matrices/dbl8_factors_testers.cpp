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
      rhsimhihilo[i] = 0.0; rhsimlohilo[i] = 0.0;
      rhsimhilolo[i] = 0.0; rhsimlololo[i] = 0.0;
      rhsrehihihi[i] = 0.0; rhsrelohihi[i] = 0.0;
      rhsrehilohi[i] = 0.0; rhsrelolohi[i] = 0.0;
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
   return 0;
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
   return 0;
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
   return 0;
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
   return 0;
}

void test_factors_real8_houseqr ( void )
{
}

void test_factors_cmplx8_houseqr ( void )
{
}
