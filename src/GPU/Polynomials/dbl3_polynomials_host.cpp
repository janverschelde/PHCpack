/* The file dbl3_polynomials_host.cpp defines functions specified
 * in dbl3_polynomials_host.h. */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "triple_double_functions.h"
#include "dbl3_convolutions_host.h"
#include "dbl3_monomials_host.h"
#include "dbl3_polynomials_host.h"

void CPU_dbl3_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
   double **forwardhi, double **forwardmi, double **forwardlo,
   double **backwardhi, double **backwardmi, double **backwardlo,
   double **crosshi, double **crossmi, double **crosslo, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_dbl3_product(deg,inputhi[ix1],inputmi[ix1],inputlo[ix1],
                              cffhi[i],cffmi[i],cfflo[i],
                              forwardhi[0],forwardmi[0],forwardlo[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            // output[dim][j] += forward[0][j];
            tdf_inc(&outputhi[dim][j],&outputmi[dim][j],&outputlo[dim][j],
                    forwardhi[0][j],forwardmi[0][j],forwardlo[0][j]);
            // output[ix1][j] += cff[i][j];
            tdf_inc(&outputhi[ix1][j],&outputmi[ix1][j],&outputlo[ix1][j],
                    cffhi[i][j],cffmi[i][j],cfflo[i][j]);
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_dbl3_product(deg,cffhi[i],cffmi[i],cfflo[i],
                              inputhi[ix1],inputmi[ix1],inputlo[ix1],
                              forwardhi[0],forwardmi[0],forwardlo[0]);
         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[0][j];
            tdf_inc(&outputhi[ix2][j],&outputmi[ix2][j],&outputlo[ix2][j],
                    forwardhi[0][j],forwardmi[0][j],forwardlo[0][j]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_dbl3_product(deg,cffhi[i],cffmi[i],cfflo[i],
                              inputhi[ix2],inputmi[ix2],inputlo[ix2],
                              backwardhi[0],backwardmi[0],backwardlo[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) // output[ix1][j] += backward[0][j];
            tdf_inc(&outputhi[ix1][j],&outputmi[ix1][j],&outputlo[ix1][j],
                    backwardhi[0][j],backwardmi[0][j],backwardlo[0][j]);

         CPU_dbl3_product(deg,forwardhi[0],forwardmi[0],forwardlo[0],
                              inputhi[ix2],inputmi[ix2],inputlo[ix2],
                              forwardhi[1],forwardmi[1],forwardlo[1]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[1][j];
            tdf_inc(&outputhi[dim][j],&outputmi[dim][j],&outputlo[dim][j],
                    forwardhi[1][j],forwardmi[1][j],forwardlo[1][j]);
      }
      else if(nvr[i] > 2)
      {
         CPU_dbl3_speel
            (nvr[i],deg,idx[i],cffhi[i],cffmi[i],cfflo[i],
             inputhi,inputmi,inputlo,forwardhi,forwardmi,forwardlo,
             backwardhi,backwardmi,backwardlo,crosshi,crossmi,crosslo);

         ix1 = nvr[i]-1;               // update the value of the polynomial
         for(int j=0; j<=deg; j++) // output[dim][j] += forward[ix1][j];
            tdf_inc(&outputhi[dim][j],&outputmi[dim][j],&outputlo[dim][j],
                    forwardhi[ix1][j],forwardmi[ix1][j],forwardlo[ix1][j]);

         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) // output[ix2][j] += forward[ix1][j];
            tdf_inc(&outputhi[ix2][j],&outputmi[ix2][j],&outputlo[ix2][j],
                    forwardhi[ix1][j],forwardmi[ix1][j],forwardlo[ix1][j]);

         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) // output[ix2][j] += backward[ix1][j];
            tdf_inc(&outputhi[ix2][j],&outputmi[ix2][j],&outputlo[ix2][j],
                    backwardhi[ix1][j],backwardmi[ix1][j],backwardlo[ix1][j]);

         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) // output[ix2][j] += cross[k-1][j];
               tdf_inc(&outputhi[ix2][j],&outputmi[ix2][j],&outputlo[ix2][j],
                       crosshi[k-1][j],crossmi[k-1][j],crosslo[k-1][j]);
         }
      }
   }
}

void CPU_dbl3_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi, double **inputlo, 
   double **outputhi, double **outputmi, double **outputlo,
   double *elapsedsec, bool verbose )
{
   double **forwardhi = new double*[dim];
   double **forwardmi = new double*[dim];
   double **forwardlo = new double*[dim];
   double **backwardhi = new double*[dim-1]; // in case dim = 2
   double **backwardmi = new double*[dim-1];
   double **backwardlo = new double*[dim-1];
   double **crosshi = new double*[dim-1];    // in case dim = 2
   double **crossmi = new double*[dim-1];
   double **crosslo = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardhi[i] = new double[deg+1]; forwardmi[i] = new double[deg+1];
      forwardlo[i] = new double[deg+1];
      backwardhi[i] = new double[deg+1]; backwardmi[i] = new double[deg+1];
      backwardlo[i] = new double[deg+1];
      crosshi[i] = new double[deg+1]; crossmi[i] = new double[deg+1];
      crosslo[i] = new double[deg+1];
   }
   forwardhi[dim-1] = new double[deg+1];
   forwardmi[dim-1] = new double[deg+1];
   forwardlo[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputhi[dim][i] = csthi[i];
      outputmi[dim][i] = cstmi[i];
      outputlo[dim][i] = cstlo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputhi[i][j] = 0.0;
         outputmi[i][j] = 0.0;
         outputlo[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_dbl3_poly_speel
      (dim,nbr,deg,nvr,idx,cffhi,cffmi,cfflo,inputhi,inputmi,inputlo,
       outputhi,outputmi,outputlo,forwardhi,forwardmi,forwardlo,
       backwardhi,backwardmi,backwardlo,crosshi,crossmi,crosslo,verbose);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(verbose)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
   }
   for(int i=0; i<dim-1; i++)
   {
      free(forwardhi[i]); free(backwardhi[i]); free(crosshi[i]);
      free(forwardmi[i]); free(backwardmi[i]); free(crossmi[i]);
      free(forwardlo[i]); free(backwardlo[i]); free(crosslo[i]);
   }
   free(forwardhi[dim-1]); free(forwardmi[dim-1]); free(forwardlo[dim-1]);
   free(forwardhi); free(backwardhi); free(crosshi);
   free(forwardmi); free(backwardmi); free(crossmi);
   free(forwardlo); free(backwardlo); free(crosslo);
}

void CPU_dbl3_conv_job
 ( int deg, int nvr, int *idx,
   double *cffhi, double *cffmi, double *cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **forwardhi, double **forwardmi, double **forwardlo,
   double **backwardhi, double **backwardmi, double **backwardlo,
   double **crosshi, double **crossmi, double **crosslo,
   ConvolutionJob job, bool verbose )
{
   const int inp1tp = job.get_first_type();
   const int inp1ix = job.get_first_input();
   const int inp2tp = job.get_second_type();
   const int inp2ix = job.get_second_input();
   const int outptp = job.get_output_type();
   const int outidx = job.get_output_index();

   if(outptp == 1) // forward product either initializes or accumulates
   {
      if(verbose) cout << "-> computing f[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_dbl3_product(deg,
            cffhi,cffmi,cfflo,
            inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
            forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl3_product(deg,
               inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
               cffhi,cffmi,cfflo,
               forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
               forwardhi[inp2ix],forwardmi[inp2ix],forwardlo[inp2ix],
               forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_dbl3_product(deg,
            crosshi[inp1ix],crossmi[inp1ix],crosslo[inp1ix],
            inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
            forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl3_product(deg,
               inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
               cffhi,cffmi,cfflo,
               forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               forwardhi[inp1ix],forwardmi[inp1ix],forwardlo[inp1ix],
               inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
               forwardhi[outidx],forwardmi[outidx],forwardlo[outidx]);
         }
      }
   }
   else if(outptp == 2) // backward product either initializes or accumulates
   {
      if(verbose) cout << "-> computing b[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               cffhi,cffmi,cfflo,
               inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               cffhi,cffmi,cfflo,
               backwardhi[inp2ix],backwardmi[inp2ix],backwardlo[inp2ix],
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_dbl3_product(deg,
               inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
               inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
               backwardhi[inp2ix],backwardmi[inp2ix],backwardlo[inp2ix],
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_dbl3_product(deg,
               backwardhi[inp1ix],backwardmi[inp1ix],backwardlo[inp1ix],
               cffhi,cffmi,cfflo,
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl3_product(deg,
               backwardhi[inp1ix],backwardmi[inp1ix],backwardlo[inp1ix],
               inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
               backwardhi[outidx],backwardmi[outidx],backwardlo[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_dbl3_product(deg,
            cffhi,cffmi,cfflo,
            inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
            crosshi[outidx],crossmi[outidx],crosslo[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl3_product(deg,
            inputhi[inp1ix],inputmi[inp1ix],inputlo[inp1ix],
            forwardhi[inp2ix],forwardmi[inp2ix],forwardlo[inp2ix],
            crosshi[outidx],crossmi[outidx],crosslo[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_dbl3_product(deg,
              forwardhi[inp1ix],forwardmi[inp1ix],forwardlo[inp1ix],
              inputhi[inp2ix],inputmi[inp2ix],inputlo[inp2ix],
              crosshi[outidx],crossmi[outidx],crosslo[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_dbl3_product(deg,
              forwardhi[inp1ix],forwardmi[inp1ix],forwardlo[inp1ix],
              backwardhi[inp2ix],backwardmi[inp2ix],backwardlo[inp2ix],
              crosshi[outidx],crossmi[outidx],crosslo[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl3_product(deg,
            backwardhi[inp1ix],backwardmi[inp1ix],backwardlo[inp1ix],
            forwardhi[inp2ix],forwardmi[inp2ix],forwardlo[inp2ix],
            crosshi[outidx],crossmi[outidx],crosslo[outidx]);
      }
   }
}

void CPU_dbl3_add_job
 ( int deg, double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double ***forwardhi, double ***forwardmi, double ***forwardlo,
   double ***backwardhi, double ***backwardmi,  double ***backwardlo, 
   double ***crosshi, double ***crossmi, double ***crosslo,
   AdditionJob job, bool verbose )
{
   const int adtype = job.get_addition_type();
   const int intype = job.get_increment_type();
   const int updmon = job.get_update_monomial();
   const int updidx = job.get_update_index();
   const int incmon = job.get_increment_monomial();
   const int incidx = job.get_increment_index();

   if(adtype == 1)
   {
      if(incmon < 0)
      {
         if(incidx < 0)
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cst[i];
               tdf_inc(&forwardhi[updmon][updidx][i],
                       &forwardmi[updmon][updidx][i],
                       &forwardlo[updmon][updidx][i],
                       csthi[i],cstmi[i],cstlo[i]);
         else
            for(int i=0; i<=deg; i++)
               // forward[updmon][updidx][i] += cff[incidx][i];
               tdf_inc(&forwardhi[updmon][updidx][i],
                       &forwardmi[updmon][updidx][i],
                       &forwardlo[updmon][updidx][i],
                       cffhi[incidx][i],cffmi[incidx][i],cfflo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += forward[incmon][incidx][i];
            tdf_inc(&forwardhi[updmon][updidx][i],
                    &forwardmi[updmon][updidx][i],
                    &forwardlo[updmon][updidx][i],
                    forwardhi[incmon][incidx][i],
                    forwardmi[incmon][incidx][i],
                    forwardlo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += backward[incmon][incidx][i];
            tdf_inc(&forwardhi[updmon][updidx][i],
                    &forwardmi[updmon][updidx][i],
                    &forwardlo[updmon][updidx][i],
                    backwardhi[incmon][incidx][i],
                    backwardmi[incmon][incidx][i],
                    backwardlo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // forward[updmon][updidx][i] += cross[incmon][incidx][i];
            tdf_inc(&forwardhi[updmon][updidx][i],
                    &forwardmi[updmon][updidx][i],
                    &forwardlo[updmon][updidx][i],
                    crosshi[incmon][incidx][i],
                    crossmi[incmon][incidx][i],
                    crosslo[incmon][incidx][i]);
      }
   }
   else if(adtype == 2)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cff[incidx][i];
            tdf_inc(&backwardhi[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardlo[updmon][updidx][i],
                    cffhi[incidx][i],cffmi[incidx][i],cfflo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += forward[incmon][incidx][i];
            tdf_inc(&backwardhi[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardlo[updmon][updidx][i],
                    forwardhi[incmon][incidx][i],
                    forwardmi[incmon][incidx][i],
                    forwardlo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += backward[incmon][incidx][i];
            tdf_inc(&backwardhi[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardlo[updmon][updidx][i],
                    backwardhi[incmon][incidx][i],
                    backwardmi[incmon][incidx][i],
                    backwardlo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // backward[updmon][updidx][i] += cross[incmon][incidx][i];
            tdf_inc(&backwardhi[updmon][updidx][i],
                    &backwardmi[updmon][updidx][i],
                    &backwardlo[updmon][updidx][i],
                    crosshi[incmon][incidx][i],
                    crossmi[incmon][incidx][i],
                    crosslo[incmon][incidx][i]);
      }
   }
   else if(adtype == 3)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cff[incidx][i];
            tdf_inc(&crosshi[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crosslo[updmon][updidx][i],
                    cffhi[incidx][i],cffmi[incidx][i],cfflo[incidx][i]);
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += forward[incmon][incidx][i];
            tdf_inc(&crosshi[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crosslo[updmon][updidx][i],
                    forwardhi[incmon][incidx][i],
                    forwardmi[incmon][incidx][i],
                    forwardlo[incmon][incidx][i]);
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += backward[incmon][incidx][i];
            tdf_inc(&crosshi[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crosslo[updmon][updidx][i],
                    backwardhi[incmon][incidx][i],
                    backwardmi[incmon][incidx][i],
                    backwardlo[incmon][incidx][i]);
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            // cross[updmon][updidx][i] += cross[incmon][incidx][i];
            tdf_inc(&crosshi[updmon][updidx][i],
                    &crossmi[updmon][updidx][i],
                    &crosslo[updmon][updidx][i],
                    crosshi[incmon][incidx][i],
                    crossmi[incmon][incidx][i],
                    crosslo[incmon][incidx][i]);
      }
   }
}

void CPU_dbl3_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
   double ***forwardhi, double ***forwardmi, double ***forwardlo,
   double ***backwardhi, double ***backwardmi, double ***backwardlo,
   double ***crosshi, double ***crossmi, double ***crosslo )
{
   for(int i=0; i<=deg; i++)
   {
      outputhi[dim][i] = csthi[i];
      outputmi[dim][i] = cstmi[i];
      outputlo[dim][i] = cstlo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputhi[i][j] = 0.0;
         outputmi[i][j] = 0.0;
         outputlo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      int ix0 = idx[k][0];   // first variable in monomial k
      int ix1 = nvr[k]-1;    // last forward has the value
      int ix2 = nvr[k]-2;    // next to last forward has last derivative
                             // last backward has the first derivative
      int ixn = idx[k][ix1]; // index of the last variable in monomial k

      for(int i=0; i<=deg; i++) // value is last forward location
         // output[dim][i] = output[dim][i] + forward[k][ix1][i];
         tdf_inc(&outputhi[dim][i],&outputmi[dim][i],&outputlo[dim][i],
                 forwardhi[k][ix1][i],forwardmi[k][ix1][i],
                 forwardlo[k][ix1][i]);

      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            // output[ix0][i] = output[ix0][i] + cff[k][i]; 
            tdf_inc(&outputhi[ix0][i],&outputmi[ix0][i],&outputlo[ix0][i],
                    cffhi[k][i],cffmi[k][i],cfflo[k][i]);
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            // output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            tdf_inc(&outputhi[ixn][i],&outputmi[ixn][i],&outputlo[ixn][i],
                    forwardhi[k][ix2][i],forwardmi[k][ix2][i],
                    forwardlo[k][ix2][i]);
            // output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
            tdf_inc(&outputhi[ix0][i],&outputmi[ix0][i],&outputlo[ix0][i],
                    backwardhi[k][ix2][i],backwardmi[k][ix2][i],
                    backwardlo[k][ix2][i]);
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  // output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
                  tdf_inc(&outputhi[ix0][i],&outputmi[ix0][i],
                          &outputlo[ix0][i],crosshi[k][j-1][i],
                          crossmi[k][j-1][i],crosslo[k][j-1][i]);
            }
         }
      }
   }
}

void CPU_dbl3_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
   double ***forwardhi, double ***forwardmi, double ***forwardlo,
   double ***backwardhi, double ***backwardmi, double ***backwardlo,
   double ***crosshi, double ***crossmi, double ***crosslo,
   AdditionJobs jobs, bool verbose )
{
   for(int k=0; k<jobs.get_depth(); k++)
   {
      if(verbose) cout << "executing addition jobs at layer "
                       << k << " :" << endl;
      for(int i=0; i<jobs.get_layer_count(k); i++)
      {
         AdditionJob job = jobs.get_job(k,i);
         if(verbose) cout << "job " << i << " : " << job << endl;

         CPU_dbl3_add_job(deg,csthi,cstmi,cstlo,cffhi,cffmi,cfflo,
            forwardhi,forwardmi,forwardlo,backwardhi,backwardmi,backwardlo,
            crosshi,crossmi,crosslo,job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {  // output[dim][i] = forward[lastmon][lastidx][i];
      outputhi[dim][i] = forwardhi[lastmon][lastidx][i];
      outputmi[dim][i] = forwardmi[lastmon][lastidx][i];
      outputlo[dim][i] = forwardlo[lastmon][lastidx][i];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++)
      {
         outputhi[0][i] = 0.0;
         outputmi[0][i] = 0.0;
         outputlo[0][i] = 0.0;
      }
   }
   else
   {
      int ix0 = jobs.get_differential_index(0,cnt);
      int ix2 = nvr[ix0] - 2;
      
      if(verbose)
         cout << "Updating derivative 0, ix0 = " << ix0
              << ", ix2 = " << ix2
              << " : b[" << ix0 << "," << ix2 << "]" << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = backward[ix0][ix2][i];
      {
         outputhi[0][i] = backwardhi[ix0][ix2][i];
         outputmi[0][i] = backwardmi[ix0][ix2][i];
         outputlo[0][i] = backwardlo[ix0][ix2][i];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);
      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         for(int i=0; i<=deg; i++) 
         {
            outputhi[k][i] = 0.0;
            outputmi[k][i] = 0.0;
            outputlo[k][i] = 0.0;
         }
      }
      else
      {
         int ix0 = jobs.get_differential_index(k,cnt);

         if(idx[ix0][0] == k) // k is first variable of monomial
         {
            int ix2 = nvr[ix0] - 2;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : b[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = backward[ix0][ix2][i];
            {
               outputhi[k][i] = backwardhi[ix0][ix2][i];
               outputmi[k][i] = backwardmi[ix0][ix2][i];
               outputlo[k][i] = backwardlo[ix0][ix2][i];
            }
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0] - 2;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : f[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = forward[ix0][ix2][i];
            {
               outputhi[k][i] = forwardhi[ix0][ix2][i];
               outputmi[k][i] = forwardmi[ix0][ix2][i];
               outputlo[k][i] = forwardlo[ix0][ix2][i];
            }
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : c[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = cross[ix0][ix2][i];
            {
               outputhi[k][i] = crosshi[ix0][ix2][i];
               outputmi[k][i] = crossmi[ix0][ix2][i];
               outputlo[k][i] = crosslo[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_dbl3_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, bool verbose )
{
   double ***forwardhi = new double**[nbr];
   double ***forwardmi = new double**[nbr];
   double ***forwardlo = new double**[nbr];
   double ***backwardhi = new double**[nbr];
   double ***backwardmi = new double**[nbr];
   double ***backwardlo = new double**[nbr];
   double ***crosshi = new double**[nbr];
   double ***crossmi = new double**[nbr];
   double ***crosslo = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardhi[k] = new double*[nvrk];
      forwardmi[k] = new double*[nvrk];
      forwardlo[k] = new double*[nvrk];
      for(int i=0; i<nvrk; i++) 
      {
         forwardhi[k][i] = new double[deg+1];
         forwardmi[k][i] = new double[deg+1];
         forwardlo[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardhi[k] = new double*[nvrk-1];
         backwardmi[k] = new double*[nvrk-1];
         backwardlo[k] = new double*[nvrk-1];
         for(int i=0; i<nvrk-1; i++) 
         {
            backwardhi[k][i] = new double[deg+1];
            backwardmi[k][i] = new double[deg+1];
            backwardlo[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crosshi[k] = new double*[nvrk-2];
         crossmi[k] = new double*[nvrk-2];
         crosslo[k] = new double*[nvrk-2];
         for(int i=0; i<nvrk-2; i++)
         {
            crosshi[k][i] = new double[deg+1];
            crossmi[k][i] = new double[deg+1];
            crosslo[k][i] = new double[deg+1];
         }
      }
   }
   clock_t start = clock();
   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      if(verbose) cout << "executing convolution jobs at layer "
                       << k << " :" << endl;
      for(int i=0; i<cnvjobs.get_layer_count(k); i++)
      {
         ConvolutionJob job = cnvjobs.get_job(k,i);
         if(verbose) cout << "job " << i << " : " << job << endl;

         int monidx = job.get_monomial_index();

         CPU_dbl3_conv_job
            (deg,nvr[monidx],idx[monidx],
             cffhi[monidx],cffmi[monidx],cfflo[monidx],
             inputhi,inputmi,inputlo,
             forwardhi[monidx],forwardmi[monidx],forwardlo[monidx],
             backwardhi[monidx],backwardmi[monidx],backwardlo[monidx],
             crosshi[monidx],crossmi[monidx],crosslo[monidx],job,verbose);
      }
   }
   //CPU_dbl_poly_updates
   //   (dim,nbr,deg,nvr,idx,cst,cff,input,output,forward,backward,cross);
   CPU_dbl3_poly_addjobs
      (dim,nbr,deg,nvr,idx,csthi,cstmi,cstlo,cffhi,cffmi,cfflo,
       inputhi,inputmi,inputlo,outputhi,outputmi,outputlo,
       forwardhi,forwardmi,forwardlo,backwardhi,backwardmi,backwardlo,
       crosshi,crossmi,crosslo,addjobs,verbose);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(verbose)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
   }
   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k];

      for(int i=0; i<nvrk; i++)
      {
         free(forwardhi[k][i]);
         free(forwardmi[k][i]);
         free(forwardlo[k][i]);
      }
      if(nvrk > 1) for(int i=0; i<nvrk-1; i++)
                   {
                      free(backwardhi[k][i]);
                      free(backwardmi[k][i]);
                      free(backwardlo[k][i]);
                   }
      if(nvrk > 2) for(int i=0; i<nvrk-2; i++)
                   {
                      free(crosshi[k][i]);
                      free(crossmi[k][i]);
                      free(crosslo[k][i]);
                   }
   }
   free(forwardhi); free(backwardhi); free(crosshi);
   free(forwardmi); free(backwardmi); free(crossmi);
   free(forwardlo); free(backwardlo); free(crosslo);
}