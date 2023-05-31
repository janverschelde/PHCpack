/* The file dbl_polynomials_host.cpp defines functions specified
 * in dbl_polynomials_host.h. */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "dbl_convolutions_host.h"
#include "dbl_monomials_host.h"
#include "dbl_polynomials_host.h"

using namespace std;

void CPU_dbl_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cff, double **input, double **output,
   double **forward, double **backward, double **cross, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_dbl_product(deg,input[ix1],cff[i],forward[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            output[dim][j] += forward[0][j];
            output[ix1][j] += cff[i][j];
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_dbl_product(deg,cff[i],input[ix1],forward[0]);
         for(int j=0; j<=deg; j++) output[ix2][j] += forward[0][j];
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_dbl_product(deg,cff[i],input[ix2],backward[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) output[ix1][j] += backward[0][j];

         CPU_dbl_product(deg,forward[0],input[ix2],forward[1]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++) output[dim][j] += forward[1][j];
      }
      else if(nvr[i] > 2)
      {
         CPU_dbl_speel
            (nvr[i],deg,idx[i],cff[i],input,forward,backward,cross,verbose,i);

         ix1 = nvr[i]-1;
         for(int j=0; j<=deg; j++)     // update the value of the polynomial
            output[dim][j] += forward[ix1][j];

         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++) output[ix2][j] += forward[ix1][j];

         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++) output[ix2][j] += backward[ix1][j];

         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++) output[ix2][j] += cross[k-1][j];
         }
      }
   }
}

void CPU_cmplx_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffre, double **cffim, double **inputre, double **inputim,
   double **outputre, double **outputim,
   double **forwardre, double **forwardim,
   double **backwardre, double **backwardim,
   double **crossre, double **crossim, bool verbose )
{
   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];
         CPU_cmplx_product(deg,inputre[ix1],inputim[ix1],
                           cffre[i],cffim[i],forwardre[0],forwardim[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "input[" << ix1 << "] * cff to f[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            outputre[dim][j] += forwardre[0][j];
            outputim[dim][j] += forwardim[0][j];
            outputre[ix1][j] += cffre[i][j];
            outputim[ix1][j] += cffim[i][j];
         }
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         CPU_cmplx_product(deg,cffre[i],cffim[i],inputre[ix1],inputim[ix1],
                           forwardre[0],forwardim[0]);
         for(int j=0; j<=deg; j++) 
         {
            outputre[ix2][j] += forwardre[0][j];
            outputim[ix2][j] += forwardim[0][j];
         }
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to f[0]" << endl;

         CPU_cmplx_product(deg,cffre[i],cffim[i],inputre[ix2],inputim[ix2],
                           backwardre[0],backwardim[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++)
         {
            outputre[ix1][j] += backwardre[0][j];
            outputim[ix1][j] += backwardim[0][j];
         }
         CPU_cmplx_product(deg,forwardre[0],forwardim[0],inputre[ix2],
                           inputim[ix2],forwardre[1],forwardim[1]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "f[0] * "
                          << "input[" << ix2 << "] to f[1]" << endl;
         for(int j=0; j<=deg; j++)
         {
            outputre[dim][j] += forwardre[1][j];
            outputim[dim][j] += forwardim[1][j];
         }
      }
      else if(nvr[i] > 2)
      {
         CPU_cmplx_speel
            (nvr[i],deg,idx[i],cffre[i],cffim[i],inputre,inputim,
             forwardre,forwardim,backwardre,backwardim,crossre,crossim);

         ix1 = nvr[i]-1;
         for(int j=0; j<=deg; j++)     // update the value of the polynomial
         {
            outputre[dim][j] += forwardre[ix1][j];
            outputim[dim][j] += forwardim[ix1][j];
         }
         ix2 = idx[i][ix1];             // derivative with respect to x[n-1]
         ix1 = nvr[i]-2;

         for(int j=0; j<=deg; j++)
         {
            outputre[ix2][j] += forwardre[ix1][j];
            outputim[ix2][j] += forwardim[ix1][j];
         }
         ix2 = idx[i][0];                 // derivative with respect to x[0]
         ix1 = nvr[i]-3;

         for(int j=0; j<=deg; j++)
         {
            outputre[ix2][j] += backwardre[ix1][j];
            outputim[ix2][j] += backwardim[ix1][j];
         }
         ix1 = nvr[i]-1;                  // derivative with respect to x[k]
         for(int k=1; k<ix1; k++)
         { 
            ix2 = idx[i][k];
            for(int j=0; j<=deg; j++)
            {
               outputre[ix2][j] += crossre[k-1][j];
               outputim[ix2][j] += crossim[k-1][j];
            }
         }
      }
   }
}

void CPU_dbl_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double **forward = new double*[dim];
   double **backward = new double*[dim-1]; // in case dim = 2
   double **cross = new double*[dim-1];    // in case dim = 2

   for(int i=0; i<dim-1; i++)
   {
      forward[i] = new double[deg+1];
      backward[i] = new double[deg+1];
      cross[i] = new double[deg+1];
   }
   forward[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++) output[dim][i] = cst[i];
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) output[i][j] = 0.0;

   clock_t start = clock();
   CPU_dbl_poly_speel
      (dim,nbr,deg,nvr,idx,cff,input,output,forward,backward,cross,vrb);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(vrblvl > 0)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
      cout << scientific << setprecision(16);
   }
   for(int i=0; i<dim-1; i++)
   {
      free(forward[i]); free(backward[i]); free(cross[i]);
   }
   free(forward[dim-1]);
   free(forward); free(backward); free(cross);
}

void CPU_cmplx_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstre, double *cstim, double **cffre, double **cffim,
   double **inputre, double **inputim, double **outputre, double **outputim,
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double **forwardre = new double*[dim];
   double **forwardim = new double*[dim];
   double **backwardre = new double*[dim-1]; // in case dim = 2
   double **backwardim = new double*[dim-1];
   double **crossre = new double*[dim-1];    // in case dim = 2
   double **crossim = new double*[dim-1];

   for(int i=0; i<dim-1; i++)
   {
      forwardre[i] = new double[deg+1]; forwardim[i] = new double[deg+1];
      backwardre[i] = new double[deg+1]; backwardim[i] = new double[deg+1];
      crossre[i] = new double[deg+1]; crossim[i] = new double[deg+1];
   }
   forwardre[dim-1] = new double[deg+1];
   forwardim[dim-1] = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      outputre[dim][i] = cstre[i];
      outputim[dim][i] = cstim[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputre[i][j] = 0.0;
         outputim[i][j] = 0.0;
      }

   clock_t start = clock();
   CPU_cmplx_poly_speel
      (dim,nbr,deg,nvr,idx,cffre,cffim,inputre,inputim,outputre,outputim,
       forwardre,forwardim,backwardre,backwardim,crossre,crossim,vrb);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(vrblvl > 0)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
      cout << scientific << setprecision(16);
   }
   for(int i=0; i<dim-1; i++)
   {
      free(forwardre[i]); free(backwardre[i]); free(crossre[i]);
      free(forwardim[i]); free(backwardim[i]); free(crossim[i]);
   }
   free(forwardre[dim-1]); free(forwardim[dim-1]);
   free(forwardre); free(backwardre); free(crossre);
   free(forwardim); free(backwardim); free(crossim);
}

void CPU_dbl_conv_job
 ( int deg, int nvr, int *idx, double *cff, double **input,
   double **forward, double **backward, double **cross,
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
         CPU_dbl_product(deg,cff,input[inp2ix],forward[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl_product(deg,input[inp1ix],cff,forward[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_dbl_product
               (deg,input[inp1ix],forward[inp2ix],forward[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_dbl_product(deg,cross[inp1ix],input[inp2ix],forward[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_dbl_product(deg,input[inp1ix],cff,forward[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl_product
               (deg,forward[inp1ix],input[inp2ix],forward[outidx]);
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
            CPU_dbl_product(deg,cff,input[inp2ix],backward[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_dbl_product(deg,cff,backward[inp2ix],backward[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_dbl_product(deg,input[inp1ix],input[inp2ix],backward[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_dbl_product
               (deg,input[inp1ix],backward[inp2ix],backward[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_dbl_product(deg,backward[inp1ix],cff,backward[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_dbl_product
               (deg,backward[inp1ix],input[inp2ix],backward[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_dbl_product(deg,cff,input[inp2ix],cross[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl_product(deg,input[inp1ix],forward[inp2ix],cross[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_dbl_product(deg,forward[inp1ix],input[inp2ix],cross[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_dbl_product(deg,forward[inp1ix],backward[inp2ix],cross[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_dbl_product(deg,backward[inp1ix],forward[inp2ix],cross[outidx]);
      }
   }
}

void CPU_cmplx_conv_job
 ( int deg, int nvr, int *idx, double *cffre, double *cffim,
   double **inputre, double **inputim,
   double **forwardre, double **forwardim,
   double **backwardre, double **backwardim,
   double **crossre, double **crossim,
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
         CPU_cmplx_product(deg,cffre,cffim,inputre[inp2ix],inputim[inp2ix],
                           forwardre[outidx],forwardim[outidx]);
      }
      else if(inp1tp == 0)
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_cmplx_product(deg,inputre[inp1ix],inputim[inp1ix],cffre,
                              cffim,forwardre[outidx],forwardim[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * f[" << inp2ix << "]" << endl;
            CPU_cmplx_product(deg,inputre[inp1ix],inputim[inp1ix],
                              forwardre[inp2ix],forwardim[inp2ix],
                              forwardre[outidx],forwardim[outidx]);
         }
      }
      else if(inp1tp == 3)
      {
         if(verbose) cout << "c[" << inp1ix
                          << "] * input[" << inp2ix << "]" << endl;
         CPU_cmplx_product(deg,crossre[inp1ix],crossim[inp1ix],
                           inputre[inp2ix],inputim[inp2ix],
                           forwardre[outidx],forwardim[outidx]);
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "input[" << inp1ix << "] * cff" << endl;
            CPU_cmplx_product(deg,inputre[inp1ix],inputim[inp1ix],cffre,
                              cffim,forwardre[outidx],forwardim[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "f[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_cmplx_product(deg,forwardre[inp1ix],forwardim[inp1ix],
                              inputre[inp2ix],inputim[inp2ix],
                              forwardre[outidx],forwardim[outidx]);
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
            CPU_cmplx_product(deg,cffre,cffim,
                              inputre[inp2ix],inputim[inp2ix],
                              backwardre[outidx],backwardim[outidx]);
         }
         else
         {
            if(verbose) cout << "cff * b[" << inp2ix << "]" << endl;
            CPU_cmplx_product(deg,cffre,cffim,
                              backwardre[inp2ix],backwardim[inp2ix],
                              backwardre[outidx],backwardim[outidx]);
         }
      }
      else if(inp1tp == 0)
      {
         if(inp2tp == 0)
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * input[" << inp2ix << endl;
            CPU_cmplx_product(deg,inputre[inp1ix],inputim[inp1ix],
                              inputre[inp2ix],inputim[inp2ix],
                              backwardre[outidx],backwardim[outidx]);
         }
         else
         {
            if(verbose) cout << "input[" << inp1ix
                             << "] * b[" << inp2ix << "]" << endl;
            CPU_cmplx_product(deg,inputre[inp1ix],inputim[inp1ix],
                              backwardre[inp2ix],backwardim[inp2ix],
                              backwardre[outidx],backwardim[outidx]);
         }
      }
      else
      {
         if(inp2tp < 0)
         {
            if(verbose) cout << "b[" << inp1ix << "] * cff" << endl;
            CPU_cmplx_product(deg,backwardre[inp1ix],backwardim[inp1ix],
                              cffre,cffim,
                              backwardre[outidx],backwardim[outidx]);
         }
         else if(inp2tp == 0)
         {
            if(verbose) cout << "b[" << inp1ix
                             << "] * input[" << inp2ix << "]" << endl;
            CPU_cmplx_product(deg,backwardre[inp1ix],backwardim[inp1ix],
                              inputre[inp2ix],inputim[inp2ix],
                              backwardre[outidx],backwardim[outidx]);
         }
      }
   }
   else if(outptp == 3) // cross product either initializes or accumulates
   {
      if(verbose) cout << "-> computing c[" << outidx << "] = ";
      if(inp1tp < 0)
      {
         if(verbose) cout << "cff * input[" << inp2ix << "]" << endl;
         CPU_cmplx_product(deg,cffre,cffim,inputre[inp2ix],inputim[inp2ix],
                           crossre[outidx],crossim[outidx]);
      }
      if(inp1tp == 0)
      {
         if(verbose) cout << "input[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_cmplx_product(deg,inputre[inp1ix],inputim[inp1ix],
                           forwardre[inp2ix],forwardim[inp2ix],
                           crossre[outidx],crossim[outidx]);
      }
      else if(inp1tp == 1)
      {
        if(inp2tp == 0)
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * input[" << inp2ix << "]" << endl;
           CPU_cmplx_product(deg,forwardre[inp1ix],forwardim[inp1ix],
                             inputre[inp2ix],inputim[inp2ix],
                             crossre[outidx],crossim[outidx]);
        }
        else
        {
           if(verbose) cout << "f[" << inp1ix
                            << "] * b[" << inp2ix << "]" << endl;
           CPU_cmplx_product(deg,forwardre[inp1ix],forwardim[inp1ix],
                             backwardre[inp2ix],backwardim[inp2ix],
                             crossre[outidx],crossim[outidx]);
        }
      }
      else if(inp1tp == 2)
      {
         if(verbose) cout << "b[" << inp1ix
                          << "] * f[" << inp2ix << "]" << endl;
         CPU_cmplx_product(deg,backwardre[inp1ix],backwardim[inp1ix],
                           forwardre[inp2ix],forwardim[inp2ix],
                           crossre[outidx],crossim[outidx]);
      }
   }
}

void CPU_dbl_add_job
 ( int deg, double *cst, double **cff,
   double ***forward, double ***backward, double ***cross,
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
               forward[updmon][updidx][i] += cst[i];
         else
            for(int i=0; i<=deg; i++)
               forward[updmon][updidx][i] += cff[incidx][i];
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            forward[updmon][updidx][i] += forward[incmon][incidx][i];
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            forward[updmon][updidx][i] += backward[incmon][incidx][i];
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            forward[updmon][updidx][i] += cross[incmon][incidx][i];
      }
   }
   else if(adtype == 2)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            backward[updmon][updidx][i] += cff[incidx][i];
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            backward[updmon][updidx][i] += forward[incmon][incidx][i];
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            backward[updmon][updidx][i] += backward[incmon][incidx][i];
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            backward[updmon][updidx][i] += cross[incmon][incidx][i];
      }
   }
   else if(adtype == 3)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
            cross[updmon][updidx][i] += cff[incidx][i];
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
            cross[updmon][updidx][i] += forward[incmon][incidx][i];
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
            cross[updmon][updidx][i] += backward[incmon][incidx][i];
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
            cross[updmon][updidx][i] += cross[incmon][incidx][i];
      }
   }
}

void CPU_cmplx_add_job
 ( int deg, double *cstre, double *cstim,
   double **cffre, double **cffim,
   double ***forwardre, double ***forwardim,
   double ***backwardre, double ***backwardim,
   double ***crossre, double ***crossim,
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
            {
               forwardre[updmon][updidx][i] += cstre[i];
               forwardim[updmon][updidx][i] += cstim[i];
            }
         else
            for(int i=0; i<=deg; i++)
            {
               forwardre[updmon][updidx][i] += cffre[incidx][i];
               forwardim[updmon][updidx][i] += cffim[incidx][i];
            }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
         {
            forwardre[updmon][updidx][i] += forwardre[incmon][incidx][i];
            forwardim[updmon][updidx][i] += forwardim[incmon][incidx][i];
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
         {
            forwardre[updmon][updidx][i] += backwardre[incmon][incidx][i];
            forwardim[updmon][updidx][i] += backwardim[incmon][incidx][i];
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
         {
            forwardre[updmon][updidx][i] += crossre[incmon][incidx][i];
            forwardim[updmon][updidx][i] += crossim[incmon][incidx][i];
         }
      }
   }
   else if(adtype == 2)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
         {
            backwardre[updmon][updidx][i] += cffre[incidx][i];
            backwardim[updmon][updidx][i] += cffim[incidx][i];
         }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
         {
            backwardre[updmon][updidx][i] += forwardre[incmon][incidx][i];
            backwardim[updmon][updidx][i] += forwardim[incmon][incidx][i];
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
         {
            backwardre[updmon][updidx][i] += backwardre[incmon][incidx][i];
            backwardim[updmon][updidx][i] += backwardim[incmon][incidx][i];
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
         {
            backwardre[updmon][updidx][i] += crossre[incmon][incidx][i];
            backwardim[updmon][updidx][i] += crossim[incmon][incidx][i];
         }
      }
   }
   else if(adtype == 3)
   {
      if(incmon < 0)
      {
         for(int i=0; i<=deg; i++)
         {
            crossre[updmon][updidx][i] += cffre[incidx][i];
            crossim[updmon][updidx][i] += cffim[incidx][i];
         }
      }
      else if(intype == 1)
      {
         for(int i=0; i<=deg; i++)
         {
            crossre[updmon][updidx][i] += forwardre[incmon][incidx][i];
            crossim[updmon][updidx][i] += forwardim[incmon][incidx][i];
         }
      }
      else if(intype == 2)
      {
         for(int i=0; i<=deg; i++)
         {
            crossre[updmon][updidx][i] += backwardre[incmon][incidx][i];
            crossim[updmon][updidx][i] += backwardim[incmon][incidx][i];
         }
      }
      else if(intype == 3)
      {
         for(int i=0; i<=deg; i++)
         {
            crossre[updmon][updidx][i] += crossre[incmon][incidx][i];
            crossim[updmon][updidx][i] += crossim[incmon][incidx][i];
         }
      }
   }
}

void CPU_dbl_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   double ***forward, double ***backward, double ***cross )
{
   for(int i=0; i<=deg; i++) output[dim][i] = cst[i];

   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) output[i][j] = 0.0;

   for(int k=0; k<nbr; k++)
   {
      int ix0 = idx[k][0];   // first variable in monomial k
      int ix1 = nvr[k]-1;    // last forward has the value
      int ix2 = nvr[k]-2;    // next to last forward has last derivative
                             // last backward has the first derivative
      int ixn = idx[k][ix1]; // index of the last variable in monomial k

      for(int i=0; i<=deg; i++) // value is last forward location
         output[dim][i] = output[dim][i] + forward[k][ix1][i];

      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
            output[ix0][i] = output[ix0][i] + cff[k][i]; 
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            output[ixn][i] = output[ixn][i] + forward[k][ix2][i];
            output[ix0][i] = output[ix0][i] + backward[k][ix2][i];
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
                  output[ix0][i] = output[ix0][i] + cross[k][j-1][i];
            }
         }
      }
   }
}

void CPU_cmplx_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstre, double *cstim,
   double **cffre, double **cffim,
   double **inputre, double **inputim,
   double **outputre, double **outputim,
   double ***forwardre, double ***forwardim,
   double ***backwardre, double ***backwardim,
   double ***crossre, double ***crossim )
{
   for(int i=0; i<=deg; i++)
   {
      outputre[dim][i] = cstre[i];
      outputim[dim][i] = cstim[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputre[i][j] = 0.0;
         outputim[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      int ix0 = idx[k][0];   // first variable in monomial k
      int ix1 = nvr[k]-1;    // last forward has the value
      int ix2 = nvr[k]-2;    // next to last forward has last derivative
                             // last backward has the first derivative
      int ixn = idx[k][ix1]; // index of the last variable in monomial k

      for(int i=0; i<=deg; i++) // value is last forward location
      {
         outputre[dim][i] = outputre[dim][i] + forwardre[k][ix1][i];
         outputim[dim][i] = outputim[dim][i] + forwardim[k][ix1][i];
      }
      if(ix1 == 0)           // monomial has only one variable
      {
         for(int i=0; i<=deg; i++)
         {
            outputre[ix0][i] = outputre[ix0][i] + cffre[k][i]; 
            outputim[ix0][i] = outputim[ix0][i] + cffim[k][i]; 
         }
      }
      else if(ix2 >= 0)      // update first and last derivative
      {
         for(int i=0; i<=deg; i++)
         {
            outputre[ixn][i] = outputre[ixn][i] + forwardre[k][ix2][i];
            outputim[ixn][i] = outputim[ixn][i] + forwardim[k][ix2][i];
            outputre[ix0][i] = outputre[ix0][i] + backwardre[k][ix2][i];
            outputim[ix0][i] = outputim[ix0][i] + backwardim[k][ix2][i];
         }
         if(ix2 > 0)         // update all other derivatives
         {
            for(int j=1; j<ix1; j++) // j-th variable in monomial k
            {
               ix0 = idx[k][j];
               for(int i=0; i<=deg; i++)
               {
                  outputre[ix0][i] = outputre[ix0][i] + crossre[k][j-1][i];
                  outputim[ix0][i] = outputim[ix0][i] + crossim[k][j-1][i];
               }
            }
         }
      }
   }
}

void CPU_dbl_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   double ***forward, double ***backward, double ***cross,
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

         CPU_dbl_add_job(deg,cst,cff,forward,backward,cross,job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
      output[dim][i] = forward[lastmon][lastidx][i];

   int cnt = jobs.get_differential_count(0);

   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      const int difidx = jobs.get_differential_index(0,0);

      if(verbose)
         cout << "Differential index for variable 0 : " << difidx << endl;

      if(difidx < 0)
      {
         for(int i=0; i<=deg; i++) output[0][i] = 0.0;
      }
      else
      {
         if(verbose)
            cout << "updating derivative 0 with coefficient "
                 << difidx << endl;

         for(int i=0; i<=deg; i++) output[0][i] = cff[difidx][i];
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

      for(int i=0; i<=deg; i++) output[0][i] = backward[ix0][ix2][i];
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);

      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         const int difidx = jobs.get_differential_index(k,0);

         if(verbose)
            cout << "Differential index for variable " << k 
                 << " : " << difidx << endl;

         if(difidx < 0)
         {
            for(int i=0; i<=deg; i++) output[k][i] = 0.0;
         }
         else
         {
            if(verbose)
               cout << "updating derivative " << k 
                    << " with coefficient " << difidx << endl;

            for(int i=0; i<=deg; i++) output[k][i] = cff[difidx][i];
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

            for(int i=0; i<=deg; i++) output[k][i] = backward[ix0][ix2][i];
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0] - 2;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : f[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) output[k][i] = forward[ix0][ix2][i];
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : c[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++) output[k][i] = cross[ix0][ix2][i];
         }
      }
   }
}

void CPU_cmplx_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstre, double *cstim,
   double **cffre, double **cffim,
   double **inputre, double **inputim,
   double **outputre, double **outputim,
   double ***forwardre, double ***forwardim,
   double ***backwardre, double ***backwardim,
   double ***crossre, double ***crossim,
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

         CPU_cmplx_add_job(deg,cstre,cstim,cffre,cffim,
                           forwardre,forwardim,backwardre,backwardim,
                           crossre,crossim,job,verbose);
      }
   }
   int lastmon = nbr-1;
   int lastidx = nvr[lastmon]-1;
   for(int i=0; i<=deg; i++) // value is last forward location
   {
      outputre[dim][i] = forwardre[lastmon][lastidx][i];
      outputim[dim][i] = forwardim[lastmon][lastidx][i];
   }
   int cnt = jobs.get_differential_count(0);

   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      const int difidx = jobs.get_differential_index(0,0);

      if(verbose)
         cout << "Differential index for variable 0 : " << difidx << endl;

      if(difidx < 0)
      {
         for(int i=0; i<=deg; i++)
         {
            outputre[0][i] = 0.0;
            outputim[0][i] = 0.0;
         }
      }
      else
      {
         if(verbose)
            cout << "updating derivative 0 with coefficient "
                 << difidx << endl;

         for(int i=0; i<=deg; i++)
         {
            outputre[0][i] = cffre[difidx][i];
            outputim[0][i] = cffim[difidx][i];
         }
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

      for(int i=0; i<=deg; i++)
      {
         outputre[0][i] = backwardre[ix0][ix2][i];
         outputim[0][i] = backwardim[ix0][ix2][i];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);

      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         const int difidx = jobs.get_differential_index(k,0);

         if(verbose)
            cout << "Differential index for variable " << k 
                 << " : " << difidx << endl;

         if(difidx < 0)
         {
            for(int i=0; i<=deg; i++)
            {
               outputre[k][i] = 0.0;
               outputim[k][i] = 0.0;
            }
         }
         else
         {
            if(verbose)
               cout << "updating derivative " << k 
                    << " with coefficient " << difidx << endl;

            for(int i=0; i<=deg; i++)
            {
               outputre[k][i] = cffre[difidx][i];
               outputim[k][i] = cffim[difidx][i];
            }
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

            for(int i=0; i<=deg; i++)
            {
               outputre[k][i] = backwardre[ix0][ix2][i];
               outputim[k][i] = backwardim[ix0][ix2][i];
            }
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0] - 2;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : f[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++)
            {
               outputre[k][i] = forwardre[ix0][ix2][i];
               outputim[k][i] = forwardim[ix0][ix2][i];
            }
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;

            if(verbose)
               cout << "Updating derivative " << k 
                    << ", ix0 = " << ix0 << ", ix2 = " << ix2
                    << " : c[" << ix0 << "," << ix2 << "]" << endl;

            for(int i=0; i<=deg; i++)
            {
               outputre[k][i] = crossre[ix0][ix2][i];
               outputim[k][i] = crossim[ix0][ix2][i];
            }
         }
      }
   }
}

void CPU_dbl_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);

   double ***forward = new double**[nbr];
   double ***backward = new double**[nbr];
   double ***cross = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forward[k] = new double*[nvrk];
      for(int i=0; i<nvrk; i++) forward[k][i] = new double[deg+1];

      if(nvrk > 1)
      {
         backward[k] = new double*[nvrk-1];
         for(int i=0; i<nvrk-1; i++) backward[k][i] = new double[deg+1];
      }
      if(nvrk > 2)
      {
         cross[k] = new double*[nvrk-2];
         for(int i=0; i<nvrk-2; i++) cross[k][i] = new double[deg+1];
      }
   }
   clock_t start = clock();
   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      if(vrb) cout << "executing convolution jobs at layer "
                   << k << " :" << endl;
      for(int i=0; i<cnvjobs.get_layer_count(k); i++)
      {
         ConvolutionJob job = cnvjobs.get_job(k,i);
         if(vrb) cout << "job " << i << " : " << job << endl;

         int monidx = job.get_monomial_index();

         CPU_dbl_conv_job
            (deg,nvr[monidx],idx[monidx],cff[monidx],input,
             forward[monidx],backward[monidx],cross[monidx],job,vrb);
      }
   }
   //CPU_dbl_poly_updates
   //   (dim,nbr,deg,nvr,idx,cst,cff,input,output,forward,backward,cross);
   CPU_dbl_poly_addjobs
      (dim,nbr,deg,nvr,idx,cst,cff,input,output,forward,backward,cross,
       addjobs,vrb);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(vrblvl > 0)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
      cout << scientific << setprecision(16);
   }
   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k];

      for(int i=0; i<nvrk; i++) free(forward[k][i]);
      if(nvrk > 1) for(int i=0; i<nvrk-1; i++) free(backward[k][i]);
      if(nvrk > 2) for(int i=0; i<nvrk-2; i++) free(cross[k][i]);
   }
   free(forward); free(backward); free(cross);
}

void CPU_cmplx_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstre, double *cstim,
   double **cffre, double **cffim,
   double **inputre, double **inputim,
   double **outputre, double **outputim,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl )
{
   const bool vrb = (vrblvl > 1);
   double ***forwardre = new double**[nbr];
   double ***forwardim = new double**[nbr];
   double ***backwardre = new double**[nbr];
   double ***backwardim = new double**[nbr];
   double ***crossre = new double**[nbr];
   double ***crossim = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k]; // number of variables in monomial k

      forwardre[k] = new double*[nvrk];
      forwardim[k] = new double*[nvrk];
      for(int i=0; i<nvrk; i++)
      {
         forwardre[k][i] = new double[deg+1];
         forwardim[k][i] = new double[deg+1];
      }
      if(nvrk > 1)
      {
         backwardre[k] = new double*[nvrk-1];
         backwardim[k] = new double*[nvrk-1];
         for(int i=0; i<nvrk-1; i++)
         {
            backwardre[k][i] = new double[deg+1];
            backwardim[k][i] = new double[deg+1];
         }
      }
      if(nvrk > 2)
      {
         crossre[k] = new double*[nvrk-2];
         crossim[k] = new double*[nvrk-2];
         for(int i=0; i<nvrk-2; i++) 
         {
            crossre[k][i] = new double[deg+1];
            crossim[k][i] = new double[deg+1];
         }
      }
   }
   clock_t start = clock();
   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      if(vrb) cout << "executing convolution jobs at layer "
                   << k << " :" << endl;
      for(int i=0; i<cnvjobs.get_layer_count(k); i++)
      {
         ConvolutionJob job = cnvjobs.get_job(k,i);
         if(vrb) cout << "job " << i << " : " << job << endl;

         int monidx = job.get_monomial_index();

         CPU_cmplx_conv_job
            (deg,nvr[monidx],idx[monidx],cffre[monidx],cffim[monidx],
             inputre,inputim,forwardre[monidx],forwardim[monidx],
             backwardre[monidx],backwardim[monidx],crossre[monidx],
             crossim[monidx],job,vrb);
      }
   }
   //CPU_dbl_poly_updates
   //   (dim,nbr,deg,nvr,idx,cst,cff,input,output,forward,backward,cross);
   CPU_cmplx_poly_addjobs
      (dim,nbr,deg,nvr,idx,cstre,cstim,cffre,cffim,inputre,outputim,
       outputre,outputim,forwardre,forwardim,backwardre,backwardim,
       crossre,crossim,addjobs,vrb);
   clock_t end = clock();
   *elapsedsec = double(end - start)/CLOCKS_PER_SEC;

   if(vrblvl > 0)
   {
      cout << fixed << setprecision(3);
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << *elapsedsec << " seconds." << endl;
      cout << scientific << setprecision(16);
   }
   for(int k=0; k<nbr; k++)
   {
      int nvrk = nvr[k];

      for(int i=0; i<nvrk; i++)
      {
         free(forwardre[k][i]); free(forwardim[k][i]);
      }
      if(nvrk > 1) for(int i=0; i<nvrk-1; i++)
                   {
                      free(backwardre[k][i]); free(backwardim[k][i]);
                   }
      if(nvrk > 2) for(int i=0; i<nvrk-2; i++)
                   {
                      free(crossre[k][i]); free(crossim[k][i]);
                   }
   }
   free(forwardre); free(backwardre); free(crossre);
   free(forwardim); free(backwardim); free(crossim);
}
