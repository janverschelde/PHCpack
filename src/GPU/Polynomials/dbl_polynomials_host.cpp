/* The file dbl_polynomials_host.cpp defines functions specified
 * in dbl_polynomials_host.h. */

#include <cstdlib>
#include <iostream>
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

         CPU_dbl_product(deg,cff[i],input[ix1],cross[0]);
         for(int j=0; j<=deg; j++) output[ix2][j] += cross[0][j];
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix1 << "] to c[0]" << endl;

         CPU_dbl_product(deg,cff[i],input[ix2],backward[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "cff * "
                          << "input[" << ix2 << "] to b[0]" << endl;
         for(int j=0; j<=deg; j++) output[ix1][j] += backward[0][j];

         CPU_dbl_product(deg,cross[0],input[ix2],forward[0]);
         if(verbose) cout << "monomial " << i << " : ";
         if(verbose) cout << "c[0] * "
                          << "input[" << ix2 << "] to f[0]" << endl;
         for(int j=0; j<=deg; j++) output[dim][j] += forward[0][j];
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

void CPU_dbl_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output, bool verbose )
{
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

   CPU_dbl_poly_speel
      (dim,nbr,deg,nvr,idx,cff,input,output,forward,backward,cross,verbose);

   for(int i=0; i<dim-1; i++)
   {
      free(forward[i]); free(backward[i]); free(cross[i]);
   }
   free(forward[dim-1]);
   free(forward); free(backward); free(cross);
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

void CPU_dbl_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   ConvolutionJobs jobs, bool verbose )
{
   double ***forward = new double**[nbr];
   double ***backward = new double**[nbr];
   double ***cross = new double**[nbr];

   for(int k=0; k<nbr; k++)
   {
      forward[k] = new double*[dim];
      backward[k] = new double*[dim-1];
      cross[k] = new double*[dim-1];
      for(int i=0; i<dim-1; i++)
      {
         forward[k][i] = new double[deg+1];
         backward[k][i] = new double[deg+1];
         cross[k][i] = new double[deg+1];
      }
      forward[k][dim-1] = new double[deg+1];
   }
   for(int i=0; i<=deg; i++) output[dim][i] = cst[i];
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) output[i][j] = 0.0;

   for(int k=0; k<jobs.get_depth(); k++)
   {
      if(verbose) cout << "running jobs at layer " << k << " :" << endl;
      for(int i=0; i<jobs.get_layer_count(k); i++)
      {
         ConvolutionJob job = jobs.get_job(k,i);
         if(verbose) cout << "job " << i << " : " << job << endl;

         int monidx = job.get_monomial_index();

         CPU_dbl_conv_job
           (deg,nvr[monidx],idx[monidx],cff[monidx],input,
            forward[monidx],backward[monidx],cross[monidx],job,verbose);
      }
   }

   for(int k=0; k<nbr; k++)
   {
      for(int i=0; i<=deg; i++)
         output[dim][i] = output[dim][i] + forward[k][nvr[k]-1][i];
   }

/*
   for(int k=0; k<nbr; k++)
   {
      for(int i=0; i<dim-1; i++)
      {
         free(forward[k][i]); free(backward[k][i]); free(cross[k][i]);
      }
      free(forward[k][dim-1]);
   }
   free(forward); free(backward); free(cross);
 */
}
