// The file job_coordinates.cpp defines the functions specified
// in job_coordinates.h.

#include <iostream>
#include "job_coordinates.h"

using namespace std;

int coefficient_count ( int dim, int nbr, int deg, int *nvr )
{
   int count = 1 + nbr + dim;

   for(int i=0; i<nbr; i++)
   {
      count = count + nvr[i];
      if(nvr[i] == 2)
         count = count + 1;
      else if(nvr[i] > 2)
         count = count + 2*(nvr[i]-2);
   }
   count = count*(deg+1);

   return count;
}

int complex_coefficient_count ( int dim, int nbr, int deg, int *nvr )
{
   int count = 1 + nbr + dim;
   // constant, coefficients of nbr monomials, dim input variables

   for(int i=0; i<nbr; i++)
   {
      count = count + nvr[i]; // forward products
      if(nvr[i] == 2)
         count = count + 1;   // one backward product
      else if(nvr[i] > 2)
         count = count + (nvr[i]-1) + (nvr[i]-2);
         // add backward and cross products
   }
   count = count*(deg+1);

   return count;
}

void coefficient_indices
 ( int dim, int nbr, int deg, int *nvr,
   int *fsums, int *bsums, int *csums,
   int *fstart, int *bstart, int *cstart )
{
   fsums[0] = nvr[0]; bsums[0] = 0; csums[0] = 0;

   if(nvr[0] == 2)
   {
      bsums[0] = 1;
   }
   else if(nvr[0] > 2)
   {
      bsums[0] = nvr[0] - 2;
      csums[0] = nvr[0] - 2;
   }
   for(int i=1; i<nbr; i++)
   {
      fsums[i] = fsums[i-1] + nvr[i];
      if(nvr[i] < 2)
      {
         bsums[i] = bsums[i-1];
         csums[i] = csums[i-1];
      }
      else if(nvr[i] == 2)
      {
         bsums[i] = bsums[i-1] + 1;
         csums[i] = csums[i-1];
      }
      else // nvr[i] > 2
      {
         bsums[i] = bsums[i-1] + nvr[i] - 2;
         csums[i] = csums[i-1] + nvr[i] - 2;
      }
   }
   fstart[0] = (1+nbr+dim)*(deg+1);
   for(int i=1; i<nbr; i++) fstart[i] = fstart[0] + fsums[i-1]*(deg+1);

   bstart[0] = fstart[0] + fsums[nbr-1]*(deg+1);
   for(int i=1; i<nbr; i++) bstart[i] = bstart[0] + bsums[i-1]*(deg+1);

   cstart[0] = bstart[0] + bsums[nbr-1]*(deg+1);
   for(int i=1; i<nbr; i++) cstart[i] = cstart[0] + csums[i-1]*(deg+1);
}

void complex_coefficient_indices
 ( int dim, int nbr, int deg, int *nvr,
   int *fsums, int *bsums, int *csums,
   int *fstart, int *bstart, int *cstart )
{
   fsums[0] = nvr[0]; bsums[0] = 0; csums[0] = 0;

   if(nvr[0] == 2)
   {
      bsums[0] = 1;                    // need one backward product
   }
   else if(nvr[0] > 2)
   {
      bsums[0] = nvr[0] - 1;
      csums[0] = nvr[0] - 2;
   }
   for(int i=1; i<nbr; i++)
   {
      fsums[i] = fsums[i-1] + nvr[i];  // as many forward as #variables
      if(nvr[i] < 2)
      {
         bsums[i] = bsums[i-1];        // no backward product
         csums[i] = csums[i-1];        // no cross product
      }
      else if(nvr[i] == 2)
      {
         bsums[i] = bsums[i-1] + 1;   // need one backward product
         csums[i] = csums[i-1];       // no cross product needed
      }
      else // nvr[i] > 2
      {
         bsums[i] = bsums[i-1] + nvr[i] - 1;
         csums[i] = csums[i-1] + nvr[i] - 2;
      }
   }
   fstart[0] = (1+nbr+dim)*(deg+1);
   for(int i=1; i<nbr; i++) fstart[i] = fstart[0] + fsums[i-1]*(deg+1);

   bstart[0] = fstart[0] + fsums[nbr-1]*(deg+1);
   for(int i=1; i<nbr; i++) bstart[i] = bstart[0] + bsums[i-1]*(deg+1);

   cstart[0] = bstart[0] + bsums[nbr-1]*(deg+1);
   for(int i=1; i<nbr; i++) cstart[i] = cstart[0] + csums[i-1]*(deg+1);
}

void write_coefficient_indices
 ( int cnt, int nbr, int *fsums, int *fstart, int *bsums, int *bstart,
   int *csums, int *cstart )
{
   cout << "The output coefficient count : " << cnt << endl;
   cout << "fsums :";
   for(int i=0; i<nbr; i++) cout << " " << fsums[i]; cout << endl;
   cout << "fstart :";
   for(int i=0; i<nbr; i++) cout << " " << fstart[i]; cout << endl;
   cout << "bsums :";
   for(int i=0; i<nbr; i++) cout << " " << bsums[i]; cout << endl;
   cout << "bstart :";
   for(int i=0; i<nbr; i++) cout << " " << bstart[i]; cout << endl;
   cout << "csums :";
   for(int i=0; i<nbr; i++) cout << " " << csums[i]; cout << endl;
   cout << "cstart :";
   for(int i=0; i<nbr; i++) cout << " " << cstart[i]; cout << endl;
}

void convjob_indices
 ( ConvolutionJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int monidx = job.get_monomial_index();
   const int jobinp1tp = job.get_first_type();
   const int jobinp1ix = job.get_first_input();
   const int jobinp2tp = job.get_second_type();
   const int jobinp2ix = job.get_second_input();
   const int joboutptp = job.get_output_type();
   const int joboutidx = job.get_output_index();
   const int deg1 = deg+1;

   if(verbose)
   {
      cout << "  in1 type : " << jobinp1tp << ", idx : " << jobinp1ix;
      cout << "  in2 type : " << jobinp2tp << ", idx : " << jobinp2ix;
      cout << "  out type : " << joboutptp << ", idx : " << joboutidx;
      cout << endl;
   }
   if(jobinp1tp < 0)           // first input is coefficient
      *inp1ix = (1 + monidx)*deg1;
   else if(jobinp1tp == 0)     // first input is input series
      *inp1ix = (1 + nbr + jobinp1ix)*deg1;
   else if(jobinp1tp == 1)     // first input is forward product
      *inp1ix = fstart[monidx] + jobinp1ix*deg1;
   else if(jobinp1tp == 2)     // first input is backward product
      *inp1ix = bstart[monidx] + jobinp1ix*deg1;
   else if(jobinp1tp == 3)     // first input is cross product
      *inp1ix = cstart[monidx] + jobinp1ix*deg1;

   if(jobinp2tp < 0)           // second input is coefficient
      *inp2ix = (1 + monidx)*deg1;
   else if(jobinp2tp == 0)     // second input is input series
      *inp2ix = (1 + nbr + jobinp2ix)*deg1;
   else if(jobinp2tp == 1)     // second input is forward product
      *inp2ix = fstart[monidx] + jobinp2ix*deg1;
   else if(jobinp2tp == 2)     // second input is backward product
      *inp2ix = bstart[monidx] + jobinp2ix*deg1;
   else if(jobinp2tp == 3)     // second input is cross product
      // *outidx = cstart[monidx] + jobinp2ix*deg1;
      *inp2ix = cstart[monidx] + jobinp2ix*deg1;

   if(joboutptp == 1)          // output is forward product
      *outidx = fstart[monidx] + joboutidx*deg1;
   else if(joboutptp == 2)     // output is backward product
   {
      if(joboutidx < nvr[monidx]-2) // last backward product is special ...
         *outidx = bstart[monidx] + joboutidx*deg1;
      else
      {
         if(nvr[monidx] == 2)
            *outidx = bstart[monidx];
         else
            *outidx = bstart[monidx] + (joboutidx-1)*deg1;

      }
   }
   else if(joboutptp == 3)    // output is cross product
      *outidx = cstart[monidx] + joboutidx*deg1;

   if(verbose)
   {
      cout << "-> inp1ix : " << *inp1ix
           << ", inp2ix : " << *inp2ix
           << ", outidx : " << *outidx << endl;
   }
}

void complex_convjob_indices
 ( ComplexConvolutionJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr, int totcff, int offset,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int monidx = job.get_monomial_index();
   const int jobinp1tp = job.get_first_type();
   const int jobinp1ix = job.get_first_input();
   const int jobinp2tp = job.get_second_type();
   const int jobinp2ix = job.get_second_input();
   const int joboutptp = job.get_output_type();
   const int joboutidx = job.get_output_index();
   const int deg1 = deg+1;
   const int totcffoffset = totcff + offset; // start of imag data array

   if(verbose)
   {
      cout << "  in1 type : " << jobinp1tp << ", idx : " << jobinp1ix;
      cout << "  in2 type : " << jobinp2tp << ", idx : " << jobinp2ix;
      cout << "  out type : " << joboutptp << ", idx : " << joboutidx;
      cout << endl;
   }

   if(jobinp1tp == -2)       // first input is imag part of coefficient
      *inp1ix = (1 + monidx)*deg1 + totcffoffset;
   else if(jobinp1tp == -1)  // first input is real part of coefficient
      *inp1ix = (1 + monidx)*deg1;
   else if(jobinp1tp == 0)   // first input is real part of input series
      *inp1ix = (1 + nbr + jobinp1ix)*deg1;
   else if(jobinp1tp == 4)   // first input is imag part of input series
      *inp1ix = (1 + nbr + jobinp1ix)*deg1 + totcffoffset;
   else if(jobinp1tp == 1)   // first input is real part of forward
      *inp1ix = fstart[monidx] + jobinp1ix*deg1;
   else if(jobinp1tp == 5)   // first input is real part of forward
      *inp1ix = fstart[monidx] + jobinp1ix*deg1 + totcffoffset;
   else if(jobinp1tp == 2)   // first input is real part of backward
      *inp1ix = bstart[monidx] + jobinp1ix*deg1;
   else if(jobinp1tp == 6)   // first input is imag part of backward
      *inp1ix = bstart[monidx] + jobinp1ix*deg1 + totcffoffset;
   else if(jobinp1tp == 3)   // first input is real part of cross
      *inp1ix = cstart[monidx] + jobinp1ix*deg1;
   else if(jobinp1tp == 7)   // first input is imag part of cross
      *inp1ix = cstart[monidx] + jobinp1ix*deg1 + totcffoffset;
   else
      cout << "Invalid job first input type!" << endl;

   if(jobinp2tp == -2)       // second input is imag part of coefficient
      *inp2ix = (1 + monidx)*deg1 + totcffoffset;
   else if(jobinp2tp == -1)  // second input is real part of coefficient
      *inp2ix = (1 + monidx)*deg1;
   else if(jobinp2tp == 0)   // second input is real part of input series
      *inp2ix = (1 + nbr + jobinp2ix)*deg1;
   else if(jobinp2tp == 4)   // second input is imag part of input series
      *inp2ix = (1 + nbr + jobinp2ix)*deg1 + totcffoffset;
   else if(jobinp2tp == 1)   // second input is real part of forward
      *inp2ix = fstart[monidx] + jobinp2ix*deg1;
   else if(jobinp2tp == 5)   // second input is imag part of forward
      *inp2ix = fstart[monidx] + jobinp2ix*deg1 + totcffoffset;
   else if(jobinp2tp == 2)   // second input is real part of backward
      *inp2ix = bstart[monidx] + jobinp2ix*deg1;
   else if(jobinp2tp == 6)   // second input is imag part of backward
      *inp2ix = bstart[monidx] + jobinp2ix*deg1 + totcffoffset;
   else if(jobinp2tp == 3)   // second input is real part of cross
      // *outidx = cstart[monidx] + jobinp2ix*deg1;
      *inp2ix = cstart[monidx] + jobinp2ix*deg1;
   else if(jobinp2tp == 7)   // second input is imag part of cross
      // *outidx = cstart[monidx] + jobinp2ix*deg1 + totcffoffset;
      *inp2ix = cstart[monidx] + jobinp2ix*deg1 + totcffoffset;
   else
      cout << "Invalid job second input type!" << endl;

   if(joboutptp == 1)          // first operand of real forward product
      *outidx = fstart[monidx] + joboutidx*deg1;
   else if(joboutptp == 4)     // second operand of real forward product
      *outidx = fstart[monidx] + joboutidx*deg1 + offset;
   else if(joboutptp == 7)     // first operand of imag forward product
      *outidx = fstart[monidx] + joboutidx*deg1 + totcffoffset;
   else if(joboutptp == 10)    // second operand of imag forward product
      *outidx = fstart[monidx] + joboutidx*deg1 + totcffoffset + offset;
   else if(joboutptp == 2)     // first operand of real backward
   {
      if(nvr[monidx] == 2)
         *outidx = bstart[monidx]; // we have one backward product
      else
         *outidx = bstart[monidx] + joboutidx*deg1;
   }
   else if(joboutptp == 5)     // second operand of real backward
   {
      if(nvr[monidx] == 2)
         *outidx = bstart[monidx] + offset;
      else
         *outidx = bstart[monidx] + joboutidx*deg1 + offset;
   }
   else if(joboutptp == 8)     // first operand of imag backward
   {
      if(nvr[monidx] == 2)
         *outidx = bstart[monidx] + totcffoffset;
      else
         *outidx = bstart[monidx] + joboutidx*deg1 + totcffoffset;
   }
   else if(joboutptp == 11)    // second operand of imag backward
   {
      if(nvr[monidx] == 2)
         *outidx = bstart[monidx] + totcffoffset + offset;
      else
         *outidx = bstart[monidx] + joboutidx*deg1 + totcffoffset + offset;
   }
   else if(joboutptp == 3)    // first operand of real cross product
      *outidx = cstart[monidx] + joboutidx*deg1;
   else if(joboutptp == 6)    // second operand of real cross product
      *outidx = cstart[monidx] + joboutidx*deg1 + offset;
   else if(joboutptp == 9)    // first operand of imag cross product
      *outidx = cstart[monidx] + joboutidx*deg1 + totcffoffset;
   else if(joboutptp == 12)   // second operand of imag cross product
      *outidx = cstart[monidx] + joboutidx*deg1 + totcffoffset + offset;
   else
      cout << "Invalid job output type!" << endl;

   if(verbose)
   {
      cout << "-> inp1ix : " << *inp1ix
           << ", inp2ix : " << *inp2ix
           << ", outidx : " << *outidx << endl;
   }
}

void complex_incjob_indices
 ( ComplexIncrementJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr, int totcff, int offset,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int monidx = job.get_monomial_index();
   const int jobinkind = job.get_increment_kind();
   const int jobincidx = job.get_increment_index();
   const int jobisreal = job.get_isreal();
   const int deg1 = deg+1;
   const int totcffoffset = totcff + offset; // start of imag data array

   if(verbose)
   {
      cout << "  kind : " << jobinkind << "  index : " << jobincidx 
           << "  real or not : " << jobisreal << endl;
   }
   if(jobinkind == 1)   // incrementing forward product
   {
      if(jobisreal)
      {
         *inp1ix = fstart[monidx] + jobincidx*deg1;
         *inp2ix = fstart[monidx] + jobincidx*deg1 + offset;
      }
      else
      {
         *inp1ix = fstart[monidx] + jobincidx*deg1 + totcffoffset;
         *inp2ix = fstart[monidx] + jobincidx*deg1 + totcffoffset + offset;
      }
   }
   if(jobinkind == 2)   // incrementing backward product
   {
      if(jobisreal)
      {
         *inp1ix = bstart[monidx] + jobincidx*deg1;
         *inp2ix = bstart[monidx] + jobincidx*deg1 + offset;
      }
      else
      {
         *inp1ix = bstart[monidx] + jobincidx*deg1 + totcffoffset;
         *inp2ix = bstart[monidx] + jobincidx*deg1 + totcffoffset + offset;
      }
   }
   if(jobinkind == 3)   // incrementing cross product
   {
      if(jobisreal)
      {
         *inp1ix = cstart[monidx] + jobincidx*deg1;
         *inp2ix = cstart[monidx] + jobincidx*deg1 + offset;
      }
      else
      {
         *inp1ix = cstart[monidx] + jobincidx*deg1 + totcffoffset;
         *inp2ix = cstart[monidx] + jobincidx*deg1 + totcffoffset + offset;
      }
   }
   *outidx = *inp1ix;

   if(verbose)
   {
      cout << "-> inp1ix : " << *inp1ix
           << ", inp2ix : " << *inp2ix
           << ", outidx : " << *outidx << endl;
   }
}

void convjobs_coordinates
 ( ConvolutionJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose )
{ 
   for(int i=0; i<jobs.get_layer_count(layer); i++)
      convjob_indices(jobs.get_job(layer,i),&inp1ix[i],&inp2ix[i],&outidx[i],
                      dim,nbr,deg,nvr,fstart,bstart,cstart,verbose);
}

void complex_convjobs_coordinates
 ( ComplexConvolutionJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr, int totcff, int offset,
   int *fstart, int *bstart, int *cstart, bool verbose )
{ 
   for(int i=0; i<jobs.get_layer_count(layer); i++)
      complex_convjob_indices
         (jobs.get_job(layer,i),&inp1ix[i],&inp2ix[i],&outidx[i],
          dim,nbr,deg,nvr,totcff,offset,fstart,bstart,cstart,verbose);
}

void complex_incjobs_coordinates
 ( ComplexIncrementJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr, int totcff, int offset,
   int *fstart, int *bstart, int *cstart, bool verbose )
{ 
   for(int i=0; i<jobs.get_layer_count(layer); i++)
      complex_incjob_indices
         (jobs.get_job(layer,i),&inp1ix[i],&inp2ix[i],&outidx[i],
          dim,nbr,deg,nvr,totcff,offset,fstart,bstart,cstart,verbose);
}

void addjob_indices
 ( AdditionJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int adtype = job.get_addition_type();
   const int intype = job.get_increment_type();
   const int updmon = job.get_update_monomial();
   const int updidx = job.get_update_index();
   const int incmon = job.get_increment_monomial();
   const int incidx = job.get_increment_index();
   const int deg1 = deg+1;

   if(verbose)
   {
      cout << "  add type : " << adtype
           << ", mon : " << updmon << ", idx : " << updidx;
      cout << "  inc type : " << intype
           << ", mon : " << incmon << ", idx : " << incidx;
      cout << endl;
   }
   if(adtype == 1)
   {
      *inp1ix = fstart[updmon] + updidx*deg1;
   }
   else if(adtype == 2)  // on GPU, one backward item less
   {
      if(updidx == 0)
         *inp1ix = bstart[updmon];
      else
         *inp1ix = bstart[updmon] + (updidx-1)*deg1;
   }
   else if(adtype == 3)
   {
      *inp1ix = cstart[updmon] + updidx*deg1;
   }
   *outidx = *inp1ix;
   if(incmon < 0)
   {
      if(incidx < 0)
         *inp2ix = 0; // start with constant coefficient
      else
         *inp2ix = (1 + incidx)*deg1;
   }
   else
   {
      if(intype == 1)
      {
         *inp2ix = fstart[incmon] + incidx*deg1;
      }
      else if(intype == 2)
      {                                // on GPU, on backward item less
         if(incidx == 0)
            *inp2ix = bstart[incmon];
         else
            *inp2ix = bstart[incmon] + (incidx-1)*deg1;
      }
      else if(intype == 3)
      {
         *inp2ix = cstart[incmon] + incidx*deg1;
      }
   }
   if(verbose)
   {
      cout << "-> inp1ix : " << *inp1ix
           << ", inp2ix : " << *inp2ix
           << ", outidx : " << *outidx << endl;
   }
}

void complex_addjob_indices
 ( ComplexAdditionJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr, int totcff, int offset,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int adtype = job.get_addition_type();
   const int intype = job.get_increment_type();
   const int updmon = job.get_update_monomial();
   const int updidx = job.get_update_index();
   const int incmon = job.get_increment_monomial();
   const int incidx = job.get_increment_index();
   const int deg1 = deg+1;
   const int totcffoffset = totcff + offset; // start of imag data array

   if(verbose)
   {
      cout << "  add type : " << adtype
           << ", mon : " << updmon << ", idx : " << updidx;
      cout << "  inc type : " << intype
           << ", mon : " << incmon << ", idx : " << incidx;
      cout << endl;
   }
   if(adtype == 1)       // real part of forward
   {
      *inp1ix = fstart[updmon] + updidx*deg1;
   }
   else if(adtype == 4)  // imaginary part of forward
   {
      *inp1ix = fstart[updmon] + updidx*deg1 + totcffoffset;
   }
   else if(adtype == 2)  // real part of backward
   {
      if(updidx == 0)
         *inp1ix = bstart[updmon];
      else
         *inp1ix = bstart[updmon] + updidx*deg1;
   }
   else if(adtype == 5)  // imaginary part of backward
   {
      if(updidx == 0)
         *inp1ix = bstart[updmon] + totcffoffset;
      else
         *inp1ix = bstart[updmon] + updidx*deg1 + totcffoffset;
   }
   else if(adtype == 3)  // real part of cross
   {
      *inp1ix = cstart[updmon] + updidx*deg1;
   }
   else if(adtype == 6)  // imaginary part of cross
   {
      *inp1ix = cstart[updmon] + updidx*deg1 + totcffoffset;
   }
   else
      cout << adtype << " is an invalid addition type!" << endl;

   *outidx = *inp1ix;

   if(incmon == -1)   // real part of the constant or coefficient
   {
      if(incidx == -1)
         *inp2ix = 0;             // start with constant coefficient
      else if(incidx == -2)
         *inp2ix = totcffoffset;  // imag part of constant coefficient
      else                        // incidx >= 0 is index of monomial
         *inp2ix = (1 + incidx)*deg1;
   }
   else if(incmon == -2)   // imaginary part of the constant or coefficient
   {
      if(incidx < 0)
         *inp2ix = totcffoffset; // start with constant coefficient
      else                       // incidx >= 0 is index of monomial
         *inp2ix = totcffoffset + (1 + incidx)*deg1;
   }
   else
   {
      if(intype == 1)      // first real operand of forward
      {
         *inp2ix = fstart[incmon] + incidx*deg1;
      }
      else if(intype == 4)      // second real operand of forward
      {
         *inp2ix = fstart[incmon] + incidx*deg1 + offset;
      }
      else if(intype == 7)      // first imaginary operand of forward
      {
         *inp2ix = fstart[incmon] + incidx*deg1 + totcffoffset;
      }
      else if(intype == 10)     // second imaginary operand of forward
      {
         *inp2ix = fstart[incmon] + incidx*deg1 + totcffoffset + offset;
      }
      else if(intype == 2) // first real operand of backward
      {
         if(incidx == 0)
            *inp2ix = bstart[incmon];
         else
            *inp2ix = bstart[incmon] + incidx*deg1;
      }
      else if(intype == 5) // second real operand of backward
      {
         if(incidx == 0)
            *inp2ix = bstart[incmon] + offset;
         else
            *inp2ix = bstart[incmon] + incidx*deg1 + offset;
      }
      else if(intype == 8) // first imaginary operand of backward
      {
         if(incidx == 0)
            *inp2ix = bstart[incmon] + totcffoffset;
         else
            *inp2ix = bstart[incmon] + incidx*deg1 + totcffoffset;
      }
      else if(intype == 11) // second imaginary operand of backward
      { 
         if(incidx == 0)
            *inp2ix = bstart[incmon] + totcffoffset + offset;
         else
            *inp2ix = bstart[incmon] + incidx*deg1
                                     + totcffoffset + offset;
      }
      else if(intype == 3)  // first real operand of cross
      {
         *inp2ix = cstart[incmon] + incidx*deg1;
      }
      else if(intype == 6)  // second real operand of cross
      {
         *inp2ix = cstart[incmon] + incidx*deg1 + offset;
      }
      else if(intype == 9)  // first imaginary operand of cross
      {
         *inp2ix = cstart[incmon] + incidx*deg1 + totcffoffset;
      }
      else if(intype == 12) // second imaginary operand of cross
      {
         *inp2ix = cstart[incmon] + incidx*deg1 + totcffoffset + offset;
      }
      else
         cout << intype << " is an invalid increment type!" << endl;
   }
   if(verbose)
   {
      cout << "-> inp1ix : " << *inp1ix
           << ", inp2ix : " << *inp2ix
           << ", outidx : " << *outidx << endl;
   }
}

void addjobs_coordinates
 ( AdditionJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   for(int i=0; i<jobs.get_layer_count(layer); i++)
      addjob_indices(jobs.get_job(layer,i),&inp1ix[i],&inp2ix[i],&outidx[i],
                     dim,nbr,deg,nvr,fstart,bstart,cstart,verbose);
}

void complex_addjobs_coordinates
 ( ComplexAdditionJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr, int totcff, int offset,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   for(int i=0; i<jobs.get_layer_count(layer); i++)
      complex_addjob_indices
         (jobs.get_job(layer,i),&inp1ix[i],&inp2ix[i],&outidx[i],
          dim,nbr,deg,nvr,totcff,offset,fstart,bstart,cstart,verbose);
}
