// The file job_coordinates.h specifies functions to compute the positions
// of data to executed convolution and addition jobs on the device.

#ifndef __job_coordinates_h__
#define __job_coordinates_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

int coefficient_count ( int dim, int nbr, int deg, int *nvr );
/*
 * DESCRIPTION :
 *   Returns the total number of coefficients in all convolutions.
 *   The count includes also the constant coefficient.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k. */

void coefficient_indices
 ( int dim, int nbr, int deg, int *nvr,
   int *fsums, int *bsums, int *csums,
   int *fstart, int *bstart, int *cstart );
/*
 * DESCRIPTION :
 *   Computes the sums of coefficients for forward, backward, cross products,
 *   and the start positions for each monomial.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k.
 *
 * ON RETURN :
 *   fsums    fsums[k] holds the sum of coefficients for the forward
 *            products of the k-th monomial and of all monomials before k;
 *   bsums    fsums[k] holds the sum of coefficients for the backward
 *            products of the k-th monomial and of all monomials before k;
 *   csums    fsums[k] holds the sum of coefficients for the cross
 *            products of the k-th monomial and of all monomials before k;
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial. */

void convjob_indices
 ( ConvolutionJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose );
/*
 * DESCRIPTION :
 *   Computes the indices of the two inputs and the output of a job.
 *
 * ON ENTRY :
 *   job      defines a convolution job;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   verbose  if true, writes extra information about the job.
 *
 * ON RETURN :
 *   inp1ix   index of the first input;
 *   inp2ix   index of the second input;
 *   outidx   index of the output. */

void convjobs_coordinates
 ( ConvolutionJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose );
/*
 * DESCRIPTION :
 *   Defines the coordinates of all jobs in the same layer.
 *
 * ON ENTRY :
 *   jobs     defines convolution jobs;
 *   layer    the index of one layer of jobs;
 *   inp1ix   space for as many integers as the jobs on the layer;
 *   inp2ix   space for as many integers as the jobs on the layer;
 *   outidx   space for as many integers as the jobs on the layer;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   verbose  if true, writes extra information about the jobs.
 *
 * ON RETURN :
 *   inp1ix   inp1ix[i] is the index of the first input of job i;
 *   inp2ix   inp2ix[i] is the index of the second input of job i;
 *   outidx   outidx[i] is the index of the output of job i. */

void addjob_indices
 ( AdditionJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose );
/*
 * DESCRIPTION :
 *   Computes the indices of the two inputs and the output of a job.
 *
 * ON ENTRY :
 *   job      defines an addition job;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   verbose  if true, writes extra information about the job.
 *
 * ON RETURN :
 *   inp1ix   index of the first input;
 *   inp2ix   index of the second input;
 *   outidx   index of the output. */

void addjobs_coordinates
 ( AdditionJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose );
/*
 * DESCRIPTION :
 *   Defines the coordinates of all jobs in the same layer.
 *
 * ON ENTRY :
 *   jobs     defines addition jobs;
 *   layer    the index of one layer of jobs;
 *   inp1ix   space for as many integers as the jobs on the layer;
 *   inp2ix   space for as many integers as the jobs on the layer;
 *   outidx   space for as many integers as the jobs on the layer;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   verbose  if true, writes extra information about the jobs.
 *
 * ON RETURN :
 *   inp1ix   inp1ix[i] is the index of the first input of job i;
 *   inp2ix   inp2ix[i] is the index of the second input of job i;
 *   outidx   outidx[i] is the index of the output of job i. */

#endif
