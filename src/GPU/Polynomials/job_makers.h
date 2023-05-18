/* The job_makers.h contains the prototypes of functions to make
 * all jobs for the polynomial evaluation and differentiation. */

#ifndef __job_makers_h__
#define __job_makers_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"
#include "complexadd_jobs.h"

void make_all_jobs
 ( int dim, int nbr, int *nvr, int **idx,
   ConvolutionJobs *cnvjobs, AdditionJobs *addjobs, bool verbose );
/*
 * DESCRIPTION :
 *   Defines all convolution and addition jobs for a polynomial.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nvr      nvr[k] has the number of variables in monomial k;
 *   idx      idx[k] holds nvr[k] indices to variables in monomial k;
 *   verbose  is the verbose flag, if true, then information about
 *            all jobs is written to screen.
 *
 * ON RETURN :
 *   cnvjobs  are the convolution jobs;
 *   addjobs  are the addition jobs. */

void make_all_complex_jobs
 ( int dim, int nbr, int *nvr, int **idx,
   ComplexConvolutionJobs *cnvjobs, ComplexIncrementJobs *incjobs,
   ComplexAdditionJobs *addjobs, bool verbose );
/*
 * DESCRIPTION :
 *   Defines all convolution and addition jobs for a polynomial,
 *   using complex vectorized arithmetic.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nvr      nvr[k] has the number of variables in monomial k;
 *   idx      idx[k] holds nvr[k] indices to variables in monomial k;
 *   verbose  is the verbose flag, if true, then information about
 *            all jobs is written to screen.
 *
 * ON RETURN :
 *   cnvjobs  are the convolution jobs;
 *   incjobs  are the increment jobs;
 *   addjobs  are the addition jobs. */

void write_jobs_report
 ( int dim, int nva, int nbr, int deg,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs );
/*
 * DESCRIPTION :
 *   Writes the dimensions of the jobs.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nva      number of variables per monomial (if fixed);
 *   nbr      number of monomials, excluding the constant;
 *   deg      truncation degree of the series;
 *   cnvjobs  are the convolution jobs;
 *   addjobs  are the addition jobs. */

void write_complex_jobs_report
 ( int dim, int nva, int nbr, int deg,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   ComplexAdditionJobs addjobs );
/*
 * DESCRIPTION :
 *   Writes the dimensions of the jobs,
 *   using vectorized complex arithmetic.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nva      number of variables per monomial (if fixed);
 *   nbr      number of monomials, excluding the constant;
 *   deg      truncation degree of the series;
 *   cnvjobs  are the convolution jobs;
 *   incjobs  are the increment jobs;
 *   addjobs  are the addition jobs. */

void write_CPU_timings ( double lapsec1, double lapsec2 );
/*
 * DESCRIPTION :
 *   Writes the timings.
 *
 * ON ENTRY :
 *   lapsec1  is the CPU time for the first algorithm, without jobs
 *   lapsec2  is the CPU time for the second algorithm, with jobs. */

#endif
