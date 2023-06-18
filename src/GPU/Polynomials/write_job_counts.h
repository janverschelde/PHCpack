/* Utilities to write the counts of convolution and addition jobs. */

#ifndef __write_job_counts_h__
#define __write_job_counts_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"
#include "complexadd_jobs.h"

void write_convolution_counts ( ConvolutionJobs jobs );
/*
 * DESCRIPTION :
 *   Writes the counts of convolution jobs in each layer. */

void write_addition_counts ( AdditionJobs jobs );
/*
 * DESCRIPTION :
 *   Writes the counts of addition jobs in each layer. */

void write_complexconv_counts ( ComplexConvolutionJobs jobs );
/*
 * DESCRIPTION :
 *   Writes the counts of convolution jobs in each layer. */

void write_complexadd_counts ( ComplexAdditionJobs jobs );
/*
 * DESCRIPTION :
 *   Writes the counts of addition jobs in each layer. */

void convolution_operation_counts
 ( int deg, ConvolutionJobs cnvjobs,
   long long int *addcnt, long long int *mulcnt, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the total number of arithmetical operations defined
 *   by the convolution jobs for series truncated at degree deg.
 *   Writes the floating-point operation counts if vrblvl > 0.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series;
 *   cnvjobs  defines the convolution jobs;
 *   vrblvl   if zero then no output, else writes the counts.
 *
 * ON RETURN :
 *   addcnt   total number of additions;
 *   mulcnt   total number of multiplications. */

void complexconv_operation_counts
 ( int deg, ComplexConvolutionJobs cnvjobs,
   long long int *addcnt, long long int *mulcnt, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the total number of arithmetical operations defined
 *   by the convolution jobs for series truncated at degree deg.
 *   Writes the floating-point operation counts if vrblvl > 0.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series;
 *   cnvjobs  defines the convolution jobs;
 *   vrblvl   if zero then no output, else writes the counts.
 *
 * ON RETURN :
 *   addcnt   total number of additions;
 *   mulcnt   total number of multiplications. */

void complexinc_operation_counts
 ( int deg, ComplexIncrementJobs incjobs,
   long long int *addcnt, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the total number of addition operations defined
 *   by the increment jobs for series truncated at degree deg.
 *   Writes the floating-point operation counts if vrblvl > 0.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series;
 *   incjobs  defines the increment jobs;
 *   vrblvl   if zero then no output, else writes the counts.
 *
 * ON RETURN :
 *   addcnt   total number of additions. */

void write_operation_counts
 ( int deg, ConvolutionJobs cnvjobs, AdditionJobs addjobs );
/*
 * DESCRIPTION :
 *   Writes the number of operations for power series truncated
 *   to degree deg, for all convolution and addition jobs. */

void write_complexop_counts
 ( int deg, ComplexConvolutionJobs cnvjobs, ComplexAdditionJobs addjobs );
/*
 * DESCRIPTION :
 *   Writes the number of operations for power series truncated
 *   to degree deg, for all convolution and addition jobs,
 *   using complex vectorized arithmetic. */

#endif
