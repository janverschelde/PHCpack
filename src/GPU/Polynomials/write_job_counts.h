/* Utilities to write the counts of convolution and addition jobs. */

#ifndef __write_job_counts_h__
#define __write_job_counts_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

void write_convolution_counts ( ConvolutionJobs jobs );
/*
 * DESCRIPTION :
 *   Writes the counts of convolution jobs in each layer. */

void write_addition_counts ( AdditionJobs jobs );
/*
 * DESCRIPTION :
 *   Writes the counts of addition jobs in each layer. */

void convolution_operation_counts
 ( int deg, ConvolutionJobs cnvjobs, long int *addcnt, long int *mulcnt,
   int vrblvl );
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

void write_operation_counts
 ( int deg, ConvolutionJobs cnvjobs, AdditionJobs addjobs );
/*
 * DESCRIPTION :
 *   Writes the number of operations for power series truncated
 *   to degree deg, for all convolution and addition jobs. */

#endif
