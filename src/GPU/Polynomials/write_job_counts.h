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

void write_operation_counts
 ( int deg, ConvolutionJobs cnvjobs, AdditionJobs addjobs );
/*
 * DESCRIPTION :
 *   Writes the number of operations for power series truncated
 *   to degree deg, for all convolution and addition jobs. */

#endif
