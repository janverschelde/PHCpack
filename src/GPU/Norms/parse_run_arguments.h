// Contains the prototype to parse the command line arguments
// in all manual test functions.

#ifndef __parse_run_arguments_h__
#define __parse_run_arugments_h__

int parse_run_arguments
 ( int argc, char *argv[], int *blocksize, int *dim, int *freq, int *mode );
/*
   Parses the argc arguments on the command line in argv[]
   for the run_ manual test functions.
   Returns 0 if four numbers are given, 1 otherwise for failure.

   ON RETURN :
     blocksize  the block size, number of threads in a block;
     dim        dimension of the vectors;
     freq       frequency of the runs;
     mode       execution mode is 0, 1, or 2
                0 : GPU run, no output,
                1 : CPU run, no output,
                2 : both GPU and CPU run with output. */

#endif
