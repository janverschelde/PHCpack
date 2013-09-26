#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "../Lib/phcpack.h"
#include "../Lib/syscon.h"
#include "../Lib/solcon.h"
#include "parallel_phcpack.h"

#define SEND_SSOL 100 /* tag for sending start solution */
#define SEND_SMUL 101 /* tag for sending multiplicity of start solution */
#define SEND_TSOL 102 /* tag for sending target solution */
#define SEND_TMUL 103 /* tag for sending multiplicity of target solution */

#define v 0 

int read_two_witness_sets
      ( int *n1, int *n2, int *dim1, int *dim2, int *deg1, int *deg2,
        int *cd )
{
   int fail;

   printf("\n");
   fail = read_a_witness_set(1,n1,dim1,deg1);
   printf("  n = %d  dimension = %d  degree = %d\n",*n1,*dim1,*deg1);
   printf("\n");
   fail = read_a_witness_set(2,n2,dim2,deg2);
   printf("  n = %d  dimension = %d  degree = %d\n",*n2,*dim2,*deg2);

   fail = extrinsic_top_diagonal_dimension(*n1,*n2,*dim1,*dim2,cd);
   printf("The top dimension of the extrinsic diagonal cascade : %d.\n",*cd);

   return fail;
}

int track_one_path ( int n, int i, int *m, double *sol )
{
   int fail,nbstep,nbfail,nbiter,nbsyst;

   fail = silent_path_tracker
            (n,m,sol,&nbstep,&nbfail,&nbiter,&nbsyst);
   if(v>0) printf("%d : #step : %3d #fail : %3d #iter : %3d #syst : %3d\n",
                  i,nbstep,nbfail,nbiter,nbsyst);

   return fail;
}

void read_dimension_of_system ( int nc, char name[nc], int *n )
{
   FILE *fp;

   *n = 0;
   
   if(v>0) printf("opening file %s\nto read the dimension ...\n",name);

   fp = fopen(name,"r");
   fscanf(fp,"%d",n);

   if(v>0) printf("found the dimension n = %d\n",*n);

   fclose(fp);
}
