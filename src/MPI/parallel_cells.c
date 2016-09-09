/* The file "parallel_cells.c" collects the definition of the routines
 * with prototypes in "parallel_cells.h". */

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "parallel_phcpack.h"
#include "../Lib/celcon.h"

#define v 2  /* verbose flag:
                  0 no output during the computations,
                  1 only one line messages on progress of computations
                  2 display of data, like systems and solutions */

void retrieve_dimensions ( int myid, int *nspt, int *dim )
{
   int fail, *mix;
   
   if(myid == 0)
   {  
      fail = celcon_read_mixed_cell_configuration();
      printf("\nReading a system to initialize the symbol table...");
      fail = read_standard_target_system();
      fail = define_output_file();
      fail = celcon_dimension_of_points(dim);
      mix = (int*)calloc((*dim-1), sizeof(int));
      fail = celcon_type_of_mixture(nspt,mix);
      if(v>0) printf("\nThe number of different support is %d\n",*nspt);
      if(v>0) printf("The dimension of lifted points is %d\n",*dim);
   }
   
   MPI_Bcast(nspt,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(dim,1,MPI_INT,0,MPI_COMM_WORLD);   
}

void supports_broadcast ( int myid, int nspt, int dim )
{
   int mix[nspt];
   int fail,i,j;
   double point[dim];
 
   if(myid == 0)
   { 
      fail = celcon_type_of_mixture(&nspt,mix);  
      if(v>0) printf("The number of occurrences of supports is: ");
      if(v>0) Print_Integer_Array(nspt,mix);
   }
    
   MPI_Bcast(mix,nspt,MPI_INT,0,MPI_COMM_WORLD);
   if(myid != 0)
      fail = celcon_set_type_of_mixture(nspt,mix);  

   if(myid == 0)
   {
      fail = celcon_length_of_supports(&nspt,mix); 
      /* get #different supports and #points of each support. */
      if(v>0) printf("The number of points in root supports is: ");
      if(v>0) Print_Integer_Array(nspt,mix);
   }
   MPI_Bcast(mix,nspt,MPI_INT,0,MPI_COMM_WORLD);
    
   for(i=0;i<nspt;i++)
      for(j=0;j<mix[i];j++)
      {
         if(myid==0) fail = celcon_get_lifted_point(dim,i+1,j+1,point);
         MPI_Bcast(point,dim,MPI_DOUBLE,0,MPI_COMM_WORLD);
         if(myid!=0) fail = celcon_append_lifted_point(dim,i+1,point);
      }
}

void system_broadcast ( int myid, int n )
{
   int fail;

   if(myid == 0)
   {
      fail = celcon_standard_random_coefficient_system(); 
   /* if(v>1) fail = celcon_write_standard_random_coefficient_system(); */
      fail = celcon_copy_into_standard_systems_container();
   }
   dimension_broadcast(myid,&n);
   monomials_broadcast(myid,n);
   if(myid != 0)
   {
      fail = celcon_copy_from_standard_systems_container();
   /* if(v>1) if(myid==1)
         fail = celcon_write_standard_random_coefficient_system(); */
   }
}
