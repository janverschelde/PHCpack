int read_two_witness_sets 
      ( int *n1, int *n2, int *dim1, int *dim2, int *deg1, int *deg2,
        int *cd );
/*
 * DESCRIPTION :
 *   The user is prompted for two witness sets.
 *
 * ON RETURN :
 *   n1        ambient dimension for the first witness set;
 *   n2        ambient dimension for the second witness set;
 *   dim1      dimension of the solutions represented by the 1st witness set;
 *   dim2      dimension of the solutions represented by the 2nd witness set;
 *   deg1      degree of the 1st witness set, i.e.: #solutions in 1st set;
 *   deg2      degree of the 2nd witness set, i.e.: #solutions in 2nd set;
 *   cd        cascade dimension: #vars and #eqs in the homotopy. */

 
int track_one_path ( int n, int i, int *m, double *sol );
/*
 * DESCRIPTION :
 *   Tracks one path starting at the given start solution.
 *
 * ON ENTRY :
 *   n        dimension of the solution;
 *   i        number of the solution;
 *   m        label of the solution, multiplicity flag;
 *   sol      coordinates of the solution vector. */

void read_dimension_of_system ( int nc, char name[nc], int *n );
/*
 * DESCRIPTION :
 *   Opens the file with the given name and reads the first number
 *   which is returned in n.
 *
 * ON ENTRY :
 *   nc       number of characters in the name;
 *   name     name of an input file which contains a system.
 *
 * ON RETURN :
 *   n        the dimension of the system on the file. */

