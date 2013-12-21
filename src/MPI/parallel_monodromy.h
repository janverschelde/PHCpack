/* The file "parallel_monodromy.h" collects prototypes for functions to
 * decompose a pure dimensional solution set into irreducible factors. */

void slices_broadcast ( int myid, int k, int n );
/*
 * DESCRIPTION :
 *   The root node generates k random slides in n-space,
 *   which are then broadcasted to all nodes.  */

void trace_slice_broadcast ( int myid, int first );
/* 
 * DESCRIPTION :
 *   Sets the coefficient of the slice used in the linear trace,
 *   for the first time, first must have the value 1. */

void gamma_broadcast ( int myid, int n );
/* 
 * DESCRIPTION :
 *   The root node generates n random complex numbers,
 *   which are broadcasted to all nodes.  */

int witness_set_distribute ( int myid, int n, int dim, int deg, int np );
/*
 * DESCRIPTION :
 *   The system defining a witness set in the systems container is
 *   broadcasted to all nodes and everybody initializes the sampler.
 *   The solution set is distributed evenly among the processors.
 *
 * ON ENTRY :
 *   myid          identification of the node;
 *   n             ambient dimension;
 *   dim           dimension of the solution set;
 *   deg           degree of the solution set;
 *   np            number of processors.
 *
 * ON RETURN :
 *   the fail code of the initialize_sampler routine. */

void f_swap_slices ( int myid );
/* 
 * DESCRIPTION :
 *   Swaps the current slices with new slices and takes new solutions
 *   as start to turn back. */

void f_track_paths ( int myid, int deg, int n, int numprocs, int mysolnum,
                     double *trackwtime );
/* 
 * DESCRIPTION :
 *   All nodes track the paths defined by the homotopy and the start
 *   solutions.  The root node collects the target solutions. */

int build_trace_grid ( int myid, int n, int dim, int deg,
                       int numprocs, double *trackwtime );
/* 
 * DESCRIPTION :
 *   A trace grid consists of 3 witness sets on parallel slices,
 *   need to provide a linear trace certificate for a factor.
 *
 * ON ENTRY :
 *   myid          identification of the node;
 *   n             ambient dimension;
 *   dim           dimension of the solution set;
 *   deg           degree of the solution set;
 *   numprocs      number of processors.
 *
 * ON RETURN :
 *   trackwtime    time node myid has been tracking paths;
 *   the return value of the function is either 0 or 1:
 *     0           if the trace grid is accurate;
 *     1           otherwise: failures occurred. */

void trace_grid_broadcast ( int myid, int d, int n );
/*
 * DESCRIPTION :
 *   The trace grid is sent to every node in the network. 
 *
 * REQUIRED : the parameters of this function have been broadcasted.
 *
 * ON ENTRY :
 *   myid          identification of the node;
 *   d             degree of the solution set;
 *   n             ambient dimension, length of solution vectors. */

void all_slices_broadcast ( int myid, int nb_slices, int dim, int n );
/*
 * DESCRIPTION :
 *   The root node generates a number of slices for broadcasting.
 * 
 * ON ENTRY :
 *   myid          identification of the node;
 *   nb_slices     number of new slices used (in addition to original);
 *   dim           dimension of the solution set;
 *   n             ambient dimension. */

int track_paths_to_new_slices ( int myid, int nb_slices, int d, int n, int k,
                                int numprocs, double *trackwtime );
/*
 * DESCRIPTION :
 *   Computes all solutions on as many slices as the value of nb_slices.
 *   If this routine returns no failure, then every node has d solution
 *   sets on as many hyperplane sections as the value of nb_slices.
 *
 * ON ENTRY :
 *   myid          identification of the node;
 *   nb_slices     number of hyperplane sections;
 *   d             number of solutions in every slice.
 *   k             dimension of the solution set, #planes in every slice;
 *   n             ambient dimension, length of every solution vector;
 *   numprocs      number of processors. 
 *
 * ON RETURN :
 *   trackwtime    time node myid has been tracking paths. */

void interactive_trace_test( int myid );
/* 
 * DESCRIPTION :
 *   Interactive test on the function trace_sum_difference;
 *   if this can be done on any node, then trace_grid_broadcast is okay. */

void interactive_in_slice_test ( int myid );
/*
 * DESCRIPTION :
 *   User can test the function in_slice to see if every node have
 *   the same solutions on every slice. */

void interactive_track_to_make_loops ( int myid );
/*
 * DESCRIPTION :
 *   Interactive program to track paths to make loops. */
