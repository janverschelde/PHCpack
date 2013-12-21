/* The file "parallel_phcpack.h" collects protypes for use
 * in parallel path tracking with PHCpack. */

#define SEND_SSOL 100 /* tag for sending start solution */
#define SEND_SMUL 101 /* tag for sending multiplicity of start solution */
#define SEND_TSOL 102 /* tag for sending target solution */
#define SEND_TMUL 103 /* tag for sending multiplicity of target solution */

void dimension_broadcast ( int myid, int *n );
/*
 * DESCRIPTION :
 *   After broadcasting the ambient dimension n,
 *   every node intializes its system container with n.
 *
 * ON ENTRY :
 *   myid     number of the node;
 *   *n       ambient dimension of the problem. */

void monomials_broadcast ( int myid, int n );
/*
 * DESCRIPTION :
 *   The system container at the root node is broadcasted to all nodes.
 *
 * ON ENTRY :
 *   myid     number of the node;
 *   n        ambient dimension. */

int start_system_broadcast ( int myid, int n, int *nbsols );
/*
 * DESCRIPTION :
 *   The manager prompts the user for a start system with solutions.
 *   The start system is broadcasted to all the nodes.
 *
 * ON ENTRY :
 *   myid     number of the node;
 *   n        number of polynomials in the system.
 *
 * ON RETURN :
 *   nbsols   number of start solutions. */

int start_system_broadcast_without_solutions ( int myid, int n, int *nbsols );
/*
 * DESCRIPTION :
 *   The manager prompts the user for a start system and for a file for the
 *   start solutions.  The start system is broadcasted to all the nodes.
 *   Only the number of solutions is written from file.
 *
 * ON ENTRY :
 *   myid     number of the node;
 *   n        number of polynomials in the system.
 *
 * ON RETURN :
 *   nbsols   number of start solutions. */

int named_start_system_broadcast
    ( int myid, int n, int *nbsols, int kind, int nc, char name[nc] );
/*
 * DESCRIPTION :
 *   The start system on file with the given name is read and broadcasted.
 *
 * ON ENTRY :
 *    myid    id label of the node;
 *    n       number of polynomials in the system;
 *    kind    type of start system or homotopy;
 *    nc      number of characters in the file name;
 *    name    file name.
 *
 * ON RETURN :
 *    nbsols  number of start solution, in case kind == 3. */

int start_in_container_broadcast ( int myid, int n );
/* 
 * DESCRIPTION :
 *   The manager node has a start system which gets copied to the systems
 *   container and broadcasted along the network.  All nodes copy then the
 *   system from the systems container to their start system. */

int target_in_container_broadcast ( int myid, int n );
/* 
 * DESCRIPTION :
 *   The manager node has a target system which gets copied to the systems
 *   container and broadcasted along the network.  All nodes copy then the
 *   system from the systems container to their target system. */

int homotopy_broadcast ( int myid, int n );
/*
 * DESCRIPTION :
 *   Start and target system at the manager node are broadcasted 
 *   to all nodes.  All nodes then initialize their homotopy. */

void solutions_broadcast ( int myid, int nbsols, int n );
/*
 * DESCRIPTION :
 *   The solutions container at node 0 are broadcasted to all other nodes.
 *
 * ON ENTRY :
 *   myid     number of the node;
 *   nbsols   number of solutions in the container at node 0;
 *   n        dimension of the vectors in the solution list. */

void solutions_distribute ( int myid, int nbsols, int n, int nprocs,
                            int *solnum );
/*
 * DESCRIPTION :
 *   Distributes the solutions in the container at the root
 *   to the solutions container at each node.
 *
 * ON ENTRY : 
 *   myid     number of the node;
 *   nbsols   total number of solutions to distribute;
 *   n        length of the solution vectors;
 *   nprocs   the number of nodes;
 *
 * ON RETURN : 
 *   solnum   number of solutions assigned to the node. */

void solutions_collect ( int myid, int nbsols, int n, 
                         int numprocs, int mysolnum );
/* 
 * DESCRIPTION :
 *   Collects all target solutions in the container at the root
 *   from the solutions container at each node. */

void print_monomials ( void );
/*
 * DESCRIPTION :
 *   Writes the monomials in the container to screen. */

void write_solution_banner_to_defined_output_file ( int nbsols, int n );
/*
 * DESCRIPTION :
 *   Writes the banner "THE SOLUTIONS" followed by the number of solutions
 *   in nbsols and the dimension n to the defined output file.
 *   One extra banner is then written to the defined output file. */

void print_solutions ( int myid );
/* 
 * DESCRIPTION :
 *   Prints the solutions in the container. */

void print_time ( double *time, int numprocs );
/* 
 * DESCRIPTION :
 *   Prints the time spent on each of the processor. */

void write_time_and_paths_to_defined_output_file
       ( int p, double time[p], int paths[p] );
/*
 * DESCRIPTION :
 *   Writes the time and the number of paths tracked by every processor
 *   to the defined output file.
 *
 * ON ENTRY :
 *   p        number of processors used for the path tracking;
 *   time     wall time for every processor;
 *   paths    number of paths tracked by every processor. */
