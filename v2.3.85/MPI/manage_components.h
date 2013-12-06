/* The file "manage_components.h" collects macro and type definitions,
 * along with prototypes of functions to manage components for the
 * parallel monodromy breakup algorithm using dynamic load balancing. */

#define SEND_TRACK_JOB 104
#define SEND_QUIT 105
#define SEND_END_POINT 106

int MIN(int X, int Y);

int SWAP(int* X, int* Y);
/*
 * DESCRIPTION :
 *   Swaps the content of X with Y. */

typedef int POINT; /* D: this is a dummy type */ 
typedef int bool; 

struct COMPONENT_STRUCT
{
  int num;
  int size;
  int* pts;
  bool is_irred;
  int n_nodes_tracking;
};
typedef struct COMPONENT_STRUCT COMPONENT;

void makeCOMPONENT(COMPONENT* c, int n, int s);
/*
 * DESCRIPTION :
 *   Creates a component with the content of the parameters. */

void killCOMPONENT(COMPONENT* c);
/*
 * DESCRIPTION :
 *   Destroys the component c. */

bool is_in_COMPONENT(int a, COMPONENT* c);
/*
 * DESCRIPTION :
 *   Returns 1 if the point a belongs to the component c,
 *   returns 0 otherwise. */

void mergeCOMPONENT(COMPONENT* a, COMPONENT* b, int* ind);
/*
 * DESCRIPTION :
 *   Merges "a" and "b", puts the result in "a", kills "b", 
 *   modifes the index "ind". */

/******* state of a slave ***********************************************/

struct SLAVE_STATE_STRUCT
{
   int num;            /* slave number */
   int  status; 
#define FREE 0
#define TRACKING 1 
#define NEW_POINT 2 /* ... NEW_POINT = tracking the label 
		       that is not in the target slice */
   int sliceA, sliceB; /* if yes, then slices it works with */ 
   int start_label;    /* point it is tracing */
};
typedef struct SLAVE_STATE_STRUCT SLAVE_STATE;

int free_slave(SLAVE_STATE* state, int np);
/*
 * DESCRIPTION :
 *   Returns the number of the first free slave,
 *   or 0 if all the nodes are busy. */

/******* interslice statistics ******************************************/

#define N_PATHS_TO_REMEMBER 100
struct ISLICE_STRUCT
{
  int n_paths;
  int n_rejected;
  int path_outcome[N_PATHS_TO_REMEMBER];
#define PATH_LED_TO_MERGER 0
#define PATH_DIVERGED 1
#define PATH_STAYED_WITHIN_COMPONENT 2 
#define PATH_CROSSED 3 
};
typedef struct ISLICE_STRUCT ISLICE_STATS;

void record_path_outcome(ISLICE_STATS* s, int o);
/*
 * DESCRIPTION :
 *   Records in "s" the outcome "o" of a path. */ 

 
#define MAX_PERCENT_REJECTED 100
#define MIN_ATTEMPTS 10

bool isGood(ISLICE_STATS* s); 
/*
 * DESCRIPTION :
 *   Determines whether the pair of slices should be used. */ 

void print_component(COMPONENT* c);
/*
 * DESCRIPTION :
 *   Print the component "c" of the witness set. */

void print_partition(int n, COMPONENT* part);
/*
 * DESCRIPTION :
 *   Print the current partition of the witness set. */

/*********** functions that depend on PHC ******************************/

bool is_in_slice(int pt, int slice);
/*
 * DESCRIPTION :
 *   Is label "pt" assigned to a point in the "slice". */

#define TRACE_TOLERANCE 0.00000001
bool is_irreducible(COMPONENT* c);
/*
 * DESCRIPTION :
 *   Returns true if the trace sum difference is less than the
 *   value of TRACE_TOLERANCE. */
