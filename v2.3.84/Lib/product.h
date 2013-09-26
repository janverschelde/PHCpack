/* file product.h contains prototypes for the linear-product
 * root counts and random linear-product systems */

int supporting_set_structure ( void );
/*
 * DESCRIPTION :
 *   Creates a supporting set structure for the system in the container.
 *
 * REQUIRED : the systems container contains a valid polynomial system. */

int write_set_structure ( void );
/*
 * DESCRIPTION :
 *   Writes the supporting set structure. */

int linear_product_root_count ( int *r );
/*
 * DESCRIPTION :
 *   Returns in r the root count based on the supporting set structure.
 *
 * REQUIRED : supporting_set_structure() was executed. */

int random_linear_product_system ( void );
/*
 * DESCRIPTION :
 *   Replaces the system in the systems container with 
 *   a random linear-product system based on the supporting set structure.
 *
 * REQUIRED :
 *   supporting_set_structure() was executed and the systems container
 *   still contains the original polynomial system. */

int solve_linear_product_system ( void );
/*
 * DESCRIPTION :
 *   Puts the solution of the random linear-product system in the
 *   solutions container.
 *
 * REQUIRED :
 *   random_linear_product_system() was executed. */

int clear_set_structure ( void );
/*
 * DESCRIPTION :
 *   Clears the set structure constructed with supporting_set_structure. */
