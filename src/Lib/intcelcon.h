/* This file "intcelcon.h" contains the prototypes of the operations
 * on the cells container for integer lifting functions.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __INTCELCON_H__
#define __INTCELCON_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int intcelcon_read_mixed_cell_configuration ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name, reads a mixed-cell configuration,
 *   which is then stored in the container. */

int intcelcon_write_mixed_cell_configuration ( void );
/*
 * DESCRIPTION :
 *   Writes the mixed-cell configuration in the container to screen. */

int intcelcon_number_of_cells ( int *length );
/*
 * DESCRIPTION :
 *   Returns in length the number of cells in the container. */

int intcelcon_dimension_of_points ( int *dimension );
/*
 * DESCRIPTION :
 *   Returns the dimension of the lifted points in the cells container. */

int intcelcon_type_of_mixture ( int *r, int *mix );
/*
 * DESCRIPTION :
 *   Returns in r the number of different supports and in mix[i] the number
 *   of occurrences of the (i+1)-th support.  By default, call this routine
 *   with mix[] the size of the dimension of the points minus one, as the
 *   size of mix[] on entry must be as large as the number of supports. */

int intcelcon_length_of_supports ( int *r, int *length );
/*
 * DESCRIPTION :
 *   Returns in r the number of different supports and in length[i] the number
 *   of points in the (i+1)-th support.  By default, call this routine
 *   with length[] the size of the dimension of the points minus one, as the
 *   size of length[] on entry must be as large as the number of supports. */

int intcelcon_get_lifted_point ( int n, int i, int j, int *point );
/*
 * DESCRIPTION :
 *   Returns the j-th point of the i-th support.
 *
 * ON ENTRY :
 *   n         length of the lifted vectors in the supports;
 *   i         index to a lifted support;
 *   j         index to a point in the i-th lifted support.
 *
 * ON RETURN :
 *   point     coordinates of the j-th point in the i-th lifted support,
 *             of dimension n. */

int intcelcon_get_inner_normal ( int n, int i, int *normal );
/*
 * DESCRIPTION :
 *   Returns the inner normal to the i-th mixed cell in the container.
 *
 * ON ENTRY :
 *   n         length of the lifted vectors in the supports;
 *   i         index to a cell in the container.
 *
 * ON RETURN :
 *   normal    coordinates of the inner normal to the i-th cell,
 *             of dimension n. */

int intcelcon_number_of_points_in_cell ( int i, int *length );
/*
 * DESCRIPTION :
 *   Returns the number of points in each support of the i-th cell.
 *
 * ON ENTRY :
 *   i         index to a cell in the container.
 *
 * ON RETURN :
 *   length    length[j] is the number of points in the (j+1)-th list
 *             of the supports of cell i, for j in 0..r-1,
 *             where r equals the number of distinct supports. */

int intcelcon_get_point_in_cell ( int n, int i, int j, int k, int *point );
/*
 * DESCRIPTION :
 *   Returns the k-th point from the j-th list of the i-th cell.
 *
 * ON ENTRY :
 *   n         length of the lifted vectors in the supports;
 *   i         index to a cell in the container;
 *   j         index to a support of the i-th cell;
 *   k         index to a point in the j-th support of the i-th cell.
 *
 * ON RETURN :
 *   point     coordinates of the k-th point of the j-th support of cell i,
 *             of dimension n. */

int intcelcon_mixed_volume ( int i, int *mv );
/*
 * DESCRIPTION :
 *   Returns in mv the mixed volume of the i-th cell. */

int intcelcon_initialize_supports ( int nbr );
/*
 * DESCRIPTION :
 *   Initializes the number of distinct supports in the cells container. */

int intcelcon_set_type_of_mixture ( int r, int *mix );
/*
 * DESCRIPTION :
 *   Sets the number of different supports to r and the number of
 *   occurrences of the (i+1)-th support to mix[i], mix has dimension r. */

int intcelcon_append_lifted_point ( int n, int i, int *point );
/*
 * DESCRIPTION :
 *   Appends the point (of dimension n) to the i-th support.
 *   Note that the first support has index 1 (not zero). */

int intcelcon_append_mixed_cell
 ( int n, int r, int k, int *labels, int *normal );
/* 
 * DESCRIPTION :
 *   Appends a mixed cell to the cells container.
 *
 * ON ENTRY :
 *   n         length of the lifted vectors in the supports;
 *   r         number of different supports;
 *   k         length of the vector of labels;
 *   labels    the total number of points in the supports in labels[0],
 *             the number of points in the j-th support is in labels[j],
 *             labels[1+r+j] contains the labels of the points,
 *             the dimension of the labels is k;
 *   normal    coordinates of the inner normal to the cell,
 *             of dimension n. */

int intcelcon_retrieve_mixed_cell
 ( int n, int r, int i, int *labels, int *normal );
/*
 * DESCRIPTION :
 *   Retrieves the i-th mixed cell from the container.
 *
 * ON ENTRY :
 *   n         length of the lifted vectors in the supports;
 *   r         number of different supports;
 *   i         index to a mixed cells in the cells container.
 * 
 * ON RETURN :
 *   labels    the total number of points in the supports in labels[0],
 *             the number of points in the j-th support is in labels[j],
 *             labels[1+r+j] contains the labels of the points;
 *   normal    coordinates of the inner normal to the cell,
 *             of dimension n. */

int intcelcon_make_subdivision ( void );
/*
 * DESCRIPTION :
 *   Computes the mixed cells for the lifted points stored in the
 *   container, with respect to the defined type of mixture. */

int intcelcon_clear_mixed_cell_configuration ( void );
/*
 * DESCRIPTION :
 *   Clears the cells container. */

#endif
