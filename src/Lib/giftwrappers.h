/* The file giftwrappers.h contains prototypes to the gift wrapping methods,
 * wrapped through the Ada code use_giftwrap of PHCpack.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __GIFTWRAPPERS_H__
#define __GIFTWRAPPERS_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int convex_hull_2d ( int nc_pts, char *pts, int *nc_hull, char *hull );
/*
 * DESCRIPTION :
 *   Computes the convex hull of a planar point configuration.
 *
 * REQUIRED :
 *   All coordinates of the points have integer values.
 *
 * ON ENTRY :
 *   nc_pts   number of characters in the string representation of
 *            a point configuration, as a list of tuples;
 *   pts      string representation of a Python-style list of tuples,
 *            e.g.: "[(7, 3), (-2, 1), (2, -5), (-3, -4), (-3, 7)]"
 *
 * ON RETURN :
 *   nc_hull  number of characters in the string representation of
 *            the convex hull, as a tuple of two lists of tuples;
 *   hull     string representation of two tuples, the first is a
 *            list of vertex points, the second element of the tuple
 *            is a list of inner normals. */  

int convex_hull ( int nc_pts, char *pts );
/*
 * DESCRIPTION :
 *   Computes the convex hull of a point configuration in 3-space or 4-space.
 *   The container is initialized with the list of facets.
 *
 * REQUIRED :
 *   All coordinates of the points have integer values.
 *
 * ON ENTRY :
 *   nc_pts   number of characters in the string representation of
 *            a point configuration, as a list of tuples;
 *   pts      string representation of a Python-style list of tuples. */

int number_of_facets ( int dim, int *nbr );
/*
 * DESCRIPTION :
 *   Returns the number of facets of a convex hull in the space
 *   of dimension dim.
 *
 * REQUIRED :
 *   The function convex_hull in the correct dimension has been called.
 *
 * ON ENTRY :
 *   dim      dimension of the ambient space.
 *
 * ON RETURN :
 *   nbr      number of facets of the convex hull. */

int retrieve_facet ( int dim, int fcn, int *nc_rep, char *fctrep );
/*
 * DESCRIPTION :
 *   Returns the string representation of a facet given the dimension
 *   and its facet number.
 *
 * REQUIRED :
 *   The function convex_hull in the correct dimension has been called
 *   and fcn is in the range 0..nbr-1, where number_of_facets(dim,&nbr).
 *
 * ON ENTRY :
 *   dim      dimension of the ambient space;
 *   fcn      facet number.
 *
 * ON RETURN :
 *   nc_rep   number of characters in the string representation of the facet;
 *   fctrep   string representation of the facet with number fcn, as a tuple
 *            of four elements, at respective positions in the tuple:
 *            0. the value of the support function at the facet,
 *            1. the inner normal to the facet, which jointly with item 0
 *               defines the hyperplane supporting the facet;
 *            2. indices of the points that span the facet,
 *               the index k will refer to the k-th point represented by
 *               the argument pts of the function convex_hull above,
 *               note that the index count starts at one instead of zero;
 *            3. number of adjacent facets, with the facet count starting
 *               at zero. */

int clear_3d_facets ( void );
/*
 * DESCRIPTION :
 *   Clears the list of facets in 3-space. */

int clear_4d_facets ( void );
/*
 * DESCRIPTION :
 *   Clears the list of facets in 4-space. */

int support_size ( int idx );
/*
 * DESCRIPTION :
 *   Returns the number of characters in the string representation of
 *   the support of the Laurent polynomial with index idx in the container.
 *
 * REQUIRED :
 *   The Laurent systems container must be initialized.
 *
 * ON ENTRY :
 *   idx     number between 1 and the number of polynomials stored. */

int support_string ( int size, char *supp );
/*
 * DESCRIPTION :
 *   Returns the string representation of the support a Laurent polynomial.
 *
 * REQUIRED :
 *   The Laurent systems container must be initialized
 *   and supp must be allocated for as many characters as size.
 *
 * ON ENTRY :
 *   size     the number of characters in the string representation
 *            of the support of a Laurent polynomial.
 *
 * ON RETURN :
 *   support  the string representation of the support of a Laurent
 *            polynomial. */

int clear_support_string ( void );
/*
 * DESCRIPTION :
 *   Deallocates the string representation of the support set
 *   as stored internally. */

int initial_form ( int dim, int nbc, char *normal );
/*
 * DESCRIPTION :
 *   Replaces the system in the Laurent system container by its initial form,
 *   defined by the inner normal with coordinates in the parameter normal,
 *   where the number of coordinates equals dim.
 *   The normal is given as a Python tuple, stored as a string. */

#endif
