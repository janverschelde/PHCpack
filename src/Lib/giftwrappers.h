/* file giftwrappers.h contains prototypes to the gift wrapping methods,
 * wrapped through the Ada code use_giftwrap of PHCpack */

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
 *   pts      string representation of a Python-style list of tuples.
 *
 * ON RETURN :
 *   nc_hull  number of characters in the string representation of
 *            the convex hull, as a tuple of two lists of tuples;
 *   hull     string representation of two tuples, the first is a
 *            list of vertex points, the second element of the tuple
 *            is a list of inner normals. */  
