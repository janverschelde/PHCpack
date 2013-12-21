with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Standard_Integer32_Vertices is

-- DESCRIPTION :
--   This function offers some functions for computing the 
--   vertices that span the polytope, given its set of points.
--   Based on the functions in this package it is possible to
--   determine the dimension of conv(L), without the actual computation
--   of the convex hull.

  function Is_In_Hull ( point : Standard_Integer_Vectors.Vector;
                        L : List ) return boolean;
  function Is_In_Hull ( point : Standard_Floating_Vectors.Vector;
                        L : List ) return boolean;

  -- DESCRIPTION :
  --   This function determines whether the point belongs to convex
  --   hull of the points in the given list.

  function Vertex_Points ( L : List ) return List;

  -- DESCRIPTION :
  --   The vertex points of conv(L) are returned.

  function Extremal_Points
              ( L : List; v : Standard_Integer_Vectors.Link_to_Vector )
              return List;

  -- DESCRIPTION :
  --   Returns the vertex points of the list L w.r.t. the direction v and -v.

  function Extremal_Points ( k,n : integer32; exL,L : List ) return List;
  function Max_Extremal_Points ( k,n : integer32; exL,L : List ) return List;

  -- DESCRIPTION :
  --   There are already k linearly independent vertex points found,
  --   given in exL.  This routine tries to take one more linearly
  --   independent vertex point out of the list L.
  --   The first function stops when a degeneracy is discovered,
  --   while the second one still proceeds and returns one point
  --   more, when possible.

  function Extremal_Points ( n : integer32; L : List ) return List;

  -- DESCRIPTION :
  --   Searches in the list L for linearly independent vertex points.
  --   If dim(conv(L)) = n then n+1 points will be returned,
  --   otherwise the number of points in the resulting list will be
  --   less but not necessarily equal to dim(conv(L))+1.
  --   The advantage of this routine lies in the fact that a degeneracy,
  --   i.e., dim(conv(L)) < n, is detected as soon as possible.

  function Max_Extremal_Points ( n : integer32; L : list ) return List;

  -- DESCRIPTION :
  --   Does the same as the function above, except for the fact that
  --   the number of points in the resulting lists equals dim(conv(L))+1.
  --   The points in the resulting list can be regarded as a basis,
  --   i.e., origin and as many linearly independent points as dim(conv(L)),
  --   for the points in L.

end Standard_Integer32_Vertices;
