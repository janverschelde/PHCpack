with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_VecVecs;           
with Standard_Floating_VecVecs;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;

package Face_Cardinalities is

-- DESCRIPTION :
--   This package provides a naive facility for counting the number of
--   k-dimensional faces of a polytope.

-- NOTE :
--   The precision is very likely to be insufficient even to count
--   the number of faces.

  function fvector ( pts : Standard_Integer_VecVecs.VecVec ) return Vector;

  function fvector ( pts : Standard_Floating_VecVecs.VecVec ) return Vector;

  -- DESCRIPTION :
  --   Computes the f-vector of a polytope.

  -- ON ENTRY :
  --   pt          a vector of integer points, the polytope is conv(pts).

  -- ON RETURN :
  --   f(-1..n)    the f-vector of conv(pts):
  --                f(-1) = 1, f(0) = #vertices, f(1) = #edges, .. ,
  --                f(n-1) = #facets, f(n) = 1, if conv(pts) is n-dimensional.

  function Face_Labels ( pts : Standard_Integer_VecVecs.VecVec;
                         k : natural32 ) return List;
  function Face_Labels ( pts : Standard_Floating_VecVecs.VecVec;
                         k : natural32 ) return List;

  -- DESCRIPTION :
  --   Returns the list of labels to the points that have been used in
  --   defining the k-faces of the polytope spanned by the points in pts.

  function Face_Labels ( pts : Standard_Integer_VecVecs.VecVec )
                       return Array_of_Lists;
  function Face_Labels ( pts : Standard_Floating_VecVecs.VecVec )
                       return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns the lists of labels to the points that have been used in
  --   defining the k-faces of the polytope spanned by the points in pts.
  --   The array on return has range 0..n, n the dimension of the points.

end Face_Cardinalities;
