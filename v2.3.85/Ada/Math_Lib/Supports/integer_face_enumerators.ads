with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;

package Integer_Face_Enumerators is

-- DESCRIPTION :
--   This package contains procedures which allow the efficient
--   enumeration of all vertices, edges and k-faces of a polytope.
--   The polytope is spanned by a given vector of points in a general
--   dimensional space.
--   In addition, a procedure has been added for enumerating all
--   facets of a certain type of the sum of a tuple of polytopes.

  generic

    with procedure Process ( i : in integer32; cont : out boolean );

  procedure Enumerate_Vertices ( pts : in VecVec );

  -- DESCRIPTION :
  --   The vector pts contains the points which span the polytope.
  --   The candidate vertices are enumerated in order of occurrence in pts.
  --   Each time a vertex has been found, the procedure Process is invoked,
  --   which returns then the entry of the vertex in pts.
  --   When cont is set to false, the enumeration stops.

  generic

    with procedure Process ( i,j : in integer32; cont : out boolean );

  procedure Enumerate_Edges ( pts : in VecVec );

  -- DESCRIPTION :
  --   The vector pts contains the points which span the polytope.
  --   The candidate edges are enumerated in lexicographic order.
  --   Each time an edge has been found, the procedure Process is
  --   invoked, which returns the entries of the edge in pts.
  --   When cont is set to be false, the enumeration stops.

  generic

    with procedure Process ( i,j : in integer32; cont : out boolean );

  procedure Enumerate_Lower_Edges ( pts : in VecVec );

  -- DESCRIPTION :
  --   Enumerates only the edges on the lower hull of the polytope.

  generic
    
    with procedure Process ( face : in Vector; cont : out boolean );

    -- REQUIRED :
    --   The range of face must be 0..k.

  procedure Enumerate_Faces ( k : in integer32; pts : in VecVec );

  -- DESCRIPTION :
  --   The vector pts contains the points which span the polytope.
  --   A k-face is a face spanned by k affine linearly independent points.
  --   The candidate k-faces are enumerated in lexicographic order.
  --   Each time an edge has been found, the procedure Process is
  --   invoked, which returns the entries of the k-face in pts.
  --   When cont is set to be false, the enumeration stops.

  generic
    
    with procedure Process ( face : in Vector; cont : out boolean );

    -- REQUIRED :
    --   The range of face must be 0..k.

  procedure Enumerate_Lower_Faces ( k : in integer32; pts : in VecVec );

  -- DESCRIPTION :
  --   Enumerates only the k-faces on the lower hull of the polytope.

  generic

    with procedure Process ( faces : in VecVec; cont : out boolean );

    -- DESCRIPTION :
    --   The parameter faces has the following meaning:
    --   faces(i) contains the entries of points in pts which span
    --   the k-face of the ith polytope, with k = typ(i).
    --   When cont is set to false, the enumeration stops.

  procedure Enumerate_Faces_of_Sum
               ( ind,typ : in vector; k : in integer32; pts : in VecVec );

  -- DESCRIPTION :
  --   Enumerates the points which span faces of a certain type
  --   of the polytope in lexicographic order.

  -- ON ENTRY :
  --   ind       ind(i) indicates where the points of the ith polytope
  --             in the vector pts begin;
  --   typ       typ(i) contains the dimension k of the k-face of the ith
  --             polytope which spans the facet of the sum;
  --   k         sum of the entries in typ, must be less than
  --             or equal to the dimension of the space;
  --   pts       contains all points of the polytopes.

  -- ON RETURN :
  --   Each time a face of the sum has been found, Process is invoked.

end Integer_Face_Enumerators;
