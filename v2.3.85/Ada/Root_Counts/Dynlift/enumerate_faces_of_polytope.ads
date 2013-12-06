with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Integer32_Simplices;       use Standard_Integer32_Simplices;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;

package Enumerate_Faces_of_Polytope is

-- DESCRIPTION :
--   The routines in this package enable the enumeration of all k-faces
--   of the lower hull of a polytope, based on a regular triangulation.

  generic
    with procedure Process ( face : in VecVec );
  procedure Enumerate_Faces_in_List
                ( L : in List; point : in Link_to_Vector; k : in integer32 );
  generic
    with procedure Process ( face : in VecVec );
  procedure Enumerate_Lower_Faces_in_List
                ( L : in List; point : in Link_to_Vector; k : in integer32 );

  -- DESCRIPTION :
  --   Enumerates all k-dimensional faces spanned by the list and
  --   the given point.  The _Lower_ indicates that only the faces of the
  --   lower hull are enumerated.

  -- ON ENTRY :
  --   L          a list of points;
  --   point      vector not in the list L;
  --   k          dimension of the faces.

  -- ON RETURN :
  --   face       a (k+1)-dimensional vector: face(0..k), face(0) = point.

  generic
    with procedure Process ( face : in VecVec );
  procedure Enumerate_Faces_in_Simplex
                ( s : in Simplex; point : in Link_to_Vector; k : in integer32 );
--  generic
--    with procedure Process ( face : in VecVec );
--  procedure Enumerate_Lower_Faces_in_Simplex
--                ( s : in Simplex; point : in Link_to_Vector; k : in natural );

  -- DESCRIPTION :
  --   Enumerates all k-dimensional faces in the simplex s that contain
  --   the given point.  The _Lower_ indicates that only the faces of the
  --   lower hull are enumerated.

  -- ON ENTRY :
  --   s          a simplex;
  --   point      vector in the simplex;
  --   k          dimension of the faces.

  -- ON RETURN :
  --   face       a (k+1)-dimensional vector: face(0..k), face(0) = point.

  generic
    with procedure Process ( face : in VecVec );
  procedure Enumerate_Faces_in_Triangulation
                ( t : in Triangulation; point : in Link_to_Vector;
                  k : in integer32 );
--  generic
--    with procedure Process ( face : in VecVec );
--  procedure Enumerate_Lower_Faces_in_Triangulation
--                ( t : in Triangulation; point : in Link_to_Vector;
--                  k : in natural );

  -- DESCRIPTION :
  --   Enumerates all k-dimensional faces in the simplices of the triangulation
  --   that contain the given point.  The _Lower_ indicates that only the faces
  --   of the lower hull are enumerated.

  -- ON ENTRY :
  --   t          a regular triangulation;
  --   point      vector in each simplex of t;
  --   k          dimension of the faces.

  -- ON RETURN :
  --   face       a (k+1)-dimensional vector: face(0..k), face(0) = point.

end Enumerate_Faces_of_Polytope;
