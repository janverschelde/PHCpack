with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;          use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;          use QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;         use QuadDobl_Complex_Matrices;

package QuadDobl_Plane_Operations is

-- DESCRIPTION :
--   This package offers some useful operations on planes,
--   represented by generators or equations in double double precision.

-- NOTE : this is a left over package, 
--   see standard_plane_representations and standard_point_coordinates;
--   The Intersect operation was moved over from standard_cascading_planes.

  procedure Random_Affine_Plane ( n,k : in integer32;
                                  b : out Vector; v : out VecVec );

  -- DESCRIPTION :
  --   Generates a random affine plane of dimension k in n-space.
  --   The directions in v are not orthogonal with respect to each other.

  function Random_Point ( b : Vector; v : VecVec ) return Vector;

  -- DESCRIPTION :
  --   Returns a random point as a sum of the basis b and a random
  --   linear combination of the direction vectors in v.

  function Evaluate ( h,x : Vector ) return Complex_Number;
  function Evaluate ( h : VecVec; x : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the evaluation of the hyperplane(s) at the point x.

  -- REQUIRED : h'range = 0..x'last; x'range = 1..h'last.

  function Orthogonalize ( v : Array_of_VecVecs ) return Array_of_VecVecs;

  -- DESCRIPTION :
  --   Returns an array of orthogonal bases for the given spaces in v.

  function In_Span ( v : VecVec; x : QuadDobl_Complex_Vectors.Vector;
                     tol : quad_double ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the projection of x onto the space spanned by v
  --   has norm less than or equal to the given tolerance tol.

  -- REQUIRED : the vectors in v form an orthonormal basis for their span.

  procedure Affine_Orthonormal_Basis
              ( n,k : in integer32; slices : in VecVec;
                b : out Vector; v,w : out VecVec );

  -- DESCRIPTION :
  --   Computes an orthonormal basis for the affine space with given slices.
  --   Note that the slices are given in the embedded space.

  -- ON ENTRY :
  --   n        dimension of the ambient space, before the embedding;
  --   k        co-dimension of the linear space;
  --   slices   coefficients of the linear slices in their embedding.

  -- ON RETURN :
  --   b        a basis point on the affine space;
  --   v        linearly independent set of vectors is a basis of the
  --            affine space minus the basis point;
  --   w        orthormal basis for the space minus the basis point.

  function Truncate ( v : in VecVec; n : in integer32 ) return VecVec;

  -- DESCRIPTION :
  --   Returns a vector of vectors, restricted to the first n entries.

  function Truncate ( v : in Array_of_VecVecs; n : in integer32 )
                    return Array_of_VecVecs;

  -- DESCRIPTION :
  --   Returns the array of vectors of vectors, where every vector is
  --   restricted to the first n entries.

  procedure Evaluate ( file : in file_type;
                       equ : in Matrix; v : in Vector;
                       res : out quad_double );
  procedure Evaluate ( equ : in Matrix; v : in Vector;
                       res : out quad_double );

  -- DESCRIPTION :
  --   Evaluates the vector v at the equations in equ.
  --   The 1-norm of the evaluation is in res.
  --   Intermediate output is written to file.

  procedure Evaluate ( file : in file_type;
                       p : in Matrix; g : in Matrix;
                       res : out quad_double );
  procedure Evaluate ( p : in Matrix; g : in Matrix;
                       res : out quad_double );

  -- DESCRIPTION :
  --   Verifies whether the generators in g satisfy the equations in p.
  --   The 1-norm of the evaluation is in res.
  --   Intermediate output is written to file.

  procedure Intersect ( e1,e2 : in Matrix; p1,p2 : in out Matrix );
  procedure Intersect ( file : in file_type;
                        e1,e2 : in Matrix; p1,p2 : in out Matrix );

  -- DESCRIPTION :
  --   Computes the intersection of two hyperplanes.

  -- REQUIRED : despite the general setup of the parameters,
  --   the code works only for hyperplanes.

  -- ON ENTRY :
  --   file      for writing diagnostics to, if file is provided;
  --   e1,e2     equations for the two planes;
  --   p1,p2     orthogonal representations of the two planes,
  --             with common offset vector.

  -- ON RETURN :
  --   p1,p2     common basis, except for one complementary vector,
  --             not in the intersection of both planes.

end QuadDobl_Plane_Operations;
