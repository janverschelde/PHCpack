with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;          use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;         use QuadDobl_Complex_Matrices;

package QuadDobl_Rescaling_Coordinates is

-- DESCRIPTION :
--   For improved numerical stability of intrinsic coordinates in
--   witness sets, we use local coordinates.  Operations exported
--   by this package help to implement a rescaling algorithm
--   in double double precision.

  function Linear_Offset_Shift
             ( a,b : Vector; t : Complex_Number ) return Vector;

  -- DESCRIPTION :
  --   Returns (1-t)*a + t*b.

  function Complement_of_Projection ( p : Matrix; v : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the complement of the projection of the vector v
  --   onto the vectors in the columns 1..p'last(2) of p.

  -- REQUIRED : the vectors in 1..p'last(2) form an orthonormal basis.

  function Distance ( p : Matrix; x : Vector ) return quad_double;

  -- DESCRIPTION :
  --   Returns the distance of the point x to the plane p.

  -- REQUIRED : the vectors in 1..p'last(2) form an orthonormal basis.

  procedure Check_Orthonormality
              ( p : in Matrix; tol : in quad_double; fail : out boolean );

  -- DESCRIPTION :
  --   Checks the orthonormality of the vectors that span the plane.
  --   The plane p is orthonormal if
  --   (1) all columns in p have normal equal to one (modulo tol); and
  --   (2) the conjugated inner product of mutually distinct columns
  --       is less than the given tolerance.
  --   If the orthonormality conditions are satified, then fail is false,
  --   otherwise the plane p fails the orthonomality test.

end QuadDobl_Rescaling_Coordinates;
