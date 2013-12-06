with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;

package Standard_Rescaling_Coordinates is

-- DESCRIPTION :
--   For improved numerical stability of intrinsic coordinates in
--   witness sets, we use local coordinates.  The operations here
--   help to implement a rescaling algorithm.

  function Linear_Offset_Shift
             ( a,b : Vector; t : Complex_Number ) return Vector;

  -- DESCRIPTION :
  --   Returns (1-t)*a + t*b.

  function Complement_of_Projection ( p : Matrix; v : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the complement of the projection of the vector v
  --   onto the vectors in the columns 1..p'last(2) of p.

  -- REQUIRED : the vectors in 1..p'last(2) form an orthonormal basis.

  function Distance ( p : Matrix; x : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the distance of the point x to the plane p.

  -- REQUIRED : the vectors in 1..p'last(2) form an orthonormal basis.

  procedure Check_Orthonormality
              ( p : in Matrix; tol : in double_float; 
                fail : out boolean );

  -- DESCRIPTION :
  --   Checks the orthonormality of the vectors that span the plane.
  --   The plane p is orthonormal if
  --   (1) all columns in p have normal equal to one (modulo tol); and
  --   (2) the conjugated inner product of mutually distinct columns
  --       is less than the given tolerance.
  --   If the orthonormality conditions are satified, then fail is false,
  --   otherwise the plane p fails the orthonomality test.

end Standard_Rescaling_Coordinates;
