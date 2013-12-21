with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;

package Standard_Moving_Planes is

-- DESCRIPTION :
--   Functions in this package evaluate a homotopy between linear spaces
--   given in parametric representations.

  function Random_Plane ( n,k : integer32 ) return Matrix;

  -- DESCRIPTION :
  --   Returns a matrix representing a k-plane in n-space.

  function One_Random_Direction ( m : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns the same k-plane as m, except for the last direction
  --   which is chosen at random.

  function Rotate ( A : Matrix; theta : double_float;
                    i1,i2 : integer32 ) return Matrix;

  -- DESCRIPTION :
  --   Rotates coordinates i1 and i2 of the plane A around theta.

  function Rotating_Plane ( A : Matrix; i1,i2 : integer32;
                            t : Complex_Number ) return Matrix;

  -- DESCRIPTION :
  --   Rotates coordinates i1 and i2 of A around the angle 2*PI*t.

  function Moving_Plane ( A,B : Matrix; t : Complex_Number ) return Matrix;

  -- DESCRIPTION :
  --   Returns the matrix (1-t)*A + t*B.

  function Moving_Plane ( A,B : Matrix; gamma,t : Complex_Number )
                        return Matrix;

  -- DESCRIPTION :
  --   Returns the matrix [(1-t) + gamma*t*(1-t)]*A + t*B.

  function Moving_Directions
              ( A,B : Matrix; t : Complex_Number; ortho : boolean )
              return Matrix;

  -- DESCRIPTION :
  --   Returns the combination (1-t)*A + t*B of the directions,
  --   that is: in the columns 1 to A'last(2) = B'last(2).

  -- REQUIRED : A'last(2) = B'last(2) = k.

  -- ON ENTRY :
  --   A        start k-plane;
  --   B        target k-plane;
  --   t        current value of the continuation parameter;
  --   ortho    if true, then the columns of the plane on return
  --            define an orthogonal matrix.

  procedure Normalize ( A : in out Matrix; col : in integer32 );
  procedure Normalize ( A : in out Matrix; col : in integer32;
                        nrm : out double_float );

  -- DESCRIPTION :
  --   Normalizes a column in A, defined by the index col.
  --   This procedure is needed for when only one column changes
  --   during the deformation.  The optional output argument
  --   returns the norm of the column before normalization.

end Standard_Moving_Planes;
