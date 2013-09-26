with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;

package Standard_Cascading_Planes is

-- DESCRIPTION :
--   Provides utilities to build a sequence of linear spaces,
--   for an intrinsic diagonal homotopy method to intersect
--   positive dimensional varieties.

  function Project ( x : Matrix; h : integer32 ) return Matrix;

  -- DESCRIPTION :
  --   Copies the first h rows of x and sets the other rows to zero.
  --   The matrix on return has the same dimensions as x.

  function Project ( x : Vector; h : integer32 ) return Vector;

  -- DESCRIPTION :
  --   Copies the first h rows of x and sets the other rows to zero.
  --   The vector on return has the same dimensions as x.

  function Double_to_Diagonal ( A : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Doubles the number of columns of A, according to the
  --   equations on the diagonal.

  function Compute_Offset ( C : Matrix; d : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the solution w to the system C*w + d = 0.

  procedure Shift_Offset ( p : in out Matrix; b : in Vector );

  -- DESCRIPTION :
  --   Updates the generators of the plane p with a new offset b.

  function Target_Space ( n,n2,apb,b : integer32; dA,BB,CC : Matrix;
                          d : Vector ) return Matrix;

  -- DESCRIPTION :
  --   Generates the equations for the target linear space.

  -- ON ENTRY :
  --   n         dimension of the ambient space before embedding;
  --   n2        equals 2*n;
  --   apb       is a + b, sum of dimension of the two components;
  --   b         smallest dimension of the two components;
  --   dA        matrix of a+b rows and n2 columns, randomized diagonal;
  --   BB        matrix of a+b rows and n2 columns, is a multiplier;
  --   CC,d      defines the system CC*w + d = 0.

  -- ON RETURN :
  --   matrix of a+b rows and columns of range 0..n2.

  function Start_Space ( g1,g2 : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns an intrinsic representation of the plane at the start,
  --   from the generators of the planes cutting the witness sets of
  --   the two components.

end Standard_Cascading_Planes;
