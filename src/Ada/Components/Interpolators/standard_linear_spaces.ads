with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;

package Standard_Linear_Spaces is

-- DESCRIPTION :
--   This package provides routines to determine the dimension and to set up
--   the equations of linear spaces for standard point configurations.

  procedure Rank ( nbvecs,n : in integer32;
                   pts : in Matrix; tol : in double_float;
                   trivec : out Matrix; rnk : out natural32 );

  -- DESCRIPTION :
  --   Determines the dimension of the vector space spanned by
  --   the points in the matrix pts.

  -- REQUIRED : trivec'range(1) contains 1..nbvecs, pts'range(1) contains
  --   1..nbvecs+1 and trivec'range(2) = pts'first(2)..pts'last(2)-1 = 1..n.

  -- ON ENTRY :
  --   nbvecs      number of vectors, #points is therefore nbvecs+1;
  --   n           dimension of the points;
  --   pts         matrix with in its rows the points;
  --   tol         tolerance to decide whether element is zero or not.

  -- ON RETURN :
  --   trivec      vector configuration in upper triangular form;
  --   rnk         rank of the matrix trivec, is dimension of space.

  function Pivots ( nbvecs,n : integer32;
                    trivec : Matrix; rnk : natural32; tol : double_float )
                  return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the column indices of the nonzero elements in the upper
  --   triangular matrix used to determine the rank.
  --   The matrix "trivec" and "rnk" are output of the procedure Rank.

  function Kernel ( pts,trivec : Matrix; rnk : natural32;
                    pivots : Standard_Integer_Vectors.Vector;
                    tol : double_float ) return VecVec;

  -- DESCRIPTION :
  --   Returns n - rnk hyperplane equations to describe the space that
  --   contains the points in pts, where n = pts'last(1).

  -- REQUIRED : trivec'range(1) = pts'range(1)
  --        and trivec'range(2) = pts'first(2)..pts'last(2)-1.

  -- ON ENTRY :
  --   pts         point configuration in the rows of the matrix;
  --   trivec      vector configuration in upper triangular form;
  --   rnk         rank of trivec, dimension of the points configuration;
  --   pivots      equals Pivots(trivec,rnk,tol);
  --   tol         tolerance to decide whether element is zero or not.

  -- ON RETURN :
  --   The space is the intersection of n - rnk hyperplanes, with
  --   coefficients hyps(i), for i in hyps'range, as follows :
  --     hyps(i)(0) + hyps(i)(1)*x(1) + .. + hyps(i)(n)*x(n) = 0,
  --   with hyps = Kernel(pts,trivec,rnk).

end Standard_Linear_Spaces;
