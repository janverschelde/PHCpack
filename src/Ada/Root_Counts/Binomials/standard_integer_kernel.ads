with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;          use Standard_Integer_Vectors;
with Standard_Integer_Matrices;         use Standard_Integer_Matrices;

package Standard_Integer_Kernel is

-- DESCRIPTION :
--   This package offers tools to compute the kernel of an integer matrix,
--   with standard 32-bit arithmetic.

  function Rank_Upper ( U : Matrix ) return integer32;

  -- DESCRIPTION :
  --   Given an upper triangular matrix, returns its rank.

  procedure Pivots_in_Upper ( U : in Matrix; rank : out integer32;
                              p : out Vector );

  -- DESCRIPTION :
  --   Given an upper triangular matrix U, computes its rank and 
  --   returns the first nonzero element on each row in p.

  -- ON ENTRY :
  --   U        an upper triangular matrix.

  -- ON RETURN :
  --   rank     the rank of U;
  --   p        p(i) is the index of the first nonzero on row i of U,
  --            p(i) equals 0 if all elements on the i-th row of U are zero.

  -- REQUIRED : p'first = U'first(1) and p'last >= U'last(1).

  procedure Normalize_Sign ( v : in out Vector );

  -- DESCRIPTION :
  --   On return the orientation of the vector v might be changed 
  --   so the sign of its first nonzero element is positive.

  function Is_In ( v : Vector; x : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if x belongs to v, false otherwise.

  function Complement ( n,r : integer32; p : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector of size n-r with the complementary pivots to p.

  -- REQUIRED :
  --   The range of p is 1..r and p contains r distinct nonzero
  --   elements in the range 1..n.
 
  procedure Upper_Kernel ( U : in Matrix; n,r : in integer32;
                           p : in Vector; V : out Link_to_Matrix );

  -- DESCRIPTION :
  --   Computes a basis for the kernel of an upper triangular matrix U.
 
  -- ON ENTRY :
  --   U          an upper triangular matrix of rank r;
  --   n          ambient dimension, number of columns of U;
  --   r          rank of the matrix U;
  --   p          p(i) is the first nonzero element of U on row i.

  -- ON RETURN :
  --   V          a basis for the kernel of U.
 
  -- REQUIRED : r < n, the dimension of the kernel is positive.
 
  procedure Kernel ( A : in Matrix;
                     r : out integer32; V : out Link_to_Matrix );

  -- DESCRIPTION :
  --   Computes a basis for the kernel of a matrix A.

  -- ON ENTRY :
  --   A          a matrix, with n columns.

  -- ON RETURN :
  --   r          the rank of the matrix A;
  --   V          if the rank r is less than A'last(2),
  --              the columns of V span the kernel of A.

end Standard_Integer_Kernel;
