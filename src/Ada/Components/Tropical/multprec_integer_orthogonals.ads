--with Multprec_Integer64_Numbers;         use Multprec_Integer64_Numbers;
--with Multprec_Integer64_Vectors;         use Multprec_Integer64_Vectors;
--with Multprec_Integer64_Matrices;        use Multprec_Integer64_Matrices;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Vectors;           use Multprec_Integer_Vectors;
with Multprec_Integer_Matrices;          use Multprec_Integer_Matrices;

package Multprec_Integer_Orthogonals is

-- DESCRIPTION :
--   Applies the Gram-Schmidt orthogonalization procedure
--   to a set of multiprecision integer vectors.

  function gcd ( A : Matrix; k : integer32 ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the greatest common divisor of the elements in
  --   the k-th column of A.

  procedure Normalize ( A : in out Matrix; k : in integer32 );

  -- DESCRIPTION :
  --   Divides all entries in the k-th column of A
  --   by their greatest common divisor.

  function Orthogonalize ( A : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns an orthogonal basis of the space spanned by 
  --   the vectors in the columns of A.

  function Complement ( A : Matrix ) return Vector;

  -- DESCRIPTION :
  --   Returns one vector perpendicular to the columns of A.

end Multprec_Integer_Orthogonals;
