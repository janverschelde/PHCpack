with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Octo_Double_Numbers;                 use Octo_Double_Numbers;
with OctoDobl_Complex_Numbers;            use OctoDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with OctoDobl_Complex_Series;             use OctoDobl_Complex_Series;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_Complex_Series_Matrices;

package OctoDobl_Series_Linear_Solvers is

-- DESCRIPTION :
--   A square linear system of series is solved by LU factorization,
--   followed by two backsubstitutions.

  function cabs ( c : Complex_Number ) return octo_double;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of real
  --   and imaginary parts of the complex number c.

  function cabs ( s : Series ) return octo_double;

  -- DESCRIPTION :
  --   Returns the cabs of the constant coefficient of the series s.

  procedure LUfac ( A : in out OctoDobl_Complex_Series_Matrices.Matrix;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 );

  -- DESCRIPTION :
  --   LUfac factors a complex matrix by gaussian elimination.

  -- ON ENTRY :
  --   A       series matrix(1..n,1..n) to be factored;
  --   n       the dimension of the matrix A.

  -- ON RETURN :
  --   A       an upper triangular matrix and the multipliers
  --           which were used to obtain it.
  --           The factorization can be written A = L*U where
  --           L is a product of permutation and unit lower
  --           triangular matrices and U is upper triangular;
  --   ipvt    an integer vector of pivot indices;
  --   info    = 0  normal value;
  --           = k  if u(k,k) = 0.0;
  --                This is not an error for this routine,
  --                but it does indicate that LUsolve will
  --                divide by zero if called.  Use rcond in
  --                lufco on the leading constants coefficients
  --                for a reliable indication of singularity.

  procedure LUsolve ( A : in OctoDobl_Complex_Series_Matrices.Matrix;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out OctoDobl_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   LUsolve solves the complex system A*x = b using the factors
  --   computed by LUfac.

  -- ON ENTRY :
  --   A       series matrix(1..n,1..n), the output from LUfac;
  --   n       the dimension of the matrix A;
  --   ipvt    the pivot vector from LUfac;
  --   b       the right hand side vector.

  -- ON RETURN :
  --   b       the solution vector x.

end OctoDobl_Series_Linear_Solvers;
