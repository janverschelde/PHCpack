with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;

package Standard_Random_Laurentials is

-- DESCRIPTION :
--   Basic facility to generate random Laurent polynomials.

  function Create_Laurent_Polynomial
             ( c : Standard_Complex_Vectors.Vector;
               e : Standard_Integer_Matrices.Matrix ) return Poly;

  -- DESCRIPTION :
  --   Returns a Laurent polynomial with coefficients in c
  --   and exponents in e.

  -- REQUIRED : c'range = e'range(2).

  function Random_Laurent_Polynomial 
             ( n,m : natural32; lower,upper : integer32 ) return Poly;

  -- DESCRIPTION :
  --   Returns a random n-variate Laurent polynomial with m monomials,
  --   with exponents uniformly generated between lower and upper.
  --   The coefficients are randomly generated on the complex unit circle.

end Standard_Random_Laurentials;
