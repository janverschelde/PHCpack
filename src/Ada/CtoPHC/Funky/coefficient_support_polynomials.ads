with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;

package Coefficient_Support_Polynomials is

-- DESCRIPTION :
--   This package contains functions to convert a polynomial in several
--   variables with complex coefficients into two vectors: its coefficients
--   and its support.  The components of these vectors are respectively
--   C doubles and C integers.  This low level representation provides a
--   protocol to pass multivariate complex polynomials from C to Ada and
--   from Ada to C.

  function Support ( p : Poly ) return C_Integer_Array;

  -- DESCRIPTION :
  --   Returns the support of the polynomial p as a one dimensional array
  --   of C integers, of range 0..m-1, where m equals the number of terms 
  --   times the number of unknowns of p.  If n is the number of unknowns
  --   of p, the the array on return can be viewed as a sequence of blocks
  --   of n numbers, where each block encodes the exponent vector of a
  --   monomial of p.

  function Coefficients ( p : Poly ) return C_Double_Array;

  -- DESCRIPTION :
  --   Returns the coefficients of the polynomial p as a one dimensional
  --   array of C doubles, of range 0..m-1, where m equals twice the
  --   number of terms in p.  The doubles in the array on return represent
  --   consecutively the real and imaginary parts of each coefficient.

  function Create ( n : natural32; c : C_Double_Array; s : C_Integer_Array ) 
                  return Poly;

  -- DESCRIPTION :
  --   Creates a polynomial in n variables with coefficients in c and
  --   the support in s.

  -- REQUIRED :
  --   If m is the number of monomials of the created polynomial,
  --   then c'length is 2*m and s'length is a multiple of n*m.

end Coefficient_Support_Polynomials;
