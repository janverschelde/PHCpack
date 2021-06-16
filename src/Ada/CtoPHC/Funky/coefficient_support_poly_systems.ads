with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;

package Coefficient_Support_Poly_Systems is

-- DESCRIPTION :
--   This package contains functions to build a coefficient-support
--   representation of a polynomial system with complex coefficients.
--   The components of the coefficient and support vectors are respectively
--   C doubles and C integers.  This low level representation provides a
--   protocol to pass multivariate complex polynomials from C to Ada and
--   from Ada to C.

  function Monomial_Count ( p : Poly_Sys ) return C_Integer_Array;

  -- DESCRIPTION :
  --   Returns a C integer array of range 0..p'length with the number of
  --   monomials in every polynomial of p.

  function Sum ( a : C_Integer_Array ) return integer32;

  -- DESCRIPTION :
  --   Returns the sum of the integers in the array a.

  function Support ( n,m : natural32; moncnt : C_Integer_Array;
                     p : Poly_Sys ) return C_Integer_Array;

  -- DESCRIPTION :
  --   Returns the support of the polynomial p as a one dimensional array
  --   of C integers, of range 0..n*m-1, where n equals the number of 
  --   unknowns in every polynomial of p and where m equals the total number 
  --   of terms in p, or the sum of all entries in moncnt.  
  --   The array on return is a sequence of blocks of n numbers,
  --   where each block encodes an exponent vector of a monomial of p.

  function Coefficients ( m : natural32; moncnt : C_Integer_Array;
                          p : Poly_Sys ) return C_Double_Array;

  -- DESCRIPTION :
  --   Returns the coefficients of the system p as a one dimensional
  --   array of C doubles, of range 0..2*m-1, where m equals the total
  --   number of terms in p.  The doubles in the array on return represent
  --   consecutively the real and imaginary parts of each coefficient.

  function Create ( n : natural32; m : C_Integer_Array;
                    c : C_Double_Array; s : C_Integer_Array ) return Poly_Sys;

  -- DESCRIPTION :
  --   Creates a polynomial system in n variables with coefficients in c and
  --   the support in s.  Note that n does not have to match with m'length.

  -- REQUIRED :
  --   If m(i) is the number of monomials of the i-th polynomial on return,
  --   then c'length is 2*m(i) and s'length is a multiple of n*m(i).

  function Concat ( n : natural32; m : C_Integer_Array;
                    c : C_Double_Array; s : C_Integer_Array )
                  return C_Double_Array;

  -- DESCRIPTION :
  --   Concatenates all information to create a polynomial system
  --   in coefficient support representation into one double array.
  --   Applying the next four function to the result of Concat allows
  --   to reconstruct the polynomial system.

  function Dimension ( x : C_Double_Array ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension from the concatenated coefficient
  --   support representation of a polynomial system.

  function Monomial_Count ( x : C_Double_Array ) return C_Integer_Array;

  -- DESCRIPTION :
  --   Returns the monomial count from the concatenated coefficient
  --   support representation of a polynomial system.

  function Support ( x : C_Double_Array ) return C_Integer_Array;

  -- DESCRIPTION :
  --   Returns the support from the concatenated coefficient
  --   support representation of a polynomial system.

  function Coefficients ( x : C_Double_Array ) return C_Double_Array;

  -- DESCRIPTION :
  --   Returns the coefficients from the concatenated coefficient
  --   support representation of a polynomial system.

end Coefficient_Support_Poly_Systems;
