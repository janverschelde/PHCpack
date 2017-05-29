with QuadDobl_Complex_Vectors;          use QuadDobl_Complex_Vectors;
with Symbol_Table;                      use Symbol_Table;
with QuadDobl_Complex_Polynomials;      use QuadDobl_Complex_Polynomials;
with QuadDobl_Pade_Approximants;        use QuadDobl_Pade_Approximants;

package QuadDobl_Pade_Approximants_io is

-- DESCRIPTION :
--   Provides operations to represents Pade approximants,
--   defined by coefficient vectors of numerator and denominator
--   in double double precision.

  function to_Poly ( c : Vector ) return Poly;

  -- DESCRIPTION :
  --   Given the coefficients in the vector c,
  --   return the polynomial representation of the polynomial
  --   c(0) + c(1)*x + c(2)*x^2 + .. + c(d)*x^d, where d = c'last.

  -- REQUIRED : c'first = 0.

  function Write ( c : Vector; s : Symbol ) return string;

  -- DESCRIPTION :
  --   Writes the string representation of the polynomial with
  --   coefficients in c using the symbol s for the variable.

  function Write ( p : Pade; s : Symbol ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the Pade approximant,
  --   using the symbol s for the variable.

end QuadDobl_Pade_Approximants_io;
