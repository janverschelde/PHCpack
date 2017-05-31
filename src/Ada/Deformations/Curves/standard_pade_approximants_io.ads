with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Symbol_Table;                      use Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Pade_Approximants;        use Standard_Pade_Approximants;

package Standard_Pade_Approximants_io is

-- DESCRIPTION :
--   Provides operations to represents Pade approximants,
--   defined by coefficient vectors of numerator and denominator
--   in standard double precision.

  function to_Poly ( c : Vector ) return Poly;

  -- DESCRIPTION :
  --   Given the coefficients in the vector c,
  --   return the polynomial representation of the polynomial
  --   c(0) + c(1)*x + c(2)*x^2 + .. + c(d)*x^d, where d = c'last.

  -- REQUIRED : c'first = 0.

  function to_System ( p : Pade ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns a tuple of two polynomials with the numerator and
  --   the denominator of p, stored as polynomials in one variable.

  function to_System ( p : Pade_Vector ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns a polynomial system of range 1..2*p'length,
  --   where the odd indexed polynomials are the numerators
  --   and the even indexed polynomials are the denominators
  --   of the Pade approximants in the vector.

  function t_symbol return Symbol;

  -- DESCRIPTION :
  --   Returns 't' as the default symbol.

  function Write ( c : Vector ) return string;
  function Write ( c : Vector; s : Symbol ) return string;

  -- DESCRIPTION :
  --   Writes the string representation of the polynomial with
  --   coefficients in c using the symbol s for the variable.
  --   By default, 't' is used as the symbol for the variable..

  function Write ( p : Pade ) return string;
  function Write ( p : Pade; s : Symbol ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of the Pade approximant,
  --   using the symbol s for the variable.
  --   By default, 't' is used as the symbol for the variable.

end Standard_Pade_Approximants_io;
