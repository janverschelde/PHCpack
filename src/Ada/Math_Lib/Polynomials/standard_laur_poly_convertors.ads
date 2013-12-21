with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Standard_Laur_Poly_Convertors is

-- DESCRIPTION :
--   This package contains routines for converting Laurent polynomials to
--   polynomials with natural exponent vectors by shifting when necessary.

  function Negative ( d : Standard_Complex_Laurentials.Degrees ) 
                    return boolean;

  -- DESCRIPTION :
  --   Returns true if there is at least one negative element in d.

  function Is_Genuine_Laurent 
              ( p : Standard_Complex_Laurentials.Poly ) return boolean;
  function Is_Genuine_Laurent ( p : Laur_Sys ) return boolean;

  -- DESCRIPTION :
  --   Returns true if p contains negative exponents.

  function Positive_Laurent_Polynomial
              ( p : Standard_Complex_Laurentials.Poly )
              return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   If the Laurent polynomial has no negative exponents anywhere,
  --   then we call it positive and converting this polynomial to
  --   a polynomial is just a type conversion on the exponents.

  function Positive_Laurent_Polynomial_System ( p : Laur_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Converts a Laurent polynomial system that has no negative exponents
  --   into a polynomial system.

  function Laurent_Polynomial_to_Polynomial
              ( p : Standard_Complex_Laurentials.Poly )
              return Standard_Complex_Polynomials.Poly;

  procedure Laurent_Polynomial_to_Polynomial
              ( L : in Standard_Complex_Laurentials.Poly;
                t : out Standard_Complex_Laurentials.Term;
                p : out Standard_Complex_Polynomials.Poly );

  function Laurent_Polynomial_to_Polynomial
              ( L : Standard_Complex_Laurentials.Poly;
                t : Standard_Complex_Laurentials.Term )
              return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Transforms a Laurent polynomial into an ordinary polynomial
  --   by multiplying by an appropriate monomial t.

  function Laurent_to_Polynomial_System ( p : Laur_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Transforms a Laurent polynomial system into a polynomial system.

end Standard_Laur_Poly_Convertors;
