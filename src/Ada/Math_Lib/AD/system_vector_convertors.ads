with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Monomials;
with Standard_Monomial_Vectors;
with Standard_Polynomial_Vectors;
with DoblDobl_Complex_Monomials;
with DoblDobl_Monomial_Vectors;
with DoblDobl_Polynomial_Vectors;
with QuadDobl_Complex_Monomials;
with QuadDobl_Monomial_Vectors;
with QuadDobl_Polynomial_Vectors;

package System_Vector_Convertors is

-- DESCRIPTION :
--   The Convert functions turn polynomials into the data structures
--   for efficient evaluation and differentiation,
--   in double, double double, and quad double precision.

  function Convert ( t : Standard_Complex_Polynomials.Term )
                   return Standard_Complex_Monomials.Monomial;
  function Convert ( t : DoblDobl_Complex_Polynomials.Term )
                   return DoblDobl_Complex_Monomials.Monomial;
  function Convert ( t : QuadDobl_Complex_Polynomials.Term )
                   return QuadDobl_Complex_Monomials.Monomial;

  -- DESCRIPTION :
  --   Returns the monomial corresponding to the term t,
  --   in double, double double, and quad double precision.

  function Is_Zero ( v : Standard_Natural_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all entries in v are zero.

  function Count_Constant
             ( p : Standard_Complex_Polynomials.Poly ) return integer32;
  function Count_Constant
             ( p : DoblDobl_Complex_Polynomials.Poly ) return integer32;
  function Count_Constant
             ( p : QuadDobl_Complex_Polynomials.Poly ) return integer32;

  -- DESCRIPITON :
  --   Returns 1 if p has a nonzero constant term, returns 0 otherwise.

  function Convert ( p : Standard_Complex_Polynomials.Poly )
                   return Standard_Monomial_Vectors.Polynomial;
  function Convert ( p : DoblDobl_Complex_Polynomials.Poly )
                   return DoblDobl_Monomial_Vectors.Polynomial;
  function Convert ( p : QuadDobl_Complex_Polynomials.Poly )
                   return QuadDobl_Monomial_Vectors.Polynomial;

  -- DESCRIPTION :
  --   Returns the polynomial corresponding to the polynomial p,
  --   in double, double double, and quad double precision.

  function Convert ( p : Standard_Complex_Polynomials.Poly )
                   return Standard_Monomial_Vectors.Link_to_Polynomial;
  function Convert ( p : DoblDobl_Complex_Polynomials.Poly )
                   return DoblDobl_Monomial_Vectors.Link_to_Polynomial;
  function Convert ( p : QuadDobl_Complex_Polynomials.Poly )
                   return QuadDobl_Monomial_Vectors.Link_to_Polynomial;

  -- DESCRIPTION :
  --   Returns the polynomial corresponding to the polynomial p,
  --   in double, double double, and quad double precision.

  function Convert ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                   return Standard_Polynomial_Vectors.System;
  function Convert ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                   return DoblDobl_Polynomial_Vectors.System;
  function Convert ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                   return QuadDobl_Polynomial_Vectors.System;

  -- DESCRIPTION :
  --   Returns the polynomial vector corresponding to the system p,
  --   in double, double double, and quad double precision.

  function Convert ( p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys )
                   return Standard_Polynomial_Vectors.Link_to_System;
  function Convert ( p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys )
                   return DoblDobl_Polynomial_Vectors.Link_to_System;
  function Convert ( p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys )
                   return QuadDobl_Polynomial_Vectors.Link_to_System;

  -- DESCRIPTION :
  --   Returns the polynomial vector corresponding to the system p,
  --   in double, double double, and quad double precision.

end System_Vector_Convertors;
