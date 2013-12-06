with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;

package Coefficient_Supported_Polynomials is

-- DESCRIPTION :
--   This package allows to define polynomials by vectors of supports
--   and corresponding coefficients.

  function Create_Standard_Polynomial
              ( e : Standard_Natural_VecVecs.VecVec )
              return Standard_Complex_Polynomials.Poly;
  function Create_Standard_Polynomial
              ( c : Standard_Complex_Vectors.Vector;
                e : Standard_Natural_VecVecs.VecVec )
              return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial which is the sum of monomials with
  --   exponents in e and with constant coefficients all equal to one
  --   if c is unspecified, otherwise the polynomial on return is the
  --   sum of the terms c(i)*x^e(i), for i in c'range = e'range.

  procedure Coefficients_and_Supports
              ( p : in Standard_Complex_Polynomials.Poly;
                c : out Standard_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Extracts the coefficients and corresponding exponent vectors
  --   of the polynomial p.

  -- REQUIRED : c'range = e'range = 1..Number_of_Terms(p).

  procedure Coefficients_and_Supports
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                c : out Standard_Complex_VecVecs.VecVec;
                e : out Standard_Natural_VecVecs.Array_of_VecVecs );

  -- DESCRIPTION :
  --   Extracts of the polynomial system p the coefficients and
  --   corresponding exponents that define the monomials.

  -- REQUIRED : c'range = e'range = p'range.

  function Create_DoblDobl_Polynomial
              ( e : Standard_Natural_VecVecs.VecVec )
              return DoblDobl_Complex_Polynomials.Poly;
  function Create_DoblDobl_Polynomial
              ( c : DoblDobl_Complex_Vectors.Vector;
                e : Standard_Natural_VecVecs.VecVec )
              return DoblDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial which is the sum of monomials with
  --   exponents in e and with constant coefficients all equal to one
  --   if c is unspecified, otherwise the polynomial on return is the
  --   sum of the terms c(i)*x^e(i), for i in c'range = e'range.

  function Create_QuadDobl_Polynomial
              ( e : Standard_Natural_VecVecs.VecVec )
              return QuadDobl_Complex_Polynomials.Poly;
  function Create_QuadDobl_Polynomial
              ( c : QuadDobl_Complex_Vectors.Vector;
                e : Standard_Natural_VecVecs.VecVec )
              return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial which is the sum of monomials with
  --   exponents in e and with constant coefficients all equal to one
  --   if c is unspecified, otherwise the polynomial on return is the
  --   sum of the terms c(i)*x^e(i), for i in c'range = e'range.

end Coefficient_Supported_Polynomials;
