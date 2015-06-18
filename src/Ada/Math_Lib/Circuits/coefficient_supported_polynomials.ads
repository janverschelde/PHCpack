with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Multprec_Complex_Polynomials;

package Coefficient_Supported_Polynomials is

-- DESCRIPTION :
--   This package allows to define polynomials by vectors of supports
--   and corresponding coefficients.

  procedure Split_Common_Factor
              ( e : in Standard_Natural_Vectors.Vector;
                f,b : out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Splits the exponent vector e into a common factor f
  --   and a bit vector b of zeroes and ones.
  --   For all i in e'range: e(i) = f(i) + b(i), b(i) is 0 or 1.

  -- REQUIRED : f'range = b'range = e'range.

  procedure Split_Common_Factors
              ( e : in Standard_Natural_VecVecs.VecVec;
                f,b : out Standard_Natural_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Splits the exponents vectors in e into common factors f
  --   and bit vectors b.

  -- REQUIRED : f'range = b'range = e'range.

  procedure Split_Common_Factors
              ( e : in Standard_Natural_VecVecs.VecVec;
                f,b : out Standard_Natural_VecVecs.VecVec;
                nof : out boolean );

  -- DESCRIPTION :
  --   Extended version of the split that checks whether there is
  --   a common factor different from zero.
  --   If on return nof is true, then all vectors in f are zero vectors,
  --   and thus: if nof, then f can be ignored.
  --   Otherwise, if nof is false, then not all vectors in f are zero
  --   and there are nontrivial common factors.

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
  procedure Coefficients_and_Supports
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                c : out DoblDobl_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec );
  procedure Coefficients_and_Supports
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                c : out QuadDobl_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec );
  procedure Coefficients_and_Supports
              ( p : in Multprec_Complex_Polynomials.Poly;
                c : out Multprec_Complex_Vectors.Vector;
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

  function Create_Multprec_Polynomial
              ( e : Standard_Natural_VecVecs.VecVec )
              return Multprec_Complex_Polynomials.Poly;
  function Create_Multprec_Polynomial
              ( c : Multprec_Complex_Vectors.Vector;
                e : Standard_Natural_VecVecs.VecVec )
              return Multprec_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial which is the sum of monomials with
  --   exponents in e and with constant coefficients all equal to one
  --   if c is unspecified, otherwise the polynomial on return is the
  --   sum of the terms c(i)*x^e(i), for i in c'range = e'range.

end Coefficient_Supported_Polynomials;
