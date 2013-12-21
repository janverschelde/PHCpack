with Standard_Natural_Numbers;          use Standard_Natural_Numbers; 
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_VecVecs;          use Multprec_Complex_VecVecs;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;

package Multprec_Polynomial_Interpolators is

-- DESCRIPTION :
--   This package provides a way to construct equations for hypersurfaces
--   by interpolation.  The accuracy is determined by the number of
--   decimal places of the incoming interpolation points.

  function Number_of_Terms ( d,n : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of terms in a dense polynomial of degree d in
  --   n variables, for which all coefficients are assumed to be nonzero.

  function Create ( d,n,cff : natural32 ) return Poly;

  -- DESCRIPTION :
  --   Returns a dense polynomial of degree d in n variables, with
  --   all monomials up to degree d have a nonzero coefficient.
  --   The choice of coefficients depends on the value of cff:
  --     cff = 1 : all coefficients are equal to one;
  --     cff = 2 : coefficients are randomly generated real numbers;
  --     otherwise : random complex numbers of modulus one are chosen.

  function Sample ( p : Poly; m,sz : natural32 ) return VecVec;

  -- DESCRIPTION :
  --   Returns m vectors x that satisfy p(x) = 0.
  --   The vectors contain numbers of the given size sz. 

  procedure Interpolate1 ( p : in Poly; v : in VecVec; ip : out Poly;
                           invcnd : out Floating_Number );
  procedure Interpolate ( p : in Poly; v : in VecVec; ip : out Poly;
                          invcnd : out Floating_Number );

  -- DESCRIPTION :
  --   Determines the coefficients of p such that p(x) = 0 for all points x
  --   in v.  The degree of freedom equals the number of terms in p minus 1.
  --   When v contains more points, they will simply be ignored.
  --   When v contains fewer points, only the highest degree
  --   terms will be specificied, the rest will remain random.
  --   Interpolate1 assumes the constant term to be one, the other
  --   Interpolate takes the sum of the coefficients equal to one.

  -- ON ENTRY :
  --   p         monomial structure of the interpolating polynomial;
  --   v         interpolation points.

  -- ON RETURN :
  --   ip        interpolating polynomial;
  --   invcnd    inverse of the estimated condition number.

  function Interpolate ( p : Poly; v : VecVec ) return Poly;

  -- DESCRIPTION :
  --   Determines the coefficients of p such that p(x) = 0
  --   for all vectors x in v.  The degree of freedom
  --   equals the number of terms in p minus one.  When v
  --   contains more points, they will simply be ignored.
  --   When v contains fewer points, only the highest degree
  --   terms will be specificied, the rest will remain random.

  function Norm ( p : Poly ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the norm of the coefficient vector of p.

  function Distance ( p,q : Poly ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the distance between the two polynomials as the
  --   norm of the difference of the coefficient vectors.

end Multprec_Polynomial_Interpolators;
