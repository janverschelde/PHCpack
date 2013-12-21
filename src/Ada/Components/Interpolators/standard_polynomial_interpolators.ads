with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;

package Standard_Polynomial_Interpolators is

-- DESCRIPTION :
--   This package provides a way to construct equations for hypersurfaces
--   by interpolation.  All arithmetic is restricted to machine precision.

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

  function Sample ( p : Poly; m : natural32 ) return VecVec;
  function Sample ( p : Poly; m : natural32;
                    scl : double_float ) return VecVec;

  -- DESCRIPTION :
  --   Returns m vectors x that satisfy p(x) = 0.
  --   If scl is parameter, then the samples will have this magnitude.

  procedure Interpolate ( p : in Poly; v : in VecVec;
                          ip : out Poly; rcond : out double_float );

  -- DESCRIPTION :
  --   Uses Gaussian elimination to solve the interpolation conditions.
  --   An estimate for the inverse of the condition number is returned
  --   in the out variable rcond.

  function Interpolate ( p : Poly; v : VecVec ) return Poly;

  -- DESCRIPTION :
  --   Determines the coefficients of p such that p(x) = 0
  --   for all vectors x in v.  The degree of freedom
  --   equals the number of terms in p minus one.  When v
  --   contains more points, they will simply be ignored.
  --   When v contains fewer points, only the highest degree
  --   terms will be specificied, the rest will remain random.

  function Norm ( p : Poly ) return double_float;

  -- DESCRIPTION :
  --   Returns the norm of the coefficient vector of p.

  function Distance ( p,q : Poly ) return double_float;

  -- DESCRIPTION :
  --   Returns the distance between the two polynomials as the
  --   norm of the difference of the coefficient vectors.

end Standard_Polynomial_Interpolators;
