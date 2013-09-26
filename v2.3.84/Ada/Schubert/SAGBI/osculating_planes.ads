with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;

package Osculating_Planes is

-- DESCRIPTION :
--   This package provides routines for generating osculating planes.

  function Standard_Basis ( n,d : natural32; s : double_float ) return Matrix;

  -- DESCRIPTION :
  --   Returns a d-plane in n-space using the standard monomial basis.
  --   The first column contains the first n standard monomials.
  --   The kth column contains the kth derivative, scaled such that the
  --   diagonal elements are all one.  All polynomials are evaluated at s.

  function Chebychev_Basis ( n,d : natural32; s : double_float ) return Matrix;

  -- DESCRIPTION :
  --   Returns a d-plane in n-space using the Chebychev polynomials.
  --   The first column contains the first n Chebychev polynomials.
  --   The kth column contains the kth derivative, scaled such that the
  --   diagonal elements are all one.  All polynomials are evaluated at s.

  function Orthogonal_Basis
             ( n,d : natural32; s : double_float ) return Matrix;

  -- DESCRIPTION :
  --   Returns a d-plane in n-space as an orthognal matrix.

  procedure Sampled_Chebychev_Basis
                ( n,d,m : in natural32; mat : out Matrix;
                  s,ratio : out double_float );

  -- DESCRIPTION :
  --   Generates m values for s and keeps the one with lowest ratio
  --   max/min for all d-minors in the matrix.

end Osculating_Planes;
