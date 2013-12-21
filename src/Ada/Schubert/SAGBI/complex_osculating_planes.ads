with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;

package Complex_Osculating_Planes is

-- DESCRIPTION :
--   This package provides routines for generating osculating planes,
--   to test a variant of the Shapiro-Shapiro conjectures.

  function Standard_Basis
             ( n,d : natural32; s : Complex_Number ) return Matrix;

  -- DESCRIPTION :
  --   Returns a d-plane in n-space using the standard monomial basis.
  --   The first column contains the first n standard monomials.
  --   The kth column contains the kth derivative, scaled such that the
  --   diagonal elements are all one.  All polynomials are evaluated at s.

--  function Chebychev_Basis
--             ( n,d : natural32; s : Complex_Number ) return Matrix;

  -- DESCRIPTION :
  --   Returns a d-plane in n-space using the Chebychev polynomials.
  --   The first column contains the first n Chebychev polynomials.
  --   The kth column contains the kth derivative, scaled such that the
  --   diagonal elements are all one.  All polynomials are evaluated at s.

 -- function Orthogonal_Basis
 --            ( n,d : natural32 ; s : Complex_Number ) return Matrix;

  -- DESCRIPTION :
  --   Returns a d-plane in n-space as an orthognal matrix.

end Complex_Osculating_Planes;
