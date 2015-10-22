with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Brackets;                           use Brackets;
with Bracket_Monomials;                  use Bracket_Monomials;
with Standard_Natural_Matrices;          use Standard_Natural_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;

package Bracket_Expansions is

-- DESCRIPTION :
--   This package provides operations to expand bracket polynomials as
--   complex multivariate polynomials in the matrix indeterminates xij.

  function Expand ( n,d : natural32; b : Bracket ) return Poly;
  function Expand ( n,d : natural32; b : Bracket_Monomial ) return Poly;
  function Expand ( n,d : natural32; b : Bracket_Term ) return Poly;
  function Expand ( n,d : natural32; b : Bracket_Polynomial ) return Poly;

  -- DESCRIPTION :
  --   On return is the expanded bracket polynomial in the xij's, 
  --   where i runs over 1..n and j over 1..d.

  function Localized_Expand ( n,d : natural32; b : Bracket ) return Poly;

  -- DESCRIPTION :
  --   For i >= n-d+1, the variable xij is either 1 or 0, depending on
  --   whether i=j+n-d or not, for i in 1..n and j in 1..d.
  --   This is the standard map to localize a d-plane, a better one
  --   is generated below.

  function Localization_Map ( n,d : natural32 ) return Matrix;

  -- DESCRIPTION :
  --   Returns a localization map for a matrix representing a d-plane
  --   in affine n-space.  The elements of the identity matrix are as
  --   usual represented by zeros and ones.  Every row will have at least
  --   one free element whose entry is marked by two.

  -- REQUIRED : n > d+1.

  function Expand ( locmap : Matrix; b : Bracket ) return Poly;

  -- DESCRIPTION :
  --   Expands a d-by-d minor of the matrix, selecting the rows with
  --   entries in b, respecting the localization map in locmap.
  --   The format of locmap must be as the output of Localization_Map.
  --   The polynomial on return has as many variables as the number of
  --   entries in the matrix locmap.
  --   The number of variables can be reduced by the procedure below.

  procedure Reduce_Variables ( locmap : in Matrix; p : in out Poly );

  -- DESCRIPTION :
  --   Reduces the #variables in the polynomial, removing all variables that
  --   correspond to zeros in the localization map.

end Bracket_Expansions;
