with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package Interpolating_Homotopies is

-- DESCRIPTION :
--   This package contains facilities for constructing interpolating
--   homotopies, based on a given m-homogeneous structure.
--   The routines are given in the order in which they should be applied.
--   Null polynomials are ignored, making scaled interpolation possible,
--   the scaling equation, used for generating the interpolating vectors,
--   can be added afterwards.  For linear scalers, the last unknown of the
--   scaling equation should be ignored in those monomials that have degree
--   one in that unknown, to avoid singular interpolation matrices.

  function Dense_Representation
              ( p : Poly_Sys; z : partition ) return Poly_Sys;
  function Dense_Representation
              ( p : Poly_Sys; z : partition; d : Matrix ) return Poly_Sys;

  -- DESCRIPTION :
  --   A dense representation of an m-homogeneous structure is returned.
  --   The coefficients of the polynomials in the returned system are all one.

  function Independent_Representation ( p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   An independent representation of a polynomial system is returned.
  --   This means that the initial term of each polynomial does not occur
  --   in every other polynomial.

  function Independent_Roots ( p : Poly_Sys ) return natural32;
  function Independent_Roots ( p : Poly_Sys; i : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of independent roots the system p can have.
  --   When the ith unknown is given as parameter, the monomials that
  --   have degree one in x_i are not counted.

  -- IMPORTANT NOTE : p must be an independent representation of a polynomial
  --                  system, otherwise the result might not be reliable.

  function Interpolate ( p : Poly_Sys; b : natural32; sols : Solution_List )
                       return Poly_Sys;
  function Interpolate ( p : Poly_Sys; i,b : natural32; sols : Solution_List )
                       return Poly_Sys;

  -- DESCRIPTION :
  --   This routine constructs a start system q with the same monomial
  --   structure as the system p.

  -- ON ENTRY :
  --  p         a polynomial system;
  --  i         monomials with degree one in x_i will be ignored;
  --  b         must equal Independent_Roots(p);
  --  sols      interpolation vectors, Length_Of(sols) = b.

  -- ON RETURN :
  --  q         system that has the given list sols as solutions.

end Interpolating_Homotopies;
