with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with TripDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with PentDobl_Complex_Polynomials;
with OctoDobl_Complex_Polynomials;
with DecaDobl_Complex_Polynomials;
with HexaDobl_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package Degrees_in_Sets_of_Unknowns is

-- DESCRIPTION :
--   This procedure provides routines for computing the degree
--   of a given complex polynomial w.r.t. to a given set.

  function Degree ( t : Standard_Complex_Polynomials.Term; s : Set )
                  return integer32;
  function Degree ( t : DoblDobl_Complex_Polynomials.Term; s : Set )
                  return integer32;
  function Degree ( t : TripDobl_Complex_Polynomials.Term; s : Set )
                  return integer32;
  function Degree ( t : QuadDobl_Complex_Polynomials.Term; s : Set )
                  return integer32;
  function Degree ( t : PentDobl_Complex_Polynomials.Term; s : Set )
                  return integer32;
  function Degree ( t : OctoDobl_Complex_Polynomials.Term; s : Set )
                  return integer32;
  function Degree ( t : DecaDobl_Complex_Polynomials.Term; s : Set )
                  return integer32;
  function Degree ( t : HexaDobl_Complex_Polynomials.Term; s : Set )
                  return integer32;

  -- DESCRIPTION :
  --   Returns the degree of the term t in the set s.

  function Degree ( p : Standard_Complex_Polynomials.Poly; s : Set )
                  return integer32;
  function Degree ( p : DoblDobl_Complex_Polynomials.Poly; s : Set )
                  return integer32;
  function Degree ( p : TripDobl_Complex_Polynomials.Poly; s : Set )
                  return integer32;
  function Degree ( p : QuadDobl_Complex_Polynomials.Poly; s : Set )
                  return integer32;
  function Degree ( p : PentDobl_Complex_Polynomials.Poly; s : Set )
                  return integer32;
  function Degree ( p : OctoDobl_Complex_Polynomials.Poly; s : Set )
                  return integer32;
  function Degree ( p : DecaDobl_Complex_Polynomials.Poly; s : Set )
                  return integer32;
  function Degree ( p : HexaDobl_Complex_Polynomials.Poly; s : Set )
                  return integer32;

  -- DESCRIPTION :
  --   Returns the degree of the polynomials p in the set s.

  function Degree_Table ( p : Poly_Sys; z : Partition )
                        return Standard_Integer_Matrices.Matrix;
  function Degree_Table ( p : Poly_Sys; z : Partition )
                        return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   The element (i,j) of the returned matrix contains Degree(p(i),z(j)).

end Degrees_in_Sets_of_Unknowns;
