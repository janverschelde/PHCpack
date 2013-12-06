with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package Degrees_in_Sets_of_Unknowns is

-- DESCRIPTION :
--   This procedure provides routines for computing the degree
--   of a given complex polynomial w.r.t. to a given set.

  function Degree ( t : Term; s : Set ) return integer32;
  function Degree ( p : Poly; s : Set ) return integer32;

  -- DESCRIPTION :
  --   Returns the degree in the given set.

  function Degree_Table ( p : Poly_Sys; z : Partition )
                        return Standard_Integer_Matrices.Matrix;
  function Degree_Table ( p : Poly_Sys; z : Partition )
                        return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   The element (i,j) of the returned matrix contains Degree(p(i),z(j)).

end Degrees_in_Sets_of_Unknowns;
