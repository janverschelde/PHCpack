with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;          use Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Monomial_Maps;            use Standard_Monomial_Maps;

package Standard_Monomial_Map_Circuits is

-- DESCRIPTION :
--   The circuits on the tropisms of monomial maps give 
--   the defining equations for the maps.

  function Equation ( map : Monomial_Map; c : Vector ) return Poly;

  -- DESCRIPTION :
  --   Returns the equation defined by the kernel vector of the circuit c
  --   using the coefficients of the map.  The equation is of the form
  --   x^c = e, where the powers of x are the coefficients of c.

  function Equations ( map : Monomial_Map; c : List ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns all the equations defined by the circuits in c.

  function Circuits ( map : Monomial_Map; d : integer32 ) return List;

  -- DESCRIPTION :
  --   Returns the list of circuits defined by the monomials in the map.
  --   Every circuit in the list is spanned by d+1 points.

end Standard_Monomial_Map_Circuits;
