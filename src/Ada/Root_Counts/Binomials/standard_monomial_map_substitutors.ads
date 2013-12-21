with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Monomial_Maps;             use Standard_Monomial_Maps;

package Standard_Monomial_Map_Substitutors is

-- DESCRIPTION :
--   This package offers functions to substitute a monomial map
--   with standard complex coefficients into Laurent polynomials.
--   Filters are provided to remove terms with small coefficients
--   in case the map makes the polynomial vanish.

  function Filter ( p : Poly; tol : double_float ) return Poly;

  -- DESCRIPTION :
  --   Copies all terms of p, except those with coefficients of
  --   magnitude less than the tolerance tol.

  function Filter ( p : Laur_Sys; tol : double_float ) return Laur_Sys;

  -- DESCRIPTION :
  --   Filters p, removing all terms with coefficients of magnitude
  --   less than the tolerance tol, omitting null polynomials as well.

  function Subs ( t : Term; map : Monomial_Map ) return Term;

  -- DESCRIPTION :
  --   Returns a term in map.d variables, obtained after substitution
  --   of the map into the term t.

  function Subs ( p : Poly; map : Monomial_Map ) return Poly;

  -- DESCRIPTION :
  --   Substitutes the map into the Laurent polynomial p.

end Standard_Monomial_Map_Substitutors;
