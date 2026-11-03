with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_rpSeries_PolySystems;       use Double_rpSeries_PolySystems;

package Test_Double_rpSeries_PolySystems is

-- DESCRIPTION :
--   Tests the operations on Laurent polynomial systems with coefficients
--   that are series with real exponents.

  procedure Write ( p : in Laur_Sys );

  -- DESCRIPTION :
  --   Writes the polynomial system p.

  procedure Test_Random_Make ( dim,tdx,nbr : in integer32 );

  -- DESCRIPTION :
  --   Tests the making of a nbr random monomials of dimension dim
  --   and coefficient series truncated at index tdx,
  --   for a system of dim equations.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the parameters and then runs the tests.

end Test_Double_rpSeries_PolySystems;
