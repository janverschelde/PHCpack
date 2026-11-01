with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_rpSeries_Monomials;         use Double_rpSeries_Monomials;

package Test_Double_rpSeries_Monomials is

-- DESCRIPTION :
--   Tests operations on monomials with real power series as coefficients.

  function Random_Monomial
             ( dim,tdx,low,upp : integer32 ) return Monomial;
  function Random_Monomial
             ( dim,tdx,low,upp : integer32 ) return Link_to_Monomial;

  -- DESCRIPTION :
  --   Returns a random monomial in dim variables,
  --   with exponents ranging between low and upp,
  --   and a coefficient series truncated at index tdx.

  function Random_Monomials
             ( nbr,dim,tdx,low,upp : integer32 ) return Monomial_Array;

  -- DESCRIPTION :
  --   Returns an array of nbr monomials in dim variables,
  --   with coefficient series truncated at index tdx,
  --   and exponents ranging between low and upp.

  procedure Write ( mon : in Monomial );
  procedure Write ( mon : in Link_to_Monomial );

  -- DESCRIPTION :
  --   Writes the coefficient series and the exponents of the monomial.

  procedure Write ( mons : in Monomial_Array );

  -- DESCRIPTION :
  --   Writes the monomials in the array mons.

  procedure Test_Gradient ( mon : in Link_to_Monomial );

  -- DESCRIPTION :
  --   Shows the components of the gradient of the monomial mon.

  procedure Test_Random_Make ( dim,tdx : in integer32 );

  -- DESCRIPTION :
  --   Tests the making of a random monomial of dimension dim
  --   and coefficient series truncated at index tdx.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the two parameters of a monomial
  --   and then runs the tests.

end Test_Double_rpSeries_Monomials;
