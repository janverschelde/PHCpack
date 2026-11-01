with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_rpSeries_Monomials;         use Double_rpSeries_Monomials;
with Double_rpSeries_Polynomials;       use Double_rpSeries_Polynomials;

package Test_Double_rpSeries_Polynomials is

-- DESCRIPTION :
--   Tests the operations on Laurent polynomials with coefficients
--   that are series with real exponents.

  procedure Add ( p : in out Poly; mon : in Monomial );
  procedure Add ( p : in out Poly; mon : in Link_to_Monomial );

  -- DESCRIPTION :
  --   Adds the monomial in mon to the polynomial p.

  function Make ( mons : Monomial_Array ) return Poly;

  -- DESCRIPTION :
  --   Makes a polynomial by adding the monomials in the array.

  function Make ( p : Poly ) return Monomial_Array;

  -- DESCRIPTION :
  --   Returns the array of the monomials in p.

  procedure Write ( p : in Poly );

  -- DESCRIPTION :
  --   Writes the monomials in p.

  procedure Test_Random_Make ( dim,tdx,nbr : in integer32 );

  -- DESCRIPTION :
  --   Tests the making of a nbr random monomials of dimension dim
  --   and coefficient series truncated at index tdx.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the two parameters of a monomial
  --   and then runs the tests.

end Test_Double_rpSeries_Polynomials;
