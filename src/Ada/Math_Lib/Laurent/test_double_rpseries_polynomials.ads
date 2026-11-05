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

  function Random_Polynomial ( dim,tdx,nbr,low,upp : integer32 ) return Poly;

  -- DESCRIPTION :
  --   Returns a random polynomial with dim variables,
  --   coefficient series truncated at index tdx,
  --   with nbr monomials and exponents ranging from low to upp.

  procedure Test_Random_Make ( dim,tdx,nbr : in integer32 );

  -- DESCRIPTION :
  --   Tests the making of a nbr random monomials of dimension dim
  --   and coefficient series truncated at index tdx.

  procedure Test_Random_Eval ( dim,tdx,nbr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random polynomial of dim variables,
  --   with nbr monomials and coefficients truncated at index tdx.
  --   Evaluates this polynomial at a random vector of series twice:
  --   once as a monomial array and once as a polynomial.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the parameters of a polynomial
  --   and then runs the tests.

end Test_Double_rpSeries_Polynomials;
