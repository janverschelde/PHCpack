with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Double_Real_Power_Series;
with Double_rpSeries_Vectors;

package Double_rpSeries_Monomials is

-- DESCRIPTION :
--   Applies the reverse mode of algorithmic differentiation
--   to monomials with real power series in double precision.

-- A monomial has two parameters:
-- (1) the number of variables, called the dimension dim; and
-- (2) the index of truncation of the coefficient series tdx.

  type Monomial ( dim,tdx : integer32 ) is record
     cff : Double_Real_Power_Series.Series(tdx);
     deg : Standard_Integer_Vectors.Vector(1..dim);
  end record;
  type Link_to_Monomial is access Monomial;

-- CONSTRUCTORS :

  function Make ( cff : Double_Real_Power_Series.Series;
                  deg : Standard_Integer_Vectors.Vector )
                return Monomial;
  function Make ( cff : Double_Real_Power_Series.Series;
                  deg : Standard_Integer_Vectors.Vector )
                return Link_to_Monomial;

  -- DESCRIPTION :
  --   Given a coefficient series and an exponent vector,
  --   returns the monomial with this coefficient and exponents.

  function Make ( x : complex_number; dim,tdx : integer32 )
                return Monomial;
  function Make ( x : complex_number; dim,tdx : integer32 )
                return Link_to_Monomial;

  -- DESCRIPTION :
  --   Returns the representation of the number x as the
  --   leading coefficient of a series truncated at index tdx,
  --   and of a monomial of dim variables, all with exponent zero.

  function Diff ( mon : Monomial; idx : integer32 ) return Monomial;
  function Diff ( mon : Link_to_Monomial;
                  idx : integer32 ) return Link_to_Monomial;

  -- DESCRIPTION :
  --   Returns the derivative of the monomial with respect to the index idx.

  -- REQUIRED : idx is in the range of mon.deg.

-- SELECTORS :

  function Eval ( mon : Monomial; x : Double_rpSeries_Vectors.Vector )
                return Double_Real_Power_Series.Link_to_Series;
  function Eval ( mon : Link_to_Monomial; x : Double_rpSeries_Vectors.Vector )
                return Double_Real_Power_Series.Link_to_Series;

  -- DESCRIPTION :
  --   Returns the value of the monomial evaluated at x.

-- DESTRUCTOR :

  procedure Clear ( mon : in out Link_to_Monomial );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the monomial mon.

end Double_rpSeries_Monomials;
