with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with TripDobl_Complex_Series;            use TripDobl_Complex_Series;

package TripDobl_Complex_Random_Series is

-- DESCRIPTION :
--   Exports functions that return random power series,
--   truncated to the given degree,
--   with coefficients in triple double precision.

  function Random_Series ( degree : integer32 ) return Series;
  function Random_Series ( degree : integer32 ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a series of the given degree,
  --   with randomly generated complex coefficients.

end TripDobl_Complex_Random_Series;
