with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with OctoDobl_Complex_Series;            use OctoDobl_Complex_Series;

package OctoDobl_Complex_Random_Series is

-- DESCRIPTION :
--   Exports functions that return random power series,
--   truncated to the given degree,
--   with coefficients in octo double precision.

  function Random_Series ( degree : integer32 ) return Series;
  function Random_Series ( degree : integer32 ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a series of the given degree,
  --   with randomly generated complex coefficients.

end OctoDobl_Complex_Random_Series;
