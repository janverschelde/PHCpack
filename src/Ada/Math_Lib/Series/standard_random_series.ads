with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Dense_Series;              use Standard_Dense_Series;

package Standard_Random_Series is

-- DESCRIPTION :
--   Exports functions that return random power series,
--   truncated to the given order.

  function Random_Series ( order : integer32 ) return Series;

  -- DESCRIPTION :
  --   Returns a series of the given order, with random coefficient,
  --   on the unit circle on the complex plane.

end Standard_Random_Series;
