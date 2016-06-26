with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Dense_Series;              use DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_Vectors;      use DoblDobl_Dense_Series_Vectors;

package DoblDobl_Random_Series is

-- DESCRIPTION :
--   Exports functions that return random power series,
--   truncated to the given order, with complex coefficients
--   in double double precision.

  function Random_Series ( order : integer32 ) return Series;

  -- DESCRIPTION :
  --   Returns a series of the given order, with random coefficients,
  --   on the unit circle on the complex plane.

  function Random_Series_Vector
             ( first,last,order : integer32 ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random series
  --   of the given order.

end DoblDobl_Random_Series;
