with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Complex_Vector_Series;

package TripDobl_Random_Series_Vectors is

-- DESCRIPTION :
--   Exports functions that return vectors of random power series,
--   truncated to the given degree,
--   with coefficients in triple double precision.

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return TripDobl_Complex_Series_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random series
  --   of the given degree.

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return TripDobl_Complex_Vector_Series.Vector;

  -- DESCRIPTION :
  --   Returns a series of degree d with have coefficient vectors
  --   of range first..last, with randomly generated complex numbers.

end TripDobl_Random_Series_Vectors;
