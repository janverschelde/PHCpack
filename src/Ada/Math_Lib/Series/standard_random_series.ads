with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Dense_Series;              use Standard_Dense_Series;
with Standard_Dense_Series_Vectors;      use Standard_Dense_Series_Vectors;
with Standard_Dense_Series_VecVecs;      use Standard_Dense_Series_VecVecs;
with Standard_Dense_Series_Matrices;     use Standard_Dense_Series_Matrices;

package Standard_Random_Series is

-- DESCRIPTION :
--   Exports functions that return random power series,
--   truncated to the given order,
--   with coefficients in standard double precision.

  function Random_Series ( order : integer32 ) return Series;

  -- DESCRIPTION :
  --   Returns a series of the given order, with random coefficients,
  --   on the unit circle on the complex plane.

  function Random_Series_Vector
             ( first,last,order : integer32 ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random series
  --   of the given order.

  function Random_Series_VecVec
             ( vvfirst,vvlast,first,last,order : integer32 ) return VecVec;

  -- DESCRIPTION :
  --   Returns a vector of vectors with random series of the given order.
  --   The outer range is vvfirst..vvlast and the inner range of
  --   the vectors is first..last.

  function Random_Series_Matrix
             ( rowfirst,rowlast,columnfirst,columnlast,order : integer32 )
             return Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range(1) rowfirst..rowlast and range(2) 
  --   columnfirst..columnlast with random series of the given order.

end Standard_Random_Series;
