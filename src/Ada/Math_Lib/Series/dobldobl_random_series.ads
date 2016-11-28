with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Dense_Series;              use DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Dense_Vector_Series;
with DoblDobl_Dense_Series_VecVecs;      use DoblDobl_Dense_Series_VecVecs;
with DoblDobl_Dense_Series_Matrices;     use DoblDobl_Dense_Series_Matrices;

package DoblDobl_Random_Series is

-- DESCRIPTION :
--   Exports functions that return random power series,
--   truncated to the given degree, with complex coefficients
--   in double double precision.

  function Random_Series ( degree : integer32 ) return Series;

  -- DESCRIPTION :
  --   Returns a series of the given degree, with random coefficients,
  --   on the unit circle on the complex plane.

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return DoblDobl_Dense_Series_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random series
  --   of the given degree.

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return DoblDobl_Dense_Vector_Series.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random series
  --   of the given degree.

  function Random_Series_VecVec
             ( vvfirst,vvlast,first,last,degree : integer32 ) return VecVec;

  -- DESCRIPTION :
  --   Returns a vector of vectors with random series of the given degree.
  --   The outer range is vvfirst..vvlast and the inner range of
  --   the vectors is first..last.

  function Random_Series_Matrix
             ( rowfirst,rowlast,columnfirst,columnlast,degree : integer32 )
             return Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range(1) rowfirst..rowlast and range(2) 
  --   columnfirst..columnlast with random series of the given degree.

end DoblDobl_Random_Series;
