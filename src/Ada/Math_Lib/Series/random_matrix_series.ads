with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Dense_Matrix_Series;
with DoblDobl_Dense_Matrix_Series;
with QuadDobl_Dense_Matrix_Series;

package Random_Matrix_Series is

-- DESCRIPTION :
--   This packages defines random matrix series for testing purposes.

  function Standard_Random_Matrix_Series
             ( deg,dim,lower,upper : integer32 )
             return Standard_Dense_Matrix_Series.Matrix;

  -- DESCRIPTION :
  --   Returns a random matrix series, with integer coefficients 
  --   in [lower, upper], of the given degree deg and dimension dim.
  --   The coefficients are stored in standard double precision.

  function Standard_Random_Matrix_Series
             ( deg,dim : integer32 )
             return Standard_Dense_Matrix_Series.Matrix;

  -- DESCRIPTION :
  --   Returns a random matrix series, with random complex coefficients
  --   of modulus one, of the given degree deg and dimension dim.
  --   The coefficients are stored in standard double precision.

  function DoblDobl_Random_Matrix_Series
             ( deg,dim,lower,upper : integer32 )
             return DoblDobl_Dense_Matrix_Series.Matrix;

  -- DESCRIPTION :
  --   Returns a random matrix series, with integer coefficients
  --   in [lower, upper], of the given degree deg and dimension dim.
  --   The coefficients are stored in double double precision.

  function DoblDobl_Random_Matrix_Series
             ( deg,dim : integer32 )
             return DoblDobl_Dense_Matrix_Series.Matrix;

  -- DESCRIPTION :
  --   Returns a random matrix series, with random complex coefficients
  --   of modulus one, of the given degree deg and dimension dim.
  --   The coefficients are stored in double double precision.

  function QuadDobl_Random_Matrix_Series
             ( deg,dim,lower,upper : integer32 )
             return QuadDobl_Dense_Matrix_Series.Matrix;

  -- DESCRIPTION :
  --   Returns a random matrix series, with integer coefficients
  --   in [lower, upper], of the given degree deg and dimension dim.
  --   The coefficients are stored in quad double precision.

  function QuadDobl_Random_Matrix_Series
             ( deg,dim : integer32 )
             return QuadDobl_Dense_Matrix_Series.Matrix;

  -- DESCRIPTION :
  --   Returns a random matrix series, with random complex coefficients
  --   of modulus one, of the given degree deg and dimension dim.
  --   The coefficients are stored in quad double precision.

end Random_Matrix_Series;
