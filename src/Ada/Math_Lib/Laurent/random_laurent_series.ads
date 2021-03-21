with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecVecs;

package Random_Laurent_Series is

-- DESCRIPTION :
--   A vector of Laurent series is represented by
--   (1) a vector of leading exponents; and
--   (2) a vector of vectors of complex number with the coefficients
--   of the Laurent series.
--   A matrix of Laurent series is represented by 
--   (1) a matrix of leading exponents; and
--   (2) a vector of vector of complex numbers with the coefficients
--   of the Laurent series.
--   The procedures below generate coefficient vector for general,
--   lower triangular, and upper triangular matrices of Laurent series.

  procedure Random_Series_Coefficients
              ( dim,deg : in integer32;
                cff : out Standard_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Returns in cff the coefficients of dim series of degree deg,
  --   randomly generated on the complex unit circle.
  --   Does all allocations for cff.

  procedure Random_Vector
              ( dim,deg,low,upp : in integer32;
                e : out Standard_Integer_Vectors.Vector;
                c : out Standard_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Returns in (e, c) a random vector of Laurent series.

  -- ON ENTRY :
  --   dim      the dimension of the vector;
  --   deg      the degree of the series;
  --   low      the lower bound on the leading exponents;
  --   upp      the upper bound on the leading exponents;
  --   e        a vector of range 1..dim;
  --   c        a vector of range 1..dim.

  -- ON RETURN :
  --   e        leading exponents of each series in the vector;
  --   c        coefficients of the Laurent series.

  procedure Random_VecVecVec
              ( v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Given a fully allocated 3-dimensional v,
  --   fills it up with random complex numbers on the unit circle.

  procedure Random_Lower_VecVecVec
              ( v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Given a fully allocated 3-dimensional v,
  --   fills it up with random complex numbers on the unit circle,
  --   but only on the lower triangular part, below the diagonal.
  --   The diagonal elements are set to one.

  procedure Random_Upper_VecVecVec
              ( v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Given a fully allocated 3-dimensional v,
  --   fills it up with random complex numbers on the unit circle,
  --   but only on the diagonal and the upper triangular part.

  procedure Random_Matrix
              ( nbrows,nbcols : in natural32; low,upp : in integer32;
                e : out Standard_Integer_Matrices.Matrix;
                v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Returns in (e, v) a random matrix of Laurent series.

  -- ON ENTRY :
  --   nbrows   the number of rows of the exponent matrix e;
  --   nbcols   the number of columns of the exponent matrix e;
  --   low      the lower bound on the leading exponents;
  --   upp      the upper bound on the leading exponents.
  --   e        e'range(1) = 1..rows and e'range(2) = 1..cols;
  --   v        a fully allocated 3-dimensional array,
  --            compatible with e.

  -- ON RETURN :
  --   e        leading exponents in the range low..upp;
  --   v        coefficients of the Laurent series.

  procedure Random_Lower_Matrix
              ( nbrows,nbcols : in natural32; low,upp : in integer32;
                e : out Standard_Integer_Matrices.Matrix;
                v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Returns in (e, v) a random lower triangular matrix of Laurent series,
  --   with ones on the diagonal.

  -- ON ENTRY :
  --   nbrows   the number of rows of the exponent matrix e;
  --   nbcols   the number of columns of the exponent matrix e;
  --   low      the lower bound on the leading exponents;
  --   upp      the upper bound on the leading exponents.
  --   e        e'range(1) = 1..rows and e'range(2) = 1..cols;
  --   v        a fully allocated 3-dimensional array,
  --            compatible with e.

  -- ON RETURN :
  --   e        leading exponents in the range low..upp;
  --   v        coefficients of the Laurent series.

  procedure Random_Upper_Matrix
              ( nbrows,nbcols : in natural32; low,upp : in integer32;
                e : out Standard_Integer_Matrices.Matrix;
                v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Returns in (e, v) a random upper triangular matrix of Laurent series,
  --   with general, nonzero elements on the diagonal.

  -- ON ENTRY :
  --   nbrows   the number of rows of the exponent matrix e;
  --   nbcols   the number of columns of the exponent matrix e;
  --   low      the lower bound on the leading exponents;
  --   upp      the upper bound on the leading exponents.
  --   e        e'range(1) = 1..rows and e'range(2) = 1..cols;
  --   v        a fully allocated 3-dimensional array,
  --            compatible with e.

  -- ON RETURN :
  --   e        leading exponents in the range low..upp;
  --   v        coefficients of the Laurent series.

end Random_Laurent_Series;
