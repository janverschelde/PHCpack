with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecMats;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Vector_Series;
with Standard_Complex_Matrix_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Vector_Series;
with DoblDobl_Complex_Matrix_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Vector_Series;
with QuadDobl_Complex_Matrix_Series;

package Series_Coefficient_Vectors is

-- DESCRIPTION :
--   This package offers functions to extract the coefficients of the
--   series in vectors into the more basic vector of vectors data type.

-- A series vector is a vector of power series.

  function Standard_Series_Coefficients
             ( s : Standard_Complex_Series_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec;
  function DoblDobl_Series_Coefficients
             ( s : DoblDobl_Complex_Series_Vectors.Vector )
             return DoblDobl_Complex_VecVecs.VecVec;
  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Complex_Series_Vectors.Vector )
             return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the coefficients of the series in the vector of vectors.
  --   The range of the k-th vector is 0..s(k).deg.

-- A vector series is power series, where the coefficients are vectors.

  function Standard_Series_Coefficients
             ( s : Standard_Complex_Vector_Series.Vector )
             return Standard_Complex_VecVecs.VecVec;
  function DoblDobl_Series_Coefficients
             ( s : DoblDobl_Complex_Vector_Series.Vector )
             return DoblDobl_Complex_VecVecs.VecVec;
  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Complex_Vector_Series.Vector )
             return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns he coefficients of the series in the vector of vectors.
  --   The range of the vector on return is 0..s.deg, and each vector
  --   has the range 1..n, where n is the dimension of the vectors in
  --   the coefficients of s.  This function returns a copy of s.cff.

-- A matrix series is a power series, where the coefficients are matrices.

  function Standard_Series_Coefficients
             ( s : Standard_Complex_Matrix_Series.Matrix )
             return Standard_Complex_VecMats.VecMat;
  function DoblDobl_Series_Coefficients
             ( s : DoblDobl_Complex_Matrix_Series.Matrix )
             return DoblDobl_Complex_VecMats.VecMat;
  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Complex_Matrix_Series.Matrix )
             return QuadDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Returns a vector of matrices, of range 0..s.deg.
  --   The matrices all have the same dimension, as in s.cff.
  --   A copy of s.cff is returned.

end Series_Coefficient_Vectors;
