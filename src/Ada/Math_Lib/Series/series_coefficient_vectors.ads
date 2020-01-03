with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;

package Series_Coefficient_Vectors is

-- DESCRIPTION :
--   This package offers functions to extract the coefficients of the
--   series in vectors into the more basic vector of vectors data type.

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

end Series_Coefficient_Vectors;
