with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_VecMats;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecMats;
with PentDobl_Complex_VecVecs;
with PentDobl_Complex_VecMats;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_VecMats;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_VecMats;
with HexaDobl_Complex_VecVecs;
with HexaDobl_Complex_VecMats;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Vector_Series;
with Standard_Complex_Matrix_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Vector_Series;
with DoblDobl_Complex_Matrix_Series;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Complex_Vector_Series;
with TripDobl_Complex_Matrix_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Vector_Series;
with QuadDobl_Complex_Matrix_Series;
with PentDobl_Complex_Series_Vectors;
with PentDobl_Complex_Vector_Series;
with PentDobl_Complex_Matrix_Series;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_Complex_Vector_Series;
with OctoDobl_Complex_Matrix_Series;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_Complex_Vector_Series;
with DecaDobl_Complex_Matrix_Series;
with HexaDobl_Complex_Series_Vectors;
with HexaDobl_Complex_Vector_Series;
with HexaDobl_Complex_Matrix_Series;

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
  function TripDobl_Series_Coefficients
             ( s : TripDobl_Complex_Series_Vectors.Vector )
             return TripDobl_Complex_VecVecs.VecVec;
  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Complex_Series_Vectors.Vector )
             return QuadDobl_Complex_VecVecs.VecVec;
  function PentDobl_Series_Coefficients
             ( s : PentDobl_Complex_Series_Vectors.Vector )
             return PentDobl_Complex_VecVecs.VecVec;
  function OctoDobl_Series_Coefficients
             ( s : OctoDobl_Complex_Series_Vectors.Vector )
             return OctoDobl_Complex_VecVecs.VecVec;
  function DecaDobl_Series_Coefficients
             ( s : DecaDobl_Complex_Series_Vectors.Vector )
             return DecaDobl_Complex_VecVecs.VecVec;
  function HexaDobl_Series_Coefficients
             ( s : HexaDobl_Complex_Series_Vectors.Vector )
             return HexaDobl_Complex_VecVecs.VecVec;

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
  function TripDobl_Series_Coefficients
             ( s : TripDobl_Complex_Vector_Series.Vector )
             return TripDobl_Complex_VecVecs.VecVec;
  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Complex_Vector_Series.Vector )
             return QuadDobl_Complex_VecVecs.VecVec;
  function PentDobl_Series_Coefficients
             ( s : PentDobl_Complex_Vector_Series.Vector )
             return PentDobl_Complex_VecVecs.VecVec;
  function OctoDobl_Series_Coefficients
             ( s : OctoDobl_Complex_Vector_Series.Vector )
             return OctoDobl_Complex_VecVecs.VecVec;
  function DecaDobl_Series_Coefficients
             ( s : DecaDobl_Complex_Vector_Series.Vector )
             return DecaDobl_Complex_VecVecs.VecVec;
  function HexaDobl_Series_Coefficients
             ( s : HexaDobl_Complex_Vector_Series.Vector )
             return HexaDobl_Complex_VecVecs.VecVec;

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
  function TripDobl_Series_Coefficients
             ( s : TripDobl_Complex_Matrix_Series.Matrix )
             return TripDobl_Complex_VecMats.VecMat;
  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Complex_Matrix_Series.Matrix )
             return QuadDobl_Complex_VecMats.VecMat;
  function PentDobl_Series_Coefficients
             ( s : PentDobl_Complex_Matrix_Series.Matrix )
             return PentDobl_Complex_VecMats.VecMat;
  function OctoDobl_Series_Coefficients
             ( s : OctoDobl_Complex_Matrix_Series.Matrix )
             return OctoDobl_Complex_VecMats.VecMat;
  function DecaDobl_Series_Coefficients
             ( s : DecaDobl_Complex_Matrix_Series.Matrix )
             return DecaDobl_Complex_VecMats.VecMat;
  function HexaDobl_Series_Coefficients
             ( s : HexaDobl_Complex_Matrix_Series.Matrix )
             return HexaDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Returns a vector of matrices, of range 0..s.deg.
  --   The matrices all have the same dimension, as in s.cff.
  --   A copy of s.cff is returned.

end Series_Coefficient_Vectors;
