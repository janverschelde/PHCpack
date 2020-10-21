with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_VecMats;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;

package QuadDobl_Complex_Matrices_cv is

-- DESCRIPTION :
--   Functions to convert matrices in quad double precision
--   to double, double double, and triple double precision.

  function to_triple_double
             ( A : QuadDobl_Complex_Matrices.Matrix )
             return TripDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the triple double equivalent to the matrix A.

  function to_double_double
             ( A : QuadDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the double double equivalent to the matrix A.

  function to_double
             ( A : QuadDobl_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the double precision equivalent to the matrix A.

  function to_triple_double
             ( A : QuadDobl_Complex_VecMats.VecMat )
             return TripDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Converts every matrix in A to triple double precision.

  function to_double_double
             ( A : QuadDobl_Complex_VecMats.VecMat )
             return DoblDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Converts every matrix in A to double double precision.

  function to_double
             ( A : QuadDobl_Complex_VecMats.VecMat )
             return Standard_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Converts every matrix in A to double precision.

end QuadDobl_Complex_Matrices_cv;
