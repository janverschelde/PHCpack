with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_VecMats;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with PentDobl_Complex_Matrices;
with PentDobl_Complex_VecMats;
with OctoDobl_Complex_Matrices;
with OctoDobl_Complex_VecMats;
with DecaDobl_Complex_Matrices;
with DecaDobl_Complex_VecMats;

package DecaDobl_Complex_Matrices_cv is

-- DESCRIPTION :
--   Functions to convert matrices in deca double precision
--   to double, double double, triple double, quad double, penta double,
--   and octo double precision.

  function to_octo_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return OctoDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the octo double equivalent to the matrix A.

  function to_penta_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return PentDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the penta double equivalent to the matrix A.

  function to_quad_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the quad double equivalent to the matrix A.

  function to_triple_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return TripDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the triple double equivalent to the matrix A.

  function to_double_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the double double equivalent to the matrix A.

  function to_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the double precision equivalent to the matrix A.

  function to_octo_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return OctoDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Converts every matrix in A to octo double precision.

  function to_penta_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return PentDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Converts every matrix in A to penta double precision.

  function to_quad_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return QuadDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Converts every matrix in A to quad double precision.

  function to_triple_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return TripDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Converts every matrix in A to triple double precision.

  function to_double_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return DoblDobl_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Converts every matrix in A to double double precision.

  function to_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return Standard_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Converts every matrix in A to double precision.

end DecaDobl_Complex_Matrices_cv;
