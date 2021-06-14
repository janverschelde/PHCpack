with Standard_Integer_Numbers;       use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;

package Multitasking_Matrix_x_Vector is

-- DESCRIPTION :
--   Provides multitasked versions of matrix-vector multiplications
--   for various types of arithmetic.

  function Silent_Multiply
              ( n : integer32;
                A : Standard_Complex_Matrices.Matrix;
                v : Standard_Complex_Vectors.Vector ) 
              return Standard_Complex_Vectors.Vector;
  function Silent_Multiply
              ( n : integer32;
                A : DoblDobl_Complex_Matrices.Matrix;
                v : DoblDobl_Complex_Vectors.Vector ) 
              return DoblDobl_Complex_Vectors.Vector;
  function Silent_Multiply
              ( n : integer32;
                A : QuadDobl_Complex_Matrices.Matrix;
                v : QuadDobl_Complex_Vectors.Vector ) 
              return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns A*v computed with n tasks without any output.

  function Reporting_Multiply
              ( n : integer32;
                A : Standard_Complex_Matrices.Matrix;
                v : Standard_Complex_Vectors.Vector ) 
              return Standard_Complex_Vectors.Vector;
  function Reporting_Multiply
              ( n : integer32;
                A : DoblDobl_Complex_Matrices.Matrix;
                v : DoblDobl_Complex_Vectors.Vector ) 
              return DoblDobl_Complex_Vectors.Vector;
  function Reporting_Multiply
              ( n : integer32;
                A : QuadDobl_Complex_Matrices.Matrix;
                v : QuadDobl_Complex_Vectors.Vector ) 
              return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns A*v computed with n tasks writing diagnostics to
  --   screen to monitor the progress of the computations.

end Multitasking_Matrix_x_Vector;
