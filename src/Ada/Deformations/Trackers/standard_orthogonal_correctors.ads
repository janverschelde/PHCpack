with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Continuation_Parameters;            use Continuation_Parameters;
with Standard_Continuation_Data;         use Standard_Continuation_Data;

package Standard_Orthogonal_Correctors is

-- DESCRIPTION :
--   This package contains implementations for a Gauss-Newton method
--   applied as corrector in an increment-and-fix continuation.
--   The corrector procedures comes in four flavors, depending on
--   Silent/Reporting: without or with intermediate output; and
--   QRLS/SVD: least squares after QR decomposition (QRLS), or
--   singular value decomposition (SVD) on the evaluated Jacobian matrix.
--   The SVD version computes the condition number.

  generic
    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Silent_QRLS_Corrector
              ( n : in integer32; s : in out Solu_Info; c : in Corr_Pars );

  -- DESCRIPTION :
  --   The predicted solution of the system H(x,t)=0 are corrected,
  --   solving the linear system for the update in the least squares
  --   sense with a QR decomposition on the Jacobian matrix.
  --   This version does not write any intermediate output.

  -- ON ENTRY :
  --   n        the number of equations in the homotopy;
  --   s        the predicted values for the solutions;
  --   c        the corrector parameters.
 
  -- ON RETURN :
  --   s        the computed solution.

  generic
    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Reporting_QRLS_Corrector
              ( file : in file_type;
                n : in integer32; s : in out Solu_Info; c : in Corr_Pars );

  -- DESCRIPTION :
  --   The predicted solution of the system H(x,t)=0 are corrected,
  --   solving the linear system for the update in the least squares
  --   sense with a QR decomposition on the Jacobian matrix.
  --   This version does write intermediate output to file.

  -- ON ENTRY :
  --   file     must be opened for output of intermediate results;
  --   n        the number of equations in the homotopy;
  --   s        the predicted values for the solutions;
  --   c        the corrector parameters.
 
  -- ON RETURN :
  --   s        the computed solution.

  generic
    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Silent_SVD_Corrector
              ( n : in integer32; s : in out Solu_Info; c : in Corr_Pars );

  -- DESCRIPTION :
  --   The predicted solution of the system H(x,t)=0 are corrected,
  --   solving the linear system for the update in the least squares
  --   sense with a singular value decomposition on the Jacobian matrix.
  --   This version does not write any intermediate output.

  -- ON ENTRY :
  --   n        the number of equations in the homotopy;
  --   s        the predicted values for the solutions;
  --   c        the corrector parameters.
 
  -- ON RETURN :
  --   s        the computed solution.

  generic
    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : Complex_Number ) return Vector;
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
  procedure Reporting_SVD_Corrector
              ( file : in file_type;
                n : in integer32; s : in out Solu_Info; c : in Corr_Pars );

  -- DESCRIPTION :
  --   The predicted solution of the system H(x,t)=0 are corrected,
  --   solving the linear system for the update in the least squares
  --   sense with a singular value decomposition on the Jacobian matrix.
  --   This version does write intermediate output to file.

  -- ON ENTRY :
  --   file     must be opened for output of intermediate results;
  --   n        the number of equations in the homotopy;
  --   s        the predicted values for the solutions;
  --   c        the corrector parameters.
 
  -- ON RETURN :
  --   s        the computed solution.

end Standard_Orthogonal_Correctors;
