with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Corrector_Convolutions is

-- DESCRIPTION :
--   The procedures in this package run Newton's method at a point,
--   on the leading coefficients of systems of convolution circuits,
--   to correct a solution after a predictor step.

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true );
  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true );
  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true );

  -- DESCRIPTION :
  --   Does one Newton step with LU factorization,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write extra output to if verbose;
  --   hom      convolution system for a homotopy;
  --   sol      an initial value for a solution at t = 0;
  --   verbose  flag to indicate if vectors need to be written.

  -- ON RETURN :
  --   sol      the updated solution;
  --   dx       the update vector applied to the solution;
  --   ipvt     pivoting information for the LU factorization;
  --   info     zero is all went well, if nonzero,
  --            then the matrix in hom.vm(0) may be singular.

  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; nrmdx : out double_float; 
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                verbose : in boolean := true );
  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; nrmdx : out double_double; 
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                verbose : in boolean := true );
  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; nrmdx : out quad_double; 
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs several Newton steps using LU factorization,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write extra output to if verbose;
  --   hom      convolution system for a homotopy;
  --   sol      an initial value for a solution at t = 0;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on the updatex dx;
  --   verbose  flag to indicate if vectors need to be written.

  -- ON RETURN :
  --   sol      the updated solution;
  --   nrmdx    max norm of the vector dx;
  --   dx       the update vector applied to the solution;
  --   ipvt     pivoting information for the LU factorization;
  --   info     zero is all went well, if nonzero,
  --            then the matrix in hom.vm(0) may be singular;
  --   fail     true if nrmdx > tol after maxit iterations,
  --            false otherwise.

end Corrector_Convolutions;
