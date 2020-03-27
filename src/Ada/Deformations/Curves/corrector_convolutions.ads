with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
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

-- STEP COEFFICIENTS :
--   The correct is applied to check whether a step size is good.
--   The Step_Coefficient functions compute the constant coefficient
--   of the power series in the convolution circuits, given a value
--   for the step size t, which is evaluated at the power series.

  function Step_Coefficient
              ( c : Standard_Complex_Vectors.Vector; t : double_float )
              return Standard_Complex_Numbers.Complex_Number;
  function Step_Coefficient
              ( c : DoblDobl_Complex_Vectors.Vector; t : double_double )
              return DoblDobl_Complex_Numbers.Complex_Number;
  function Step_Coefficient
              ( c : QuadDobl_Complex_Vectors.Vector; t : quad_double )
              return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Given in c is the coefficient vector of a series
  --   and in t is some step size, t > 0.
  --   Returns the value of the series c evaluated at t,
  --   computed in double, double double, or quad double precision.

  procedure Step_Coefficient
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                t : in double_float );
  procedure Step_Coefficient
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                t : in double_double );
  procedure Step_Coefficient
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                t : in quad_double );

  -- DESCRIPTION :
  --   Replaces the leading coefficient of all series in hom by
  --   the series coefficient evaluated at t,
  --   in double, double double, or quad double precision.

  procedure Step_Coefficient
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                t : in double_float );
  procedure Step_Coefficient
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                t : in double_double );
  procedure Step_Coefficient
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                t : in quad_double );

  -- DESCRIPTION :
  --   Replaces the leading coefficient of all series in hom by
  --   the series coefficient evaluated at t,
  --   in double, double double, and quad double precision.

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
