with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers;

package Fabry_on_Homotopy_Helpers is

-- DESCRIPTION :
--   Helper procedures to compute the Newton-Fabry convergence radius
--   for artificial or natural-parameter homotopies.

  procedure Prompt_for_Parameters
              ( maxit : in out integer32; tol : in out double_float;
                verbose : out boolean );

  -- DESCRIPTION :
  --   Interactive setting of parameters.

  -- ON ENTRY :
  --   maxit    initial number of maximum number of iterations;
  --   tol      initial value for the tolerance.

  -- ON RETURN :
  --   maxit    new value for the maximum number of iterations;
  --   tol      new value for the tolerance;
  --   verbose  true if output during the Newton steps is wanted.

  procedure Prompt_and_Write
              ( file : in file_type; nbtasks : in out natural32;
                maxit : in out integer32; tol : in out double_float;
                verbose : out boolean );

  -- DESCRIPTION :
  --   Prompts for parameters and writes their values to file.
  --   Terminates with closing statement, signaling the end of
  --   the interactive input, as this is useful for larger problems.

  -- ON ENTRY :
  --   nbtasks  initial number of tasks,
  --            if zero, then there is a prompt for the number of tasks;
  --   maxit    initial number of maximum number of iterations;
  --   tol      initial value for the tolerance.

  -- ON RETURN :
  --   nbtasks  new value for the number of tasks (if nonzero);
  --   maxit    new value for the maximum number of iterations;
  --   tol      new value for the tolerance;
  --   verbose  true if output during the Newton steps is wanted.

  procedure Write_Report
              ( file : in file_type; rad,err : in double_float;
                zpt : in Standard_Complex_Numbers.Complex_Number;
                fail : in boolean );
  procedure Write_Report
              ( file : in file_type; rad,err : in double_double;
                zpt : in DoblDobl_Complex_Numbers.Complex_Number;
                fail : in boolean );
  procedure Write_Report
              ( file : in file_type; rad,err : in triple_double;
                zpt : in TripDobl_Complex_Numbers.Complex_Number;
                fail : in boolean );
  procedure Write_Report
              ( file : in file_type; rad,err : in quad_double;
                zpt : in QuadDobl_Complex_Numbers.Complex_Number;
                fail : in boolean );
  procedure Write_Report
              ( file : in file_type; rad,err : in penta_double;
                zpt : in PentDobl_Complex_Numbers.Complex_Number;
                fail : in boolean );
  procedure Write_Report
              ( file : in file_type; rad,err : in octo_double;
                zpt : in OctoDobl_Complex_Numbers.Complex_Number;
                fail : in boolean );
  procedure Write_Report
              ( file : in file_type; rad,err : in deca_double;
                zpt : in DecaDobl_Complex_Numbers.Complex_Number;
                fail : in boolean );
  procedure Write_Report
              ( file : in file_type; rad,err : in hexa_double;
                zpt : in HexaDobl_Complex_Numbers.Complex_Number;
                fail : in boolean );

  -- DESCRIPTION :
  --   Writes the results of the computed Newton-Fabry convergence radius.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   rad      convergence radius;
  --   err      estimated error on the radius;
  --   zpt      location of the nearest singularity;
  --   fail     true if failed to reach the tolerance.

end Fabry_on_Homotopy_Helpers;
