with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;

package Scaling_Methods is

-- DESCRIPTION :
--   Equation and variable scaling is defined for double, double double,
--   quad double precision, and arbitrary multiprecision.

  procedure Display_Info;

  -- DESCRIPTION :
  --   Display information about the scaling procedures on screen.

  procedure Equation_Scaling
              ( file : in file_type;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Equation_Scaling
              ( file : in file_type;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Equation_Scaling
              ( file : in file_type;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Equation_Scaling
              ( file : in file_type;
                p : in out Multprec_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Performs equation scaling on the the system p.
  --   Writes timing information on file.

  procedure Variable_Scaling
              ( file : in file_type;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                basis : out natural32;
                scvc : out Standard_Complex_Vectors.Link_to_Vector );
  procedure Variable_Scaling
              ( file : in file_type;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                basis : out natural32;
                scvc : out DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure Variable_Scaling
              ( file : in file_type;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                basis : out natural32;
                scvc : out QuadDobl_Complex_Vectors.Link_to_Vector );
  procedure Variable_Scaling
              ( file : in file_type;
                p : in out Multprec_Complex_Poly_Systems.Poly_Sys;
                basis : out natural32;
                scvc : out Multprec_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Performs variable scaling on the system p.
  --   Writes timing information on file.

  procedure Write_Results
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scvc : in Standard_Complex_Vectors.Link_to_Vector );
  procedure Write_Results
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scvc : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure Write_Results
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scvc : in QuadDobl_Complex_Vectors.Link_to_Vector );
  procedure Write_Results
              ( file : in file_type;
                p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scvc : in Multprec_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Writes the results of the scaling procedure on file.
  --   These results are the scaled system p, and in case basis /= 0,
  --   the scaling coefficients in the vectors scvc.

  procedure Main ( file : in file_type;
                   p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                   basis : out natural32;
                   scvc : out Standard_Complex_Vectors.Link_to_Vector;
                   verbose : in integer32 := 0 );
  procedure Main ( file : in file_type;
                   p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                   basis : out natural32;
                   scvc : out DoblDobl_Complex_Vectors.Link_to_Vector;
                   verbose : in integer32 := 0 );
  procedure Main ( file : in file_type;
                   p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                   basis : out natural32;
                   scvc : out QuadDobl_Complex_Vectors.Link_to_Vector;
                   verbose : in integer32 := 0 );
  procedure Main ( file : in file_type;
                   p : in out Multprec_Complex_Poly_Systems.Poly_Sys;
                   basis : out natural32;
                   scvc : out Multprec_Complex_Vectors.Link_to_Vector;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   This is an interactive driver for phc running in full mode.

  -- ON ENTRY :
  --   file     file to write intermediate results and diagnostics on;
  --   p        a polynomial system;
  --   verbose  the verbose level.

  -- ON RETURN :
  --   p        the scaled polynomial system;
  --   basis    number basis used for scaling, used as flag:
  --            if basis /= 0, then variable scaling has been applied;
  --   scvc     scaling coefficients, only /= null when basis /= 0.

end Scaling_Methods;
