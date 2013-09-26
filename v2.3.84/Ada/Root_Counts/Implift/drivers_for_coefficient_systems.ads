with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Drivers_for_Coefficient_Systems is

-- DESCRIPTION :
--   A coefficient system is a system with the same supports as the target
--   system, but with coefficients that may be changed into random ones.
--   This package provides drivers to set up the polyhedral continuation.

  procedure Driver_for_Coefficient_System
               ( file : in file_type; p : in Poly_Sys; k : in natural32;
                 byebye : in boolean;
                 q : out Poly_Sys; qfile,solsfile : in out file_type;
                 tosolve,ranstart,contrep : out boolean );
  procedure Driver_for_Coefficient_System
               ( file : in file_type; p : in Laur_Sys; k : in natural32;
                 byebye : in boolean;
                 q : out Laur_Sys; qfile,solsfile : in out file_type;
                 tosolve,ranstart,contrep : out boolean );

  -- DESCRIPTION :
  --   Prepares the settings for polyhedral continuation.

  -- ON ENTRY :
  --   file      output file, must be opened for output;
  --   p         a polynomial system;
  --   k         if k > 0, then q will not be written on file;
  --   byebye    if true, then a bye-bye message will appear on screen,
  --             if false, then no bye-bye.

  -- ON RETURN :
  --   q         system to solve, if tosolve is true;
  --   qfile     output file for q, is opened for output;
  --   solsfile  output file for solutions of q, is opened for output;
  --   tosolve   true if polyhedral will be applied, false otherwise;
  --   ranstart  true if q is random coefficient start system,
  --             false otherwise;
  --   contrep   true if output information is wanted during continuation.

end Drivers_for_Coefficient_Systems;
