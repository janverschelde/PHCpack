with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

procedure Driver_for_Root_Refining
              ( file : in file_type; scalp,p : in Poly_Sys;
                basis : in natural32;
                scalvec : in Link_to_Vector; sols : in out Solution_List);

-- DESCRIPTION :
--   This is the driver for root refining after the continuation.

-- ON ENTRY :
--   file       to write diagnostics on;
--   scalp      the scaled system;
--   p          the original polynomial system;
--   basis      used for scaling;
--   scalvec    vector of coefficients used for scaling;
--   sols       a list of solutions of scalp.

-- ON RETURN :
--   sols       the refined solutions.
