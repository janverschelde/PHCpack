with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

procedure Driver_for_Own_Start_System
             ( file : in file_type; p : in Poly_Sys;
               q : out Poly_Sys; qsols : in out Solution_List );

-- DESCRIPTION :
--   This procedure implements an interactive driver for reading
--   a start system delivered by user.

-- ON ENTRY :
--   file      to write diagnostics on;
--   p         a polynomial system.

-- ON RETURN :
--   q         a start system based on the chosen root count;
--   qsols     the solutions of q.
