with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

procedure Driver_for_Winding_Numbers
             ( file : in file_type; p : in Poly_Sys;
               sols : in out Solution_List );

-- DESCRIPTION :
--   Interactive driver for the computation of winding numbers.
--   The user will be asked to define a homotopy.

-- ON ENTRY :
--   file       file to write intermediate results on;
--   p          a polynomial system;
--   sols       solution list, with t < target value for continuation.

-- ON RETURN :
--   sols       refined solution list with appropriate winding number.
