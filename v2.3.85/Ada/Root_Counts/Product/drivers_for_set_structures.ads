with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Drivers_for_Set_Structures is

  procedure Set_Structure_Info;

  -- DESCRIPTION :
  --   Displays information on set structures on screen.

  procedure Read_Set_Structure ( n : in natural32 );

  -- DESCRIPTION :
  --   Allows the user to give a set structure for a system of n equations.

  procedure Driver_for_Set_Structure
               ( file : in file_type; p : in Poly_Sys;
                 b : in out natural32; lpos : in out List;
                 q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Allows the interactive computation of a generalized Bezout number,
  --   with an optional construction of a start system.

  -- ON ENTRY :
  --   file      output file;
  --   p         a polynomial system.

  -- ON RETURN :
  --   b         a bound based on the set structure;
  --   lpos      a list of positions indicating the acceptable classes;
  --   q         a random product start system;
  --   qsols     the solutions of q.

end Drivers_for_Set_Structures;
