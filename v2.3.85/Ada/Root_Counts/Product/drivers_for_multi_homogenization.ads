with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Drivers_for_Multi_Homogenization is

  procedure Multi_Homogenization_Info;

  -- DESCRIPTION :
  --   Displays information about multi-homogenization on screen.

  procedure Driver_for_Multi_Homogenization
               ( file : in file_type; p : in Poly_Sys; b : in out natural32;
                 q : out Poly_Sys; rq : out Prod_Sys;
                 qsols : out Solution_List );

  -- DESCRIPTION :
  --   This is an interactive driver for multi-homogenization.

  -- ON ENTRY :
  --   file      to write diagnostics on;
  --   p         a polynomial system.

  -- ON RETURN :
  --   b         a bound based on the degree structure;
  --   q         a random product start system;
  --   rq        product format of the start system q;
  --   qsols     the solutions of q.

end Drivers_for_Multi_Homogenization;
