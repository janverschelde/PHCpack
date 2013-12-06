with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Drivers_for_m_Homogenization is

  procedure m_Homogenization_Info;

  -- DESCRIPTION :
  --   Displays information about m-homogenization on screen.

  procedure Driver_for_m_Homogenization
                ( file : in file_type; p : in Poly_Sys; b : in out natural64;
                  q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Computation of an m-homogeneous Bezout number with the option
  --   of constructing an m-homogeneous start system.

  -- ON ENTRY :
  --   file       output file;
  --   p          a polynomial system.

  -- ON RETURN :
  --   b          m-homogeneous Bezout number;
  --   q          m-homogeneous start system;
  --   qsols      solutions of q.

end Drivers_for_m_Homogenization;
