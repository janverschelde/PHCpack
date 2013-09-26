with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Drivers_for_Implicit_Lifting is

  procedure Implicit_Lifting_Info;

  -- DESCRIPTION :
  --   Displays information on the BKK Bound on screen.

  procedure Driver_for_Mixture_Bezout_BKK
                ( file : in file_type; p : in Poly_Sys; byebye : in boolean;
                  q : out Poly_Sys; qsols : out Solution_List;
                  b : out natural32 );
  procedure Driver_for_Mixture_Bezout_BKK
                ( file : in file_type; p : in Laur_Sys; byebye : in boolean;
                  q : out Laur_Sys; qsols : out Solution_List;
                  b : out natural32 );

  -- DESCRIPTION :
  --   This driver allows to compute a mixture between a generalized Bezout
  --   number and the BKK bound.

  -- ON ENTRY :
  --   file       to write diagnostics on;
  --   p          a polynomial system;
  --   byebye     if true, then a bye-bye message will appear on screen,
  --              if false, then no bye-bye.

  -- ON RETURN :
  --   q          a start system;
  --   qsols      the solutions of q;
  --   b          Bezout-BKK bound for the number of finite solutions of p.

end Drivers_for_Implicit_Lifting;
