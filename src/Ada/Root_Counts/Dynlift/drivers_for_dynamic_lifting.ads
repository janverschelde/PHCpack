with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Drivers_for_Dynamic_Lifting is

  procedure Dynamic_Lifting_Info;

  -- DESCRIPTION :
  --   Displays information on dynamic lifting on screen.

  procedure Driver_for_Dynamic_Mixed_Volume_Computation 
                ( file : in file_type; p : in Poly_Sys; byebye : in boolean;
                  q : out Poly_Sys; qsols : out Solution_List;
                  mv : out natural32 );
  procedure Driver_for_Dynamic_Mixed_Volume_Computation 
                ( file : in file_type; p : in Laur_Sys; byebye : in boolean;
                  q : out Laur_Sys; qsols : out Solution_List;
                  mv : out natural32 );
  
  -- DESCRIPTION :
  --   This procedure presents an interactive driver for the computation
  --   of the mixed volume.

  -- ON ENTRY :
  --   file       output file, must be opened for output;
  --   p          a polynomial system.

  -- ON RETURN :
  --   q          a start system with randomly choosen coefficients,
  --              which can be used in a coefficient homotopy;
  --   qsols      the solutions of q;
  --   mv         mixed volume of p and the number of solutions of q.

end Drivers_for_Dynamic_Lifting;
