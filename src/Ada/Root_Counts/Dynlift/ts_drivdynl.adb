with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Drivers_for_Dynamic_Lifting;        use Drivers_for_Dynamic_Lifting;

procedure ts_drivdynl is

-- DESCRIPTION :
--   This procedure calls the driver to dynamic lifting.

  lp,lq : Link_to_Poly_Sys;
  qsols : Solution_List;
  b : natural32;
  file : file_type;

begin
  new_line;
  put_line("Test on driver for dynamic lifting.");
  new_line;
  get(lp);
  new_line;
  put_line("Reading the name of the output file.");
  Read_Name_and_Create_File(file);
  lq := new Poly_Sys(lp'range);
  Driver_for_Dynamic_Mixed_Volume_Computation(file,lp.all,true,lq.all,qsols,b);
end ts_drivdynl;
