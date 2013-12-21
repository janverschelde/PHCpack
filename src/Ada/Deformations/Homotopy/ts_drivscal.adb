with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Drivers_for_Scaling;                use Drivers_for_Scaling;

procedure ts_drivscal is

-- DESCRIPTION :
--   Calls the driver routine to scaling a polynomial system.

  file : file_type;
  lp : Link_to_Poly_Sys;
  scvc : Link_to_Vector;
  bas : natural32 := 2;

begin
  new_line;
  put_line("Test on scaling polynomial systems.");
  new_line;
  get(lp);
  new_line;
  put_line("Reading the name of the output file.");
  Read_Name_and_Create_File(file);
  Driver_for_Scaling(file,lp.all,bas,scvc);
end ts_drivscal;
