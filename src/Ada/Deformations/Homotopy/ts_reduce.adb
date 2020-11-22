with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Main_Reduction;

procedure ts_reduce is

-- DESCRIPTION :
--   Calls the driver routine to reducing a polynomial system.

  file : file_type;
  lp : Link_to_Poly_Sys;
  d : natural32;

begin
  new_line;
  put_line("Test on reducing polynomial systems.");
  new_line;
  get(lp);
  new_line;
  put_line("Reading the name of the output file.");
  Read_Name_and_Create_File(file);
  Main_Reduction.Reduce(file,lp.all,d,false);
end ts_reduce;
