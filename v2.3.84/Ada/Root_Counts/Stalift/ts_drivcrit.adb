with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Driver_for_Criterion;

procedure ts_drivcrit is

-- DESCRIPTION :
--   This procedure calls the driver for the criterion.

  lp : Link_to_Poly_Sys;
  file : file_type;

begin
  new_line;
  put_line("Interactive testing of driver for the criterion.");
  new_line;
  get(lp);
  new_line;
  put_line("Reading the name of the output file.");
  Read_Name_and_Create_File(file);
  declare
    supports : Array_of_Lists(lp'range) := Create(lp.all);
  begin
    put_line(file,"The supports of the system : "); put(file,supports);
    Driver_for_Criterion(file,supports);
    put_line(file,"The reduced supports : "); put(file,supports);
  end;
end ts_drivcrit;
