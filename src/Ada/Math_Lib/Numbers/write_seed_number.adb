with Standard_Random_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;

procedure write_seed_number ( file : in file_type ) is
begin
  put(file,"Seed used in random number generators : ");
  put(file,Standard_Random_Numbers.Get_Seed,1);
  put_line(file,".");
end Write_Seed_Number;
