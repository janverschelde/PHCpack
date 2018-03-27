with Standard_Random_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;

procedure Write_Number_of_Tasks ( file : in file_type; nt : in natural32 ) is
begin
  if nt = 0 then
    put_line(file,"No multitasking was used in this run.");
  else
    put(file,"Number of tasks used in this run : ");
    put(file,nt,1);
    put_line(file,".");
  end if;
end Write_Number_of_Tasks;
