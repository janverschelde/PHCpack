with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Solutions_Pool;

procedure ts_solpool is

  n : integer32 := 0;

begin
  new_line;
  put_line("Testing the operations in the solutions pool...");
  new_line;
  put("-> The size of the solutions pool : ");
  put(Solutions_Pool.Size,1); new_line;
  put("Give the number of solution lists : "); get(n);
  Solutions_Pool.Initialize(n);
  put("-> The size of the solutions pool : ");
  put(Solutions_Pool.Size,1); new_line;
end ts_solpool;
