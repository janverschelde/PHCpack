with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Systems_Pool;

procedure ts_syspool is

  n : integer32 := 0;
  p : Link_to_Poly_Sys;

begin
  new_line;
  put_line("Testing the operations in the systems pool...");
  new_line;
  put("-> The size of the systems pool : ");
  put(Systems_Pool.Size,1); new_line;
  put("Give the number of systems : "); get(n);
  Systems_Pool.Initialize(n);
  put("-> The size of the systems pool : ");
  put(Systems_Pool.Size,1); new_line;
  new_line;
  put_line("Reading a polynomial system..."); get(p);
  for k in 1..n loop
    Systems_Pool.Create(k,p.all);
  end loop;
end ts_syspool;
