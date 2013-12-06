with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;

procedure ts_mixvol is

-- DESCRIPTION :
--   Tests the mixed-volume computation.

  lp : Link_to_Poly_Sys;

begin
  new_line;
  put_line("Test on the mixed-volume computation.");
  new_line;
  get(lp);
  new_line;
  declare
    supports : constant Array_of_Lists(lp'range) := Create(lp.all);
    mv : constant natural32 := Mixed_Volume(lp'last,supports);
  begin
    put("The mixed volume : "); put(mv,1); new_line;
  end;
end ts_mixvol;
