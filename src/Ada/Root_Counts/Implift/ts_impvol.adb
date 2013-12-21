with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Trees_of_Vectors;                   use Trees_of_Vectors;
with Trees_of_Vectors_io;                use Trees_of_Vectors_io;
with Volumes;                            use Volumes;

procedure ts_impvol is

-- DESCRIPTION :
--   Testing mixed-volume computation by implicit lifting.

  lp : Link_to_Poly_Sys;

begin
  new_line;
  put_line("Testing mixed-volume computation by implicit lifting.");
  new_line;
  get(lp);
  declare
    supports : constant Array_of_Lists(lp'range) := Create(lp.all);
    n : constant natural32 := natural32(lp'last);
    tv : Tree_of_Vectors;
    mv : natural32;
  begin
    put_line("The supports of the system : "); put(supports);
    Mixed_Volume(n,supports,tv,mv);
    put("The mixed volume : "); put(mv,1); new_line;
    put_line("The tree of vectors : "); put(tv);
  end;
end ts_impvol;
