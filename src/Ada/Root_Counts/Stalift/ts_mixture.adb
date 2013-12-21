with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Symbol_Table;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;

procedure ts_mixture is

  ans : character;

  procedure Compute_Mixture ( p : in Poly_Sys ) is

    supports : Array_of_Lists(p'range) := Create(p);
    mix,perms : Link_to_Vector;

  begin
    Compute_Mixture(supports,mix,perms);
    put("Type of mixture : "); put(mix.all); new_line;
  end Compute_Mixture;

  procedure Main_Test is

    file : file_type;
    lp : Link_to_Poly_Sys;
  
  begin
    put_line("Reading the name where the file is.");
    Read_Name_and_Open_File(file);
    get(file,lp);
    Compute_Mixture(lp.all);
  end Main_Test;

begin
  new_line;
  put_line("Testing the computation of the type of mixture.");
  new_line;
  loop
    Main_Test;
    put("Do you want more tests ? (y/n) "); get(ans);
    exit when ans /= 'y';
    skip_line;
    Symbol_Table.Clear;
  end loop;
end ts_mixture;
