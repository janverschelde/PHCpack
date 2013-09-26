with text_io,integer_io;                use text_io,integer_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;

procedure ts_rewrite is

  procedure Test_Rewrite ( n : in natural ) is

    p : Vector(0..n) := Random_Vector(0,n);

  begin
    put_line("Generated a random vector...");
  end Test_Rewrite;

  procedure Main is

    n : natural;

  begin
    put("Give the degree : "); get(n);
    Test_Rewrite(n);
  end Main;

begin
  new_line;
  put_line("Testing the rewriting of a random polynomial.");
  new_line;
  Main;
end ts_rewrite;