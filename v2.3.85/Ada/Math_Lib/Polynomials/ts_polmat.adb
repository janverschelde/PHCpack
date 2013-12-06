with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Symbol_Table;
with Standard_Complex_Poly_Matrices;     use Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;

procedure ts_polmat is

  procedure Test_io ( n,m : in integer32  ) is

    pm : Matrix(1..n,1..m);

  begin
    put("Give "); put(n*m,1); put_line(" polynomials : ");
    get(pm);
    put_line("The matrix of polynomials : ");
    put(pm);
  end Test_io;

  procedure Main is

    n,m : integer32 := 0;
    nb : natural32 := 0;

  begin
    put("Give number of rows : "); get(n);
    put("Give number of columns : "); get(m);
    put("Give number of symbols : "); get(nb);
    Symbol_Table.Init(nb);
    Test_io(n,m);
  end Main;

begin
  new_line;
  put_line("Testing input/output of matrices of polynomials.");
  new_line;
  Main;
end ts_polmat;
