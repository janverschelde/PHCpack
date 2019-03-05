with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Random_Polynomials;
with Standard_Mixed_Residuals;

procedure ts_mixres is

-- DESCRIPTION :
--   Test the computation of mixed residuals.

  procedure Test ( dim : in integer32; deg,nbr : in natural32 ) is

  -- DESCRIPTION :
  --   Runs tests on vectors of dimension dim.

    rpt : constant Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim);
    apt : constant Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Mixed_Residuals.AbsVal(rpt);
    nvr : constant natural32 := natural32(dim);
    pol : Standard_Complex_Polynomials.Poly
        := Standard_Random_Polynomials.Random_Sparse_Poly(nvr,deg,nbr,0);
    abp : Standard_Complex_Polynomials.Poly
        := Standard_Mixed_Residuals.AbsVal(pol);
    res : double_float;

  begin
    put_line("a random point :"); put_line(rpt);
    put_line("its absolute values :"); put_line(apt);
    Symbol_Table.Init(nvr);
    put_line("a random polynomial :"); put(pol); new_line;
    put_line("its absolute version :"); put(abp); new_line;
    res := Standard_Mixed_Residuals.Residual(pol,abp,rpt);
    put("The mixed residual : "); put(res); new_line;
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension
  --   and then launches the tests.

    dim : integer32 := 0;
    deg,nbr : natural32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the largest degree : "); get(deg);
    put("Give the number of monomials : "); get(nbr);
    Test(dim,deg,nbr);
  end Main;

begin
  Main;
end ts_mixres;
