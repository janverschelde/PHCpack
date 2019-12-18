with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Dense_Series;
--with Standard_Dense_Vector_Series2;
--with Standard_Dense_Vector_Series2_io;    use Standard_Dense_Vector_Series2_io;
with Standard_Dense_Series_Vectors;
with Standard_Random_Series;
with Standard_Series_Polynomials;
with Series_and_Polynomials_io;           use Series_and_Polynomials_io;
with Standard_Series_Poly_Functions;

procedure ts_speelser is

-- DESCRIPTION :
--   A vectorized version of the computation of Speelpenning products
--   for truncated power series.

  function Product ( dim,deg : in integer32 )
                   return Standard_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the product of the first dim variables,
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with double precision coefficients.

    res : Standard_Series_Polynomials.Poly;
    trm : Standard_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    trm.dg(1) := 1;
    trm.cf := Standard_Dense_Series.Create(1.0,deg);
    res := Standard_Series_Polynomials.Create(trm);
    Standard_Series_Polynomials.Clear(trm);
    for i in 2..dim loop
      declare
        trmtwo : Standard_Series_Polynomials.Term;
      begin
        trmtwo.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
        trmtwo.dg(i) := 1;
        trmtwo.cf := Standard_Dense_Series.Create(1.0,deg);
        Standard_Series_Polynomials.Mul(res,trmtwo);
        Standard_Series_Polynomials.Clear(trmtwo);
      end;
    end loop;
    return res;
  end Product;

  procedure Standard_Test ( dim,deg : in integer32 ) is

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg.

    prd : constant Standard_Series_Polynomials.Poly := Product(dim,deg);
    x : constant Standard_Dense_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series.Random_Series_Vector(1,dim,deg);
    y : Standard_Dense_Series.Series;

  begin
    put(dim,1); put(" random series of degree "); put(deg,1);
    put_line(" :"); put(x);
    put("The product of variables : ");
    put(prd); new_line;
    y := Standard_Series_Poly_Functions.Eval(prd,x);
    put_line("The value of the product at the random series :");
    put(y);
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree and the dimension.

    dim,deg : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    Standard_Test(dim,deg);
  end Main;

begin
  Main;
end ts_speelser;
