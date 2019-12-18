with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Dense_Series;
with Standard_Dense_Series_Vectors;
with Standard_Random_Series;
with Standard_Series_Polynomials;
with Series_and_Polynomials_io;           use Series_and_Polynomials_io;
with Standard_Series_Poly_Functions;

procedure ts_speelser is

-- DESCRIPTION :
--   A vectorized version of the computation of Speelpenning products
--   for truncated power series.

  function Standard_Product
             ( dim,deg : in integer32 )
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
  end Standard_Product;

  procedure Standard_Evaluate_Product
              ( x : in Standard_Dense_Series_Vectors.Vector;
                d : in integer32 ) is

  -- DESCRIPTION :
  --   Inline computation of the coefficients of the product
  --   of the series in x, which are all of degree d.
  --   The range of x starts at 1 and ends at 2 or a larger value.

    acc,res : Standard_Complex_Vectors.Vector(0..d);

  begin
    if x'last = 2 then
      res(0) := x(1).cff(0)*x(2).cff(0);
      for k in 1..d loop
        res(k) := x(1).cff(0)*x(2).cff(k);
        for i in 1..k loop
          res(k) := res(k) + x(1).cff(i)*x(2).cff(k-i);
        end loop;
      end loop;
    else
      acc(0) := x(1).cff(0)*x(2).cff(0);
      for k in 1..d loop
        acc(k) := x(1).cff(0)*x(2).cff(k);
        for i in 1..k loop
          acc(k) := acc(k) + x(1).cff(i)*x(2).cff(k-i);
        end loop;
      end loop;
      for j in 3..x'last loop
        res(0) := acc(0)*x(j).cff(0);
        for k in 1..d loop
          res(k) := acc(0)*x(j).cff(k);
          for i in 1..k loop
            res(k) := res(k) + acc(i)*x(j).cff(k-i);
          end loop;
        end loop;
        if j < x'last
         then acc := res;
        end if;
      end loop;
    end if;
    put_line("The coefficients of the product : ");
    put_line(res);
  end Standard_Evaluate_Product;

  procedure Standard_Test ( dim,deg : in integer32 ) is

  -- DESCRIPTION :
  --   Generates random coefficient vectors for as many series as
  --   the value of dim, series truncated to degree deg.

    prd : constant Standard_Series_Polynomials.Poly
        := Standard_Product(dim,deg);
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
    put(y); new_line;
    Standard_Evaluate_Product(x,deg);
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
