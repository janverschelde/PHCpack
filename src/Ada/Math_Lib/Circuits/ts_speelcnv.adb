with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;
with Standard_Random_Vectors;
with Exponent_Indices;
with Standard_Dense_Series;
with Standard_Dense_Series_Vectors;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Functions;
with Standard_Random_Series;
with Series_and_Polynomials_io;           use Series_and_Polynomials_io;

procedure ts_speelcnv is

-- DESCRIPTION :
--   Tests the evaluation of the gradient of a polynomial in several variables,
--   in a power series of some fixed degree.

  function Standard_Polynomial
             ( dim,deg : in integer32;
               xps : Standard_Integer_VecVecs.VecVec )
             return Standard_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the polynomial in dim variables, with exponents in xps. 
  --   as a polynomial where the coefficients are truncated power series,
  --   truncated to degree deg, with standard double precision coefficients.

    res : Standard_Series_Polynomials.Poly
        := Standard_Series_Polynomials.Null_Poly;
    trm : Standard_Series_Polynomials.Term;

  begin
    for k in xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
      for i in xps(k)'range loop
        trm.dg(xps(k)(i)) := 1;
      end loop;
      trm.cf := Standard_Dense_Series.Create(1.0,deg);
      Standard_Series_Polynomials.Add(res,trm);
      Standard_Series_Polynomials.Clear(trm);
    end loop;
    return res;
  end Standard_Polynomial;

  function Standard_Gradient
             ( p : Standard_Series_Polynomials.Poly;
               x : Standard_Dense_Series_Vectors.Vector )
             return Standard_Dense_Series_Vectors.Vector is

  -- DESCRIPTION :
  --   Evaluates the gradient of p at x, for testing purposes,
  --   in standard double precision.

    res : Standard_Dense_Series_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      declare
        dpk : Standard_Series_Polynomials.Poly
            := Standard_Series_Polynomials.Diff(p,k);
        val : constant Standard_Dense_Series.Series
            := Standard_Series_Poly_Functions.Eval(dpk,x);
      begin
        Standard_Series_Polynomials.Clear(dpk);
        res(k) := val;
      end;
    end loop;
    return res;
  end Standard_Gradient;

  function Random_Exponents
             ( dim,nbr : integer32 ) 
             return Standard_Integer_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Generates as many random 0/1 vectors of dimension dim
  --   as the value of nrb.

    res : Standard_Integer_VecVecs.VecVec(1..nbr);

  begin
    for i in 1..nbr loop
      declare
        xp : Standard_Integer_Vectors.Vector(1..dim)
           := Standard_Random_Vectors.Random_Vector(1,dim,0,1);
      begin
        res(i) := new Standard_Integer_Vectors.Vector'(xp);
      end;
    end loop;
    return res;
  end Random_Exponents;

  procedure Standard_Test ( dim,deg,nbr : in integer32 ) is

    xps : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Exponents(dim,nbr);
    idx : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Exponent_Indices.Exponent_Index(xps);
    pol : Standard_Series_Polynomials.Poly
        := Standard_Polynomial(dim,deg,idx);
    x : constant Standard_Dense_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series.Random_Series_Vector(1,dim,deg);
    y : Standard_Dense_Series.Series;
    grad : Standard_Dense_Series_Vectors.Vector(1..dim);

  begin
    put_line("Some random exponents :");
    Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :");
    Standard_Integer_VecVecs_io.put(idx);
    put_line("the polynomial :"); put(pol);
    y := Standard_Series_Poly_Functions.Eval(pol,x);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    grad := Standard_Gradient(pol,x);
    put_line("the last derivative :"); put(grad(dim)); new_line;
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree, the dimension, and
  --   the number of monomials.  Then runs the tests.

    dim,deg,nbr : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    put("Give the number of monomials : "); get(nbr);
    new_line;
    Standard_Test(dim,deg,nbr);
  end Main;

begin
  Main;
end ts_speelcnv;
