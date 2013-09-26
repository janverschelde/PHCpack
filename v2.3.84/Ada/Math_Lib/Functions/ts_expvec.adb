with text_io;                            use text_io;
with Symbol_Table;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;        use Standard_Integer_VecVecs_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Exponent_Vectors;                   use Exponent_Vectors;

procedure ts_expvec is

-- DESCRIPTION :
--   This routine provides basic testing routines for complex polynomials.

  procedure Read ( m : out natural32;
                   q : out Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Tests the input/output of a polynomial in several variables
  --   and with complex coefficients.

    n : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;

  begin
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    put_line("Your polynomial : "); put(p); new_line;
    Symbol_Table.Clear;
    q := p;
    m := n;
  end Read;

  function Coeff ( p : Standard_Complex_Polynomials.Poly;
                   e : Standard_Integer_VecVecs.VecVec )
                 return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(e'range);
    deg : Degrees := new Standard_Natural_Vectors.Vector(e(e'first)'range);

  begin
    for i in e'range loop
      for j in e(i)'range loop
        deg(j) := natural32(e(i)(j));
      end loop;
      res(i) := Coeff(p,deg);
    end loop;
    Clear(deg);
    return res;
  end Coeff;

  procedure Test_Eval ( n : in natural32;
                        p : in Standard_Complex_Polynomials.Poly ) is

    ev : constant Standard_Integer_VecVecs.VecVec := Create(p);
    cf : constant Standard_Complex_Vectors.Vector(ev'range) := Coeff(p,ev);
    x : Standard_Complex_Vectors.Vector(1..integer32(n));
    y1,y2 : Complex_Number;
    ans : character;

  begin
    put_line("The exponent vectors : "); put(ev);
    put_line("The coefficients : ");
    for i in cf'range loop
      put(cf(i)); new_line;
    end loop;
    loop
      put("Give "); put(n,1); put_line(" complex numbers : "); 
      for i in x'range loop
        get(x(i));
      end loop;
      y1 := Eval(p,x);
      y2 := Eval(ev,cf,x);
      put("Eval poly p(x) : "); put(y1); new_line;
      put("Eval cf*ev^x   : "); put(y2); new_line;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Eval;

  procedure Main is

    n : natural32;
    p : Poly;

  begin
    new_line;
    put_line("Interactive testing of the operations on exponent vectors.");
    new_line;
    Read(n,p);
    Test_Eval(n,p);
  end Main;

begin
  Main;
end ts_expvec;
