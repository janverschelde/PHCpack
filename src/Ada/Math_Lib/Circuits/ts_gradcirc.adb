with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;
with Standard_Gradient_Circuits;         use Standard_Gradient_Circuits;

procedure ts_gradcirc is

-- DESCRIPTION :
--   Testing on the evaluation and differentiation with circuits.

  procedure Compare ( x : in Standard_Complex_Vectors.Vector;
                      y : in Standard_Complex_Vectors.Vector;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in x with the values in y.
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use Standard_Complex_Numbers;

    dff : Complex_Number;
    val : double_float;
 
  begin
    for i in x'range loop
      dff := x(i) - y(i);
      val := AbsVal(dff);
      if val > tol or output then
        put("difference at "); put(i,1); put_line(" :");
        put(x(i)); new_line;
        put(y(i)); new_line;
        put("error :"); put(val,3); new_line;
      end if;
    end loop;
  end Compare;

  procedure Standard_Circuit_Test
               ( p : in Standard_Complex_Polynomials.Poly;
                 c : in Standard_Gradient_Circuits.Circuit ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the polynomial p and the circuit
  --   at some random point to see if the values match.

    n : constant integer32 := integer32(Number_of_Variables(c));
    x : constant Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    ydx,zdx : Standard_Complex_Vectors.Vector(0..n);

  begin
    ydx(0) := Standard_Complex_Poly_Functions.Eval(p,x);
    for k in 1..n loop
      declare
        q : Standard_Complex_Polynomials.Poly
          := Standard_Complex_Polynomials.Diff(p,k);
      begin
        ydx(k) := Standard_Complex_Poly_Functions.Eval(q,x);
      end;
    end loop;
    zdx := EvalDiff(c,x);
    Compare(ydx,zdx,1.0E-8,true);
  end Standard_Circuit_Test;

  procedure Standard_Test ( p : in Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Tests the creation of a circuit to evaluate and differentiate p.

    c : Circuit := Create(p);
    m : constant natural32 := Number_of_Terms(c); 
    n : constant natural32 := Number_of_Variables(c); 

    use Standard_Natural_VecVecs;
  
  begin
    put("number of coefficients : "); put(m,1); new_line;
    put("number of variables : "); put(n,1); new_line;
    put_line("The coefficients : ");
    for k in 1..integer32(m) loop
      put(Coefficient(c,k)); new_line;
    end loop;
    put_line("The positions of the variables : ");
    for k in 1..integer32(m) loop
      put(Positions(c,k)); new_line;
    end loop;
    if Factors(c) = null then
      put_line("There are no common factors.");
    else
      put_line("The common factors : ");
      for k in 1..integer32(m) loop
        put(Factors(c,k)); new_line;
      end loop;
    end if;
    Standard_Circuit_Test(p,c);
  end Standard_Test;

  procedure Main is

    n : natural32 := 0;
    p : Poly;

  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give a polynomial : "); get(p);
    Standard_Test(p);
  end Main;

begin
  Main;
end ts_gradcirc;
