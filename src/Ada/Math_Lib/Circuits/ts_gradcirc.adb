with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Vectors;
with Multprec_Random_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;
with Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Polynomials;       use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Multprec_Complex_Polynomials;       use Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Functions;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;
with Standard_Gradient_Circuits;
with Standard_Jacobian_Circuits;
with DoblDobl_Gradient_Circuits;
with DoblDobl_Jacobian_Circuits;
with QuadDobl_Gradient_Circuits;
with QuadDobl_Jacobian_Circuits;
with Multprec_Gradient_Circuits;
with Multprec_Jacobian_Circuits;

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

  procedure Compare ( x : in Standard_Complex_Matrices.Matrix;
                      y : in Standard_Complex_Matrices.Matrix;
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
    for i in x'range(1) loop
      for j in x'range(2) loop
        dff := x(i,j) - y(i,j);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference at ("); put(i,1); put(",");
          put(j,1); put_line(") :");
          put(x(i,j)); new_line;
          put(y(i,j)); new_line;
          put("error :"); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( x : in Standard_Complex_Matrices.Matrix;
                      y : in Standard_Complex_VecVecs.VecVec;
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
    for i in x'range(1) loop
      for j in x'range(2) loop
        dff := x(i,j) - y(i)(j);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference at ("); put(i,1); put(",");
          put(j,1); put_line(") :");
          put(x(i,j)); new_line;
          put(y(i)(j)); new_line;
          put("error :"); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( x : in DoblDobl_Complex_Vectors.Vector;
                      y : in DoblDobl_Complex_Vectors.Vector;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in x with the values in y.
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use DoblDobl_Complex_Numbers;

    dff : Complex_Number;
    val : double_double;
 
  begin
    for i in x'range loop
      dff := x(i) - y(i);
      val := AbsVal(dff);
      if val > tol or output then
        put("difference at "); put(i,1); put_line(" :");
        put(x(i)); new_line;
        put(y(i)); new_line;
        put("error : "); put(val,3); new_line;
      end if;
    end loop;
  end Compare;

  procedure Compare ( x : in DoblDobl_Complex_Matrices.Matrix;
                      y : in DoblDobl_Complex_Matrices.Matrix;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in x with the values in y.
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use DoblDobl_Complex_Numbers;

    dff : Complex_Number;
    val : double_double;
 
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        dff := x(i,j) - y(i,j);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference at ("); put(i,1); put(",");
          put(j,1); put_line(") :");
          put(x(i,j)); new_line;
          put(y(i,j)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( x : in DoblDobl_Complex_Matrices.Matrix;
                      y : in DoblDobl_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in x with the values in y.
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use DoblDobl_Complex_Numbers;

    dff : Complex_Number;
    val : double_double;
 
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        dff := x(i,j) - y(i)(j);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference at ("); put(i,1); put(",");
          put(j,1); put_line(") :");
          put(x(i,j)); new_line;
          put(y(i)(j)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( x : in QuadDobl_Complex_Vectors.Vector;
                      y : in QuadDobl_Complex_Vectors.Vector;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in x with the values in y.
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use QuadDobl_Complex_Numbers;

    dff : Complex_Number;
    val : quad_double;
 
  begin
    for i in x'range loop
      dff := x(i) - y(i);
      val := AbsVal(dff);
      if val > tol or output then
        put("difference at "); put(i,1); put_line(" :");
        put(x(i)); new_line;
        put(y(i)); new_line;
        put("error : "); put(val,3); new_line;
      end if;
    end loop;
  end Compare;

  procedure Compare ( x : in QuadDobl_Complex_Matrices.Matrix;
                      y : in QuadDobl_Complex_Matrices.Matrix;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in x with the values in y.
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use QuadDobl_Complex_Numbers;

    dff : Complex_Number;
    val : quad_double;
 
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        dff := x(i,j) - y(i,j);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference at ("); put(i,1); put(",");
          put(j,1); put_line(") :");
          put(x(i,j)); new_line;
          put(y(i,j)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( x : in QuadDobl_Complex_Matrices.Matrix;
                      y : in QuadDobl_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in x with the values in y.
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use QuadDobl_Complex_Numbers;

    dff : Complex_Number;
    val : quad_double;
 
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        dff := x(i,j) - y(i)(j);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference at ("); put(i,1); put(",");
          put(j,1); put_line(") :");
          put(x(i,j)); new_line;
          put(y(i)(j)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( x : in Multprec_Complex_Vectors.Vector;
                      y : in Multprec_Complex_Vectors.Vector;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in x with the values in y.
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use Multprec_Complex_Numbers;

    dff : Complex_Number;
    val : Floating_Number;
 
  begin
    for i in x'range loop
      dff := x(i) - y(i);
      val := AbsVal(dff);
      if val > tol or output then
        put("difference at "); put(i,1); put_line(" :");
        put(x(i)); new_line;
        put(y(i)); new_line;
        put("error : "); put(val,3); new_line;
      end if;
    end loop;
  end Compare;

  procedure Compare ( x : in Multprec_Complex_Matrices.Matrix;
                      y : in Multprec_Complex_Matrices.Matrix;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in x with the values in y.
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use Multprec_Complex_Numbers;

    dff : Complex_Number;
    val : Floating_Number;
 
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        dff := x(i,j) - y(i,j);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference at ("); put(i,1); put(",");
          put(j,1); put_line(") :");
          put(x(i,j)); new_line;
          put(y(i,j)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
        Clear(dff); Clear(val);
      end loop;
    end loop;
  end Compare;

  procedure Compare ( x : in Multprec_Complex_Matrices.Matrix;
                      y : in Multprec_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in x with the values in y.
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use Multprec_Complex_Numbers;

    dff : Complex_Number;
    val : Floating_Number;
 
  begin
    for i in x'range(1) loop
      for j in x'range(2) loop
        dff := x(i,j) - y(i)(j);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference at ("); put(i,1); put(",");
          put(j,1); put_line(") :");
          put(x(i,j)); new_line;
          put(y(i)(j)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
        Clear(dff); Clear(val);
      end loop;
    end loop;
  end Compare;

  procedure Standard_Circuit_Test
               ( p : in Standard_Complex_Polynomials.Poly;
                 c : in Standard_Gradient_Circuits.Circuit ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the polynomial p and the circuit
  --   at some random point to see if the values match.

    use Standard_Gradient_Circuits;

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
    zdx := Standard_Gradient_Circuits.EvalDiff(c,x);
    Compare(ydx,zdx,1.0E-8,true);
  end Standard_Circuit_Test;

  procedure DoblDobl_Circuit_Test
               ( p : in DoblDobl_Complex_Polynomials.Poly;
                 c : in DoblDobl_Gradient_Circuits.Circuit ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the polynomial p and the circuit
  --   at some random point to see if the values match.

    use DoblDobl_Gradient_Circuits;

    n : constant integer32 := integer32(Number_of_Variables(c));
    x : constant DoblDobl_Complex_Vectors.Vector(1..n)
      := DoblDobl_Random_Vectors.Random_Vector(1,n);
    ydx,zdx : DoblDobl_Complex_Vectors.Vector(0..n);

  begin
    ydx(0) := DoblDobl_Complex_Poly_Functions.Eval(p,x);
    for k in 1..n loop
      declare
        q : DoblDobl_Complex_Polynomials.Poly
          := DoblDobl_Complex_Polynomials.Diff(p,k);
      begin
        ydx(k) := DoblDobl_Complex_Poly_Functions.Eval(q,x);
      end;
    end loop;
    zdx := DoblDobl_Gradient_Circuits.EvalDiff(c,x);
    Compare(ydx,zdx,1.0E-8,true);
  end DoblDobl_Circuit_Test;

  procedure QuadDobl_Circuit_Test
               ( p : in QuadDobl_Complex_Polynomials.Poly;
                 c : in QuadDobl_Gradient_Circuits.Circuit ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the polynomial p and the circuit
  --   at some random point to see if the values match.

    use QuadDobl_Gradient_Circuits;

    n : constant integer32 := integer32(Number_of_Variables(c));
    x : constant QuadDobl_Complex_Vectors.Vector(1..n)
      := QuadDobl_Random_Vectors.Random_Vector(1,n);
    ydx,zdx : QuadDobl_Complex_Vectors.Vector(0..n);

  begin
    ydx(0) := QuadDobl_Complex_Poly_Functions.Eval(p,x);
    for k in 1..n loop
      declare
        q : QuadDobl_Complex_Polynomials.Poly
          := QuadDobl_Complex_Polynomials.Diff(p,k);
      begin
        ydx(k) := QuadDobl_Complex_Poly_Functions.Eval(q,x);
      end;
    end loop;
    zdx := QuadDobl_Gradient_Circuits.EvalDiff(c,x);
    Compare(ydx,zdx,1.0E-8,true);
  end QuadDobl_Circuit_Test;

  procedure Multprec_Circuit_Test
               ( p : in Multprec_Complex_Polynomials.Poly;
                 c : in Multprec_Gradient_Circuits.Circuit;
                 size : in natural32 ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the polynomial p and the circuit
  --   at some random point to see if the values match.

    use Multprec_Gradient_Circuits;

    n : constant integer32 := integer32(Number_of_Variables(c));
    x : constant Multprec_Complex_Vectors.Vector(1..n)
      := Multprec_Random_Vectors.Random_Vector(1,n,size);
    ydx,zdx : Multprec_Complex_Vectors.Vector(0..n);

  begin
    ydx(0) := Multprec_Complex_Poly_Functions.Eval(p,x);
    for k in 1..n loop
      declare
        q : Multprec_Complex_Polynomials.Poly
          := Multprec_Complex_Polynomials.Diff(p,k);
      begin
        ydx(k) := Multprec_Complex_Poly_Functions.Eval(q,x);
      end;
    end loop;
    zdx := Multprec_Gradient_Circuits.EvalDiff(c,x);
    Compare(ydx,zdx,1.0E-8,true);
  end Multprec_Circuit_Test;

  procedure Standard_Test ( p : in Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Tests the creation of a circuit to evaluate and differentiate p,
  --   for a polynomial p with coefficients in standard double precision.

    use Standard_Natural_VecVecs;
    use Standard_Gradient_Circuits;

    c : Circuit := Create(p);
    m : constant natural32 := Number_of_Terms(c); 
    n : constant natural32 := Number_of_Variables(c); 
  
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

  procedure DoblDobl_Test ( p : in DoblDobl_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Tests the creation of a circuit to evaluate and differentiate p,
  --   for a polynomial p with coefficients in double double precision.

    use Standard_Natural_VecVecs;
    use DoblDobl_Gradient_Circuits;

    c : Circuit := Create(p);
    m : constant natural32 := Number_of_Terms(c); 
    n : constant natural32 := Number_of_Variables(c); 
  
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
    DoblDobl_Circuit_Test(p,c);
  end DoblDobl_Test;

  procedure QuadDobl_Test ( p : in QuadDobl_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Tests the creation of a circuit to evaluate and differentiate p,
  --   for a polynomial p with coefficients in quad double precision.

    use Standard_Natural_VecVecs;
    use QuadDobl_Gradient_Circuits;

    c : Circuit := Create(p);
    m : constant natural32 := Number_of_Terms(c); 
    n : constant natural32 := Number_of_Variables(c); 
  
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
    QuadDobl_Circuit_Test(p,c);
  end QuadDobl_Test;

  procedure Multprec_Test ( p : in Multprec_Complex_Polynomials.Poly;
                            size : in natural32 ) is

  -- DESCRIPTION :
  --   Tests the creation of a circuit to evaluate and differentiate p,
  --   for a polynomial p with coefficients of the given size.

    use Standard_Natural_VecVecs;
    use Multprec_Gradient_Circuits;

    c : Circuit := Create(p);
    m : constant natural32 := Number_of_Terms(c); 
    n : constant natural32 := Number_of_Variables(c); 
  
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
    Multprec_Circuit_Test(p,c,size);
  end Multprec_Test;

  procedure Standard_Test ( n : natural32 ) is

  -- DESCRPITION :
  --  Prompts the user for a polynomial in n variables,
  --  with coefficients in standard double precision.

    p : Standard_Complex_Polynomials.Poly;

  begin
    put("Give a polynomial : "); get(p);
    Standard_Test(p);
  end Standard_Test;

  procedure DoblDobl_Test ( n : natural32 ) is

  -- DESCRPITION :
  --  Prompts the user for a polynomial in n variables,
  --  with coefficients in double double precision.

    p : DoblDobl_Complex_Polynomials.Poly;

  begin
    put("Give a polynomial : "); get(p);
    DoblDobl_Test(p);
  end DoblDobl_Test;

  procedure QuadDobl_Test ( n : natural32 ) is

  -- DESCRPITION :
  --  Prompts the user for a polynomial in n variables,
  --  with coefficients in double double precision.

    p : QuadDobl_Complex_Polynomials.Poly;

  begin
    put("Give a polynomial : "); get(p);
    QuadDobl_Test(p);
  end QuadDobl_Test;

  procedure Multprec_Test ( n,size : natural32 ) is

  -- DESCRPITION :
  --  Prompts the user for a polynomial in n variables,
  --  with coefficients of the given size.

    p : Multprec_Complex_Polynomials.Poly;

  begin
    put("Give a polynomial : "); get(p);
    Multprec_Test(p,size);
  end Multprec_Test;

  procedure Standard_Jacobian_Test 
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                c : in Standard_Jacobian_Circuits.Circuit ) is

  -- DESCRIPTION :
  --   Tests the evaluation at a random point using the circuit c
  --   defined by a polynomial system p,
  --   with comparisons to the straightforward evaluations.

    use Standard_Jacobian_Circuits;

    nm : constant integer32 := integer32(Number_of_Monomials(c));
    nq : constant integer32 := integer32(Number_of_Polynomials(c));
    nv : constant integer32 := integer32(Number_of_Variables(c));
    x : constant Standard_Complex_Vectors.Vector(1..nv)
      := Standard_Random_Vectors.Random_Vector(1,nv);
    px : constant Standard_Complex_Vectors.Vector(1..nq)
       := Standard_Complex_Poly_SysFun.Eval(p,x);
    jm : Standard_Complex_Jaco_Matrices.Jaco_Mat(1..nq,1..nv)
       := Standard_Complex_Jaco_Matrices.Create(p);
    Ap : Standard_Complex_Matrices.Matrix(jm'range(1),jm'range(2))
       := Standard_Complex_Jaco_Matrices.Eval(jm,x);
    wrk : Standard_Complex_VecVecs.VecVec(1..nm) := WorkSpace(c);
    y : Standard_Complex_Vectors.Vector(1..nq);
    A : Standard_Complex_Matrices.Matrix(1..nq,1..nv);
    B : Standard_Complex_VecVecs.VecVec(1..nv);

  begin
    EvalDiff(c,x,wrk,y,A);
    for i in 1..nv loop
      B(i) := new Standard_Complex_Vectors.Vector(1..nq);
    end loop;
    EvalDiff(c,x,wrk,y,B);
    put_line("Comparing the evaluation :");
    Compare(px,y,1.0E-8,true);
    put_line("Comparing the differentiation :");
    Compare(Ap,A,1.0E-8,true);
    put_line("Comparing the vectors of vectors :");
    Compare(Ap,B,1.0E-8,true);
    Standard_Complex_Jaco_Matrices.Clear(jm);
  end Standard_Jacobian_Test;

  procedure DoblDobl_Jacobian_Test 
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                c : in DoblDobl_Jacobian_Circuits.Circuit ) is

  -- DESCRIPTION :
  --   Tests the evaluation at a random point using the circuit c
  --   defined by a polynomial system p,
  --   with comparisons to the straightforward evaluations.

    use DoblDobl_Jacobian_Circuits;

    nm : constant integer32 := integer32(Number_of_Monomials(c));
    nq : constant integer32 := integer32(Number_of_Polynomials(c));
    nv : constant integer32 := integer32(Number_of_Variables(c));
    x : constant DoblDobl_Complex_Vectors.Vector(1..nv)
      := DoblDobl_Random_Vectors.Random_Vector(1,nv);
    px : constant DoblDobl_Complex_Vectors.Vector(1..nq)
       := DoblDobl_Complex_Poly_SysFun.Eval(p,x);
    jm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(1..nq,1..nv)
       := DoblDobl_Complex_Jaco_Matrices.Create(p);
    Ap : DoblDobl_Complex_Matrices.Matrix(jm'range(1),jm'range(2))
       := DoblDobl_Complex_Jaco_Matrices.Eval(jm,x);
    wrk : DoblDobl_Complex_VecVecs.VecVec(1..nm) := WorkSpace(c);
    y : DoblDobl_Complex_Vectors.Vector(1..nq);
    A : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nv);
    B : DoblDobl_Complex_VecVecs.VecVec(1..nv);

  begin
    EvalDiff(c,x,wrk,y,A);
    for i in 1..nv loop
      B(i) := new DoblDobl_Complex_Vectors.Vector(1..nq);
    end loop;
    EvalDiff(c,x,wrk,y,B);
    put_line("Comparing the evaluation :");
    Compare(px,y,1.0E-8,true);
    put_line("Comparing the differentiation :");
    Compare(Ap,A,1.0E-8,true);
    put_line("Comparing the vectors of vectors :");
    Compare(Ap,B,1.0E-8,true);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
  end DoblDobl_Jacobian_Test;

  procedure QuadDobl_Jacobian_Test 
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                c : in QuadDobl_Jacobian_Circuits.Circuit ) is

  -- DESCRIPTION :
  --   Tests the evaluation at a random point using the circuit c
  --   defined by a polynomial system p,
  --   with comparisons to the straightforward evaluations.

    use QuadDobl_Jacobian_Circuits;

    nm : constant integer32 := integer32(Number_of_Monomials(c));
    nq : constant integer32 := integer32(Number_of_Polynomials(c));
    nv : constant integer32 := integer32(Number_of_Variables(c));
    x : constant QuadDobl_Complex_Vectors.Vector(1..nv)
      := QuadDobl_Random_Vectors.Random_Vector(1,nv);
    px : constant QuadDobl_Complex_Vectors.Vector(1..nq)
       := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    jm : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(1..nq,1..nv)
       := QuadDobl_Complex_Jaco_Matrices.Create(p);
    Ap : QuadDobl_Complex_Matrices.Matrix(jm'range(1),jm'range(2))
       := QuadDobl_Complex_Jaco_Matrices.Eval(jm,x);
    wrk : QuadDobl_Complex_VecVecs.VecVec(1..nm) := WorkSpace(c);
    y : QuadDobl_Complex_Vectors.Vector(1..nq);
    A : QuadDobl_Complex_Matrices.Matrix(1..nq,1..nv);
    B : QuadDobl_Complex_VecVecs.VecVec(1..nv);

  begin
    EvalDiff(c,x,wrk,y,A);
    for i in 1..nv loop
      B(i) := new QuadDobl_Complex_Vectors.Vector(1..nq);
    end loop;
    EvalDiff(c,x,wrk,y,B);
    put_line("Comparing the evaluation :");
    Compare(px,y,1.0E-8,true);
    put_line("Comparing the differentiation :");
    Compare(Ap,A,1.0E-8,true);
    put_line("Comparing the vectors of vectors :");
    Compare(Ap,B,1.0E-8,true);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
  end QuadDobl_Jacobian_Test;

  procedure Multprec_Jacobian_Test 
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                c : in Multprec_Jacobian_Circuits.Circuit;
                size : in natural32 ) is

  -- DESCRIPTION :
  --   Tests the evaluation at a random point using the circuit c
  --   defined by a polynomial system p,
  --   with comparisons to the straightforward evaluations,
  --   in arbitrary multiprecision with numbers of the given size.

    use Multprec_Jacobian_Circuits;

    nm : constant integer32 := integer32(Number_of_Monomials(c));
    nq : constant integer32 := integer32(Number_of_Polynomials(c));
    nv : constant integer32 := integer32(Number_of_Variables(c));
    x : Multprec_Complex_Vectors.Vector(1..nv)
      := Multprec_Random_Vectors.Random_Vector(1,nv,size);
    px : Multprec_Complex_Vectors.Vector(1..nq)
       := Multprec_Complex_Poly_SysFun.Eval(p,x);
    jm : Multprec_Complex_Jaco_Matrices.Jaco_Mat(1..nq,1..nv)
       := Multprec_Complex_Jaco_Matrices.Create(p);
    Ap : Multprec_Complex_Matrices.Matrix(jm'range(1),jm'range(2))
       := Multprec_Complex_Jaco_Matrices.Eval(jm,x);
    wrk : Multprec_Complex_VecVecs.VecVec(1..nm) := WorkSpace(c);
    y : Multprec_Complex_Vectors.Vector(1..nq);
    A : Multprec_Complex_Matrices.Matrix(1..nq,1..nv);
    B : Multprec_Complex_VecVecs.VecVec(1..nv);

  begin
    EvalDiff(c,x,wrk,y,A);
    for i in 1..nv loop
      B(i) := new Multprec_Complex_Vectors.Vector(1..nq);
    end loop;
    EvalDiff(c,x,wrk,y,B);
    put_line("Comparing the evaluation :");
    Compare(px,y,1.0E-8,true);
    put_line("Comparing the differentiation :");
    Compare(Ap,A,1.0E-8,true);
    put_line("Comparing the vectors of vectors :");
    Compare(Ap,B,1.0E-8,true);
    Multprec_Complex_Jaco_Matrices.Clear(jm);
  end Multprec_Jacobian_Test;

  procedure Standard_System_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then
  --   tests the operations on the circuit representing the system,
  --   in standard double precision.

    use Standard_Jacobian_Circuits;
    use Standard_Natural_Vectors;

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    comfac : Standard_Natural_Vectors.Link_to_Vector;
    c : Circuit;
    nq : integer32;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(p);
    c := Create(p.all);
    put("-> number of equations : ");
    put(Number_of_Polynomials(c),1); new_line;
    put("-> number of variables : ");
    put(Number_of_Variables(c),1); new_line;
    put("-> number of monomials : ");
    put(Number_of_Monomials(c),1); new_line;
    nq := integer32(Number_of_Polynomials(c));
    for k in 1..nq loop
      put("-> polynomial "); put(k,1); put(" has ");
      put(Number_of_Terms(c,k),1); put_line(" terms");
    end loop;
    for k in 1..nq loop
      put("-> polynomial "); put(k,1); put_line(" : ");
      for i in 1..integer32(Number_of_Terms(c,k)) loop
        put(Coefficient(c,k,i));
        put(Product(c,k,i).all);
        comfac := Factor(c,k,i);
        if comfac /= null then
          put(" +"); put(comfac.all);
        end if;
        new_line;
      end loop;
    end loop;
    Standard_Jacobian_Test(p.all,c);
  end Standard_System_Test;

  procedure DoblDobl_System_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then
  --   tests the operations on the circuit representing the system,
  --   in double double precision.

    use DoblDobl_Jacobian_Circuits;
    use Standard_Natural_Vectors;

    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    comfac : Standard_Natural_Vectors.Link_to_Vector;
    c : Circuit;
    nq : integer32;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(p);
    c := Create(p.all);
    put("-> number of equations : ");
    put(Number_of_Polynomials(c),1); new_line;
    put("-> number of variables : ");
    put(Number_of_Variables(c),1); new_line;
    put("-> number of monomials : ");
    put(Number_of_Monomials(c),1); new_line;
    nq := integer32(Number_of_Polynomials(c));
    for k in 1..nq loop
      put("-> polynomial "); put(k,1); put(" has ");
      put(Number_of_Terms(c,k),1); put_line(" terms");
    end loop;
    for k in 1..nq loop
      put("-> polynomial "); put(k,1); put_line(" : ");
      for i in 1..integer32(Number_of_Terms(c,k)) loop
        put(Coefficient(c,k,i));
        put(Product(c,k,i).all);
        comfac := Factor(c,k,i);
        if comfac /= null then
          put(" +"); put(comfac.all);
        end if;
        new_line;
      end loop;
    end loop;
    DoblDobl_Jacobian_Test(p.all,c);
  end DoblDobl_System_Test;

  procedure QuadDobl_System_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then
  --   tests the operations on the circuit representing the system,
  --   in quad double precision.

    use QuadDobl_Jacobian_Circuits;
    use Standard_Natural_Vectors;

    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    comfac : Standard_Natural_Vectors.Link_to_Vector;
    c : Circuit;
    nq : integer32;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(p);
    c := Create(p.all);
    put("-> number of equations : ");
    put(Number_of_Polynomials(c),1); new_line;
    put("-> number of variables : ");
    put(Number_of_Variables(c),1); new_line;
    put("-> number of monomials : ");
    put(Number_of_Monomials(c),1); new_line;
    nq := integer32(Number_of_Polynomials(c));
    for k in 1..nq loop
      put("-> polynomial "); put(k,1); put(" has ");
      put(Number_of_Terms(c,k),1); put_line(" terms");
    end loop;
    for k in 1..nq loop
      put("-> polynomial "); put(k,1); put_line(" : ");
      for i in 1..integer32(Number_of_Terms(c,k)) loop
        put(Coefficient(c,k,i));
        put(Product(c,k,i).all);
        comfac := Factor(c,k,i);
        if comfac /= null then
          put(" +"); put(comfac.all);
        end if;
        new_line;
      end loop;
    end loop;
    QuadDobl_Jacobian_Test(p.all,c);
  end QuadDobl_System_Test;

  procedure Multprec_System_Test ( size : in natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then
  --   tests the operations on the circuit representing the system,
  --   in arbitrary multiprecision, with numbers of the given size.

    use Multprec_Jacobian_Circuits;
    use Standard_Natural_Vectors;

    p : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    comfac : Standard_Natural_Vectors.Link_to_Vector;
    c : Circuit;
    nq : integer32;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(p);
    c := Create(p.all);
    put("-> number of equations : ");
    put(Number_of_Polynomials(c),1); new_line;
    put("-> number of variables : ");
    put(Number_of_Variables(c),1); new_line;
    put("-> number of monomials : ");
    put(Number_of_Monomials(c),1); new_line;
    nq := integer32(Number_of_Polynomials(c));
    for k in 1..nq loop
      put("-> polynomial "); put(k,1); put(" has ");
      put(Number_of_Terms(c,k),1); put_line(" terms");
    end loop;
    for k in 1..nq loop
      put("-> polynomial "); put(k,1); put_line(" : ");
      for i in 1..integer32(Number_of_Terms(c,k)) loop
        put(Coefficient(c,k,i));
        put(Product(c,k,i).all);
        comfac := Factor(c,k,i);
        if comfac /= null then
          put(" +"); put(comfac.all);
        end if;
        new_line;
      end loop;
    end loop;
    Multprec_Jacobian_Test(p.all,c,size);
  end Multprec_System_Test;

  function Ask_for_Size return natural32 is

  -- DESCRIPTION :
  --   Asks the user for the number of decimal places
  --   and returns the corresponding size of the numbers.

    deci,size : natural32 := 0;

  begin
    new_line;
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    return size;
  end Ask_for_Size;

  procedure Main is

  -- DESCRPITION :
  --  Prompts the user for the precision and the dimension.

    ans,sys : character;
    n,sz : natural32 := 0;

  begin
    new_line;
    put_line("MENU to test encapsulation of gradient evaluations ...");
    put_line("  0. standard double precision; or");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision; or");
    put_line("  3. arbitrary multiprecision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    new_line;
    put("Test system or one polynomial ? (s/p) ");
    Ask_Alternative(sys,"sp");
    if sys = 'p' then
      new_line;
      put("Give the number of variables : "); get(n);
      Symbol_Table.Init(n);
      case ans is 
        when '0' => Standard_Test(n);
        when '1' => DoblDobl_Test(n);
        when '2' => QuadDobl_Test(n);
        when '3' => sz := Ask_for_Size;
                    Multprec_Test(n,sz);
        when others => null;
      end case;
    else
      case ans is
        when '0' => Standard_System_Test;
        when '1' => DoblDobl_System_Test;
        when '2' => QuadDobl_System_Test;
        when '3' => sz := Ask_for_Size;
                    Multprec_System_Test(sz);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_gradcirc;
