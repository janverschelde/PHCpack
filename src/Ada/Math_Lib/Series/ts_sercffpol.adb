with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Series;
with Standard_Complex_Series_io;         use Standard_Complex_Series_io;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_Vectors_io; use Standard_Complex_Series_Vectors_io;
with Standard_Complex_Series_VecVecs;
with Standard_Complex_Series_Matrices;
with Standard_Random_Series_Vectors;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_io;         use DoblDobl_Complex_Series_io;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors_io; use DoblDobl_Complex_Series_Vectors_io;
with DoblDobl_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_Matrices;
with DoblDobl_Random_Series_Vectors;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_io;         use QuadDobl_Complex_Series_io;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors_io; use QuadDobl_Complex_Series_Vectors_io;
with QuadDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_Matrices;
with QuadDobl_Random_Series_Vectors;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Polynomials_io;    use Standard_CSeries_Polynomials_io;
with Standard_CSeries_Poly_Functions;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_Systems_io;   use Standard_CSeries_Poly_Systems_io;
with Standard_CSeries_Poly_SysFun;
with Standard_CSeries_Jaco_Matrices;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Polynomials_io;    use DoblDobl_CSeries_Polynomials_io;
with DoblDobl_CSeries_Poly_Functions;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems_io;   use DoblDobl_CSeries_Poly_Systems_io;
with DoblDobl_CSeries_Poly_SysFun;
with DoblDobl_CSeries_Jaco_Matrices;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Polynomials_io;    use QuadDobl_CSeries_Polynomials_io;
with QuadDobl_CSeries_Poly_Functions;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems_io;   use QuadDobl_CSeries_Poly_Systems_io;
with QuadDobl_CSeries_Poly_SysFun;
with QuadDobl_CSeries_Jaco_Matrices;
with Random_Series_Polynomials;          use Random_Series_Polynomials;

procedure ts_sercffpol is

-- DESCRIPTION :
--   Test on coefficient parameter evaluation for polynomials with
--   series coefficients.

  procedure Test_Eval
              ( n,d : in integer32;
                p : in out Standard_CSeries_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Generates a random series vector to test the evaluation.
 
  -- ON ENTRY :
  --   n        number of variables in p;
  --   d        degree of the series coefficients;
  --   p        a polynomial in n variables with series coefficients.

    x : constant Standard_Complex_Series_Vectors.Vector
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    y : constant Standard_Complex_Series.Link_to_Series
      := Standard_CSeries_Poly_Functions.Eval(p,x);
    m : constant integer32
      := integer32(Standard_CSeries_Polynomials.Number_of_Terms(p));
    c : constant Standard_Complex_Series_Vectors.Vector(1..m)
      := Standard_CSeries_Poly_Functions.Coeff(p);
    f : constant Standard_CSeries_Poly_Functions.Eval_Coeff_Poly
      := Standard_CSeries_Poly_Functions.Create(p);
    z : constant Standard_Complex_Series.Link_to_Series
      := Standard_CSeries_Poly_Functions.Eval(f,c,x);

  begin
    put_line("Evaluation at a random vector ..."); put(x);
    put_line("The value at the polynomial :"); put(y);
    put_line("The value at the coefficient polynomial :"); put(z);
  end Test_Eval;

  procedure Test_Eval
              ( n,d : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Generates a random series vector to test the evaluation.
 
  -- ON ENTRY :
  --   n        number of variables in p;
  --   d        degree of the series coefficients;
  --   p        a polynomial in n variables with series coefficients.

    x : constant Standard_Complex_Series_Vectors.Vector
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    y : constant Standard_Complex_Series_Vectors.Vector
      := Standard_CSeries_Poly_SysFun.Eval(p,x);
    c : constant Standard_Complex_Series_VecVecs.VecVec(p'range)
      := Standard_CSeries_Poly_SysFun.Coeff(p);
    f : constant Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(p'range)
      := Standard_CSeries_Poly_SysFun.Create(p);
    z : constant Standard_Complex_Series_Vectors.Vector
      := Standard_CSeries_Poly_SysFun.Eval(f,c,x);
    jpm : constant Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,1..n)
        := Standard_CSeries_Jaco_Matrices.Create(p);
    ejm : Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(p'range,1..n);
    mlt : Standard_CSeries_Jaco_Matrices.Mult_Factors(p'range,1..n);
    A,B : Standard_Complex_Series_Matrices.Matrix(p'range,1..n);

  begin
    put_line("Evaluation at a random vector ..."); put(x);
    put_line("The value at the polynomial :"); put(y);
    put_line("The value at the coefficient polynomial :"); put(z);
    Standard_CSeries_Jaco_Matrices.Create(p,ejm,mlt);
    A := Standard_CSeries_Jaco_Matrices.Eval(jpm,x);
    B := Standard_CSeries_Jaco_Matrices.Eval(ejm,mlt,c,x);
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A("); put(i,1); put(","); put(j,1); put_line("):");
        put(A(i,j)); new_line;
        put("B("); put(i,1); put(","); put(j,1); put_line("):");
        put(B(i,j)); new_line;
      end loop;
    end loop;
  end Test_Eval;

  procedure Test_Eval
              ( n,d : in integer32;
                p : in DoblDobl_CSeries_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Generates a random series vector to test the evaluation.
 
  -- ON ENTRY :
  --   n        number of variables in p;
  --   d        degree of the series coefficients;
  --   p        a polynomial in n variables with series coefficients.

    x : constant DoblDobl_Complex_Series_Vectors.Vector
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    y : constant DoblDobl_Complex_Series.Link_to_Series
      := DoblDobl_CSeries_Poly_Functions.Eval(p,x);
    m : constant integer32
      := integer32(DoblDobl_CSeries_Polynomials.Number_of_Terms(p));
    c : constant DoblDobl_Complex_Series_Vectors.Vector(1..m)
      := DoblDobl_CSeries_Poly_Functions.Coeff(p);
    f : constant DoblDobl_CSeries_Poly_Functions.Eval_Coeff_Poly
      := DoblDobl_CSeries_Poly_Functions.Create(p);
    z : constant DoblDobl_Complex_Series.Link_to_Series
      := DoblDobl_CSeries_Poly_Functions.Eval(f,c,x);

  begin
    put_line("Evaluation at a random vector ..."); put(x);
    put_line("The value at the polynomial :"); put(y);
    put_line("The value at the coefficient polynomial :"); put(z);
  end Test_Eval;

  procedure Test_Eval
              ( n,d : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Generates a random series vector to test the evaluation.
 
  -- ON ENTRY :
  --   n        number of variables in p;
  --   d        degree of the series coefficients;
  --   p        a polynomial in n variables with series coefficients.

    x : constant DoblDobl_Complex_Series_Vectors.Vector
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    y : constant DoblDobl_Complex_Series_Vectors.Vector
      := DoblDobl_CSeries_Poly_SysFun.Eval(p,x);
    c : constant DoblDobl_Complex_Series_VecVecs.VecVec(p'range)
      := DoblDobl_CSeries_Poly_SysFun.Coeff(p);
    f : constant DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(p'range)
      := DoblDobl_CSeries_Poly_SysFun.Create(p);
    z : constant DoblDobl_Complex_Series_Vectors.Vector
      := DoblDobl_CSeries_Poly_SysFun.Eval(f,c,x);
    jpm : constant DoblDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,1..n)
        := DoblDobl_CSeries_Jaco_Matrices.Create(p);
    ejm : DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(p'range,1..n);
    mlt : DoblDobl_CSeries_Jaco_Matrices.Mult_Factors(p'range,1..n);
    A,B : DoblDobl_Complex_Series_Matrices.Matrix(p'range,1..n);

  begin
    put_line("Evaluation at a random vector ..."); put(x);
    put_line("The value at the polynomial :"); put(y);
    put_line("The value at the coefficient polynomial :"); put(z);
    DoblDobl_CSeries_Jaco_Matrices.Create(p,ejm,mlt);
    A := DoblDobl_CSeries_Jaco_Matrices.Eval(jpm,x);
    B := DoblDobl_CSeries_Jaco_Matrices.Eval(ejm,mlt,c,x);
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A("); put(i,1); put(","); put(j,1); put_line("):");
        put(A(i,j)); new_line;
        put("B("); put(i,1); put(","); put(j,1); put_line("):");
        put(B(i,j)); new_line;
      end loop;
    end loop;
  end Test_Eval;

  procedure Test_Eval
              ( n,d : in integer32;
                p : in QuadDobl_CSeries_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Generates a random series vector to test the evaluation.
 
  -- ON ENTRY :
  --   n        number of variables in p;
  --   d        degree of the series coefficients;
  --   p        a polynomial in n variables with series coefficients.

    x : constant QuadDobl_Complex_Series_Vectors.Vector
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    y : constant QuadDobl_Complex_Series.Link_to_Series
      := QuadDobl_CSeries_Poly_Functions.Eval(p,x);
    m : constant integer32
      := integer32(QuadDobl_CSeries_Polynomials.Number_of_Terms(p));
    c : constant QuadDobl_Complex_Series_Vectors.Vector(1..m)
      := QuadDobl_CSeries_Poly_Functions.Coeff(p);
    f : constant QuadDobl_CSeries_Poly_Functions.Eval_Coeff_Poly
      := QuadDobl_CSeries_Poly_Functions.Create(p);
    z : constant QuadDobl_Complex_Series.Link_to_Series
      := QuadDobl_CSeries_Poly_Functions.Eval(f,c,x);

  begin
    put_line("Evaluation at a random vector ..."); put(x);
    put_line("The value at the polynomial :"); put(y);
    put_line("The value at the coefficient polynomial :"); put(z);
  end Test_Eval;

  procedure Test_Eval
              ( n,d : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Generates a random series vector to test the evaluation.
 
  -- ON ENTRY :
  --   n        number of variables in p;
  --   d        degree of the series coefficients;
  --   p        a polynomial in n variables with series coefficients.

    x : constant QuadDobl_Complex_Series_Vectors.Vector
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    y : constant QuadDobl_Complex_Series_Vectors.Vector
      := QuadDobl_CSeries_Poly_SysFun.Eval(p,x);
    c : constant QuadDobl_Complex_Series_VecVecs.VecVec(p'range)
      := QuadDobl_CSeries_Poly_SysFun.Coeff(p);
    f : constant QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(p'range)
      := QuadDobl_CSeries_Poly_SysFun.Create(p);
    z : constant QuadDobl_Complex_Series_Vectors.Vector
      := QuadDobl_CSeries_Poly_SysFun.Eval(f,c,x);
    jpm : constant QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,1..n)
        := QuadDobl_CSeries_Jaco_Matrices.Create(p);
    ejm : QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(p'range,1..n);
    mlt : QuadDobl_CSeries_Jaco_Matrices.Mult_Factors(p'range,1..n);
    A,B : QuadDobl_Complex_Series_Matrices.Matrix(p'range,1..n);

  begin
    put_line("Evaluation at a random vector ..."); put(x);
    put_line("The value at the polynomial :"); put(y);
    put_line("The value at the coefficient polynomial :"); put(z);
    QuadDobl_CSeries_Jaco_Matrices.Create(p,ejm,mlt);
    A := QuadDobl_CSeries_Jaco_Matrices.Eval(jpm,x);
    B := QuadDobl_CSeries_Jaco_Matrices.Eval(ejm,mlt,c,x);
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A("); put(i,1); put(","); put(j,1); put_line("):");
        put(A(i,j)); new_line;
        put("B("); put(i,1); put(","); put(j,1); put_line("):");
        put(B(i,j)); new_line;
      end loop;
    end loop;
  end Test_Eval;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for the number of variables, number of terms,
  --   degree of a polynomial and the degree of the series coefficients,
  --   to generate a random polynomial to test in double precision.

    nq,nv,m,dp,ds : natural32 := 0;
    p : Standard_CSeries_Polynomials.Poly;

  begin
    new_line;
    put_line("Test on coefficient polynomials with series coefficients.");
    put("-> give the number of polynomials : "); get(nq);
    put("-> give the number of variables : "); get(nv);
    put("-> give the number of terms : "); get(m);
    put("-> give the degree of the polynomial(s) : "); get(dp);
    put("-> give the degree of the series coefficients : "); get(ds);
    if nq = 1 then
      p := Standard_Random_Polynomial(nv,m,dp,ds);
      put_line("A random polynomial with series coefficients :"); put(p);
      Test_Eval(integer32(nv),integer32(ds),p);
    else
      declare
        q : constant Standard_CSeries_Poly_Systems.Poly_Sys
          := Standard_Random_System(nq,nv,m,dp,ds);
      begin
        put_line("A random system with series coefficients :"); put(q);
        Test_Eval(integer32(nv),integer32(ds),q);
      end;
    end if;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for the number of variables, number of terms,
  --   degree of a polynomial and the degree of the series coefficients,
  --   to generate a random polynomial to test in double double precision.

    nq,nv,m,dp,ds : natural32 := 0;
    p : DoblDobl_CSeries_Polynomials.Poly;

  begin
    new_line;
    put_line("Test on coefficient polynomials with series coefficients.");
    put("-> give the number of polynomials : "); get(nq);
    put("-> give the number of variables : "); get(nv);
    put("-> give the number of terms : "); get(m);
    put("-> give the degree of the polynomial(s) : "); get(dp);
    put("-> give the degree of the series coefficients : "); get(ds);
    if nq = 1 then
      p := DoblDobl_Random_Polynomial(nv,m,dp,ds);
      put_line("A random polynomial with series coefficients :"); put(p);
      Test_Eval(integer32(nv),integer32(ds),p);
    else
      declare
        q : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys
          := DoblDobl_Random_System(nq,nv,m,dp,ds);
      begin
        put_line("A random system with series coefficients :"); put(q);
        Test_Eval(integer32(nv),integer32(ds),q);
      end;
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for the number of variables, number of terms,
  --   degree of a polynomial and the degree of the series coefficients,
  --   to generate a random polynomial to test in quad double precision.

    nq,nv,m,dp,ds : natural32 := 0;
    p : QuadDobl_CSeries_Polynomials.Poly;

  begin
    new_line;
    put_line("Test on coefficient polynomials with series coefficients.");
    put("-> give the number of polynomials : "); get(nq);
    put("-> give the number of variables : "); get(nv);
    put("-> give the number of terms : "); get(m);
    put("-> give the degree of the polynomial(s) : "); get(dp);
    put("-> give the degree of the series coefficients : "); get(ds);
    if nq = 1 then
      p := QuadDobl_Random_Polynomial(nv,m,dp,ds);
      put_line("A random polynomial with series coefficients :"); put(p);
      Test_Eval(integer32(nv),integer32(ds),p);
    else
      declare
        q : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys
          := QuadDobl_Random_System(nq,nv,m,dp,ds);
      begin
        put_line("A random system with series coefficients :"); put(q);
        Test_Eval(integer32(nv),integer32(ds),q);
      end;
    end if;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision and then launches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sercffpol;
