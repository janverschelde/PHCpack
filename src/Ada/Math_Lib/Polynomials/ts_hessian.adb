with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Symbol_Table;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Hessians;
with DoblDobl_Complex_Hessians;
with QuadDobl_Complex_Hessians;

procedure ts_hessian is

  procedure Write ( h : in Standard_Complex_Hessians.Link_to_Hessian ) is

  -- DESCRIPTION :
  --   Writes the Hessian to screen.
 
  begin
    for i in h'range(1) loop
      for j in h'range(2) loop
        put("H["); put(i,1); put(","); put(j,1); put("] : ");
        put(h(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Write ( h : in DoblDobl_Complex_Hessians.Link_to_Hessian ) is

  -- DESCRIPTION :
  --   Writes the Hessian to screen.
 
  begin
    for i in h'range(1) loop
      for j in h'range(2) loop
        put("H["); put(i,1); put(","); put(j,1); put("] : ");
        put(h(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Write ( h : in QuadDobl_Complex_Hessians.Link_to_Hessian ) is

  -- DESCRIPTION :
  --   Writes the Hessian to screen.
 
  begin
    for i in h'range(1) loop
      for j in h'range(2) loop
        put("H["); put(i,1); put(","); put(j,1); put("] : ");
        put(h(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Write ( m : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the matrix m to screen.
 
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        put("m["); put(i,1); put(","); put(j,1); put("] : ");
        put(m(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Write ( m : in DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the matrix m to screen.
 
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        put("m["); put(i,1); put(","); put(j,1); put("] : ");
        put(m(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Write ( m : in QuadDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the matrix m to screen.
 
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        put("m["); put(i,1); put(","); put(j,1); put("] : ");
        put(m(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Standard_Test 
              ( p : in Standard_Complex_Polynomials.Poly;
                n : in natural32 ) is

  -- DESCRIPTION :
  --   Tests the Hessian of p, a polynomial in n variables,
  --   in standard double precision.

    h : Standard_Complex_Hessians.Link_to_Hessian
      := Standard_Complex_Hessians.Create(p);
    x : Standard_Complex_Vectors.Vector(1..integer32(n));
    m : Standard_Complex_Matrices.Matrix(x'range,x'range);

  begin
    Write(h);
    new_line;
    put("Reading "); put(n,1);
    put_line(" complex numbers for evaluation : ");
    for i in x'range loop
      put(i,1); put(" : "); get(x(i));
    end loop;
    put_line("Evaluation at "); put_line(x);
    m := Standard_Complex_Hessians.Eval(h,x);
    Write(m);
  end Standard_Test;

  procedure Standard_Test 
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests the Hessian of p, a polynomial system,
  --   in standard double precision.

    h : Standard_Complex_Hessians.Array_of_Hessians(p'range)
      := Standard_Complex_Hessians.Create(p);

  begin
    for i in h'range loop
      put("The Hessian for polynomial ");
      put(i,1); put_line(" :");
      Write(h(i));
    end loop;
  end Standard_Test;

  procedure DoblDobl_Test 
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                n : in natural32 ) is

  -- DESCRIPTION :
  --   Tests the Hessian of p, a polynomial in n variables,
  --   in double double precision.

    h : DoblDobl_Complex_Hessians.Link_to_Hessian
      := DoblDobl_Complex_Hessians.Create(p);
    x : DoblDobl_Complex_Vectors.Vector(1..integer32(n));
    m : DoblDobl_Complex_Matrices.Matrix(x'range,x'range);

  begin
    Write(h);
    new_line;
    put("Reading "); put(n,1);
    put_line(" complex numbers for evaluation : ");
    for i in x'range loop
      put(i,1); put(" : "); get(x(i));
    end loop;
    put_line("Evaluation at "); put_line(x);
    m := DoblDobl_Complex_Hessians.Eval(h,x);
    Write(m);
  end DoblDobl_Test;

  procedure DoblDobl_Test 
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests the Hessian of p, a polynomial system,
  --   in double double precision.

    h : DoblDobl_Complex_Hessians.Array_of_Hessians(p'range)
      := DoblDobl_Complex_Hessians.Create(p);

  begin
    for i in h'range loop
      put("The Hessian for polynomial ");
      put(i,1); put_line(" :");
      Write(h(i));
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test 
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                n : in natural32 ) is

  -- DESCRIPTION :
  --   Tests the Hessian of p, a polynomial in n variables,
  --   in quad double precision.

    h : QuadDobl_Complex_Hessians.Link_to_Hessian
      := QuadDobl_Complex_Hessians.Create(p);
    x : QuadDobl_Complex_Vectors.Vector(1..integer32(n));
    m : QuadDobl_Complex_Matrices.Matrix(x'range,x'range);

  begin
    Write(h);
    new_line;
    put("Reading "); put(n,1);
    put_line(" complex numbers for evaluation : ");
    for i in x'range loop
      put(i,1); put(" : "); get(x(i));
    end loop;
    put_line("Evaluation at "); put_line(x);
    m := QuadDobl_Complex_Hessians.Eval(h,x);
    Write(m);
  end QuadDobl_Test;

  procedure QuadDobl_Test 
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests the Hessian of p, a polynomial system,
  --   in quad double precision.

    h : QuadDobl_Complex_Hessians.Array_of_Hessians(p'range)
      := QuadDobl_Complex_Hessians.Create(p);

  begin
    for i in h'range loop
      put("The Hessian for polynomial ");
      put(i,1); put_line(" :");
      Write(h(i));
    end loop;
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of variables and a polynomial.
  --   Then runs some tests on the Hessian in standard double precision.

    ans : character;
    n : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    s : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    put("Polynomial system ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading a polynomial system ...");
      get(s);
      Standard_Test(s.all);
    else
      put("Give the number of variables : "); get(n);
      Symbol_Table.Init(n);
      put("Give a polynomial : "); get(p);
      new_line;
      put("-> your polynomial : "); put(p); new_line;
      Standard_Test(p,n);
    end if;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of variables and a polynomial.
  --   Then runs some tests on the Hessian in double double precision.

    ans : character;
    n : natural32 := 0;
    p : DoblDobl_Complex_Polynomials.Poly;
    s : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    put("Polynomial system ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading a polynomial system ...");
      get(s);
      DoblDobl_Test(s.all);
    else
      put("Give the number of variables : "); get(n);
      Symbol_Table.Init(n);
      put("Give a polynomial : "); get(p);
      new_line;
      put("-> your polynomial : "); put(p); new_line;
      DoblDobl_Test(p,n);
    end if;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of variables and a polynomial.
  --   Then runs some tests on the Hessian in quad double precision.

    ans : character;
    n : natural32 := 0;
    p : QuadDobl_Complex_Polynomials.Poly;
    s : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    put("Polynomial system ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading a polynomial system ...");
      get(s);
      QuadDobl_Test(s.all);
    else
      put("Give the number of variables : "); get(n);
      Symbol_Table.Init(n);
      put("Give a polynomial : "); get(p);
      new_line;
      put("-> your polynomial : "); put(p); new_line;
      QuadDobl_Test(p,n);
    end if;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision and then launches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision : ");
    put_line("  0. double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_hessian;
