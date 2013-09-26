with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Symbol_Table;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Multprec_Floating_Vectors;
with Multprec_Floating_Matrices;
with Multprec_Floating_Matrices_io;      use Multprec_Floating_Matrices_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;
with Multprec_Complex_Matrices_io;       use Multprec_Complex_Matrices_io;
with Standard_Floating_Polynomials;
with Standard_Floating_Poly_Functions;
with Multprec_Floating_Polynomials;
with Multprec_Floating_Poly_Functions;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Functions;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Standard_Complex_to_Real_Poly;      use Standard_Complex_to_Real_Poly;
with Multprec_Complex_to_Real_Poly;      use Multprec_Complex_to_Real_Poly;
with Standard_Complex_Poly_Systems;
with Standard_Floating_Poly_Systems;
with Standard_Floating_Poly_Systems_io;  use Standard_Floating_Poly_Systems_io;
with Standard_Floating_Jaco_Matrices;
with Standard_Complex_Jaco_Matrices;
with Multprec_Complex_Poly_Systems;
with Multprec_Floating_Poly_Systems;
with Multprec_Floating_Poly_Systems_io;  use Multprec_Floating_Poly_Systems_io;
with Multprec_Floating_Jaco_Matrices;
with Multprec_Complex_Jaco_Matrices;

procedure ts_rpol is

-- DESCRIPTION :
--   This routine provides basic testing routines for real polynomials.

-- NOTICE the bug: type in "x*y + (2-I)*y;" and look at the output...

  procedure Test_Standard_io is

  -- DESCRIPTION :
  --   Tests the input/output of a polynomial in several variables
  --   with standard complex coefficients and the conversions.

    n : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    q : Standard_Floating_Polynomials.Poly;

  begin
    new_line;
    put_line
      ("Interactive testing of input/output of standard complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    new_line;
    put_line("Your polynomial : "); put(p); new_line;
    put_line("Dropping the imaginary parts of the coefficients ...");
    q := Convert_Complex_to_Real(p);
    Standard_Complex_Polynomials.Clear(p);
    put_line("Converting back to a complex poly for writing ...");
    p := Convert_Real_to_Complex(q);
    put_line("Your polynomial : "); put(p); new_line;
    Symbol_Table.Clear;
  end Test_Standard_io;

  procedure Test_Multprec_io is

  -- DESCRIPTION :
  --   Tests the input/output of a polynomial in several variables
  --   with multiprecision complex coefficients and the conversions.

    n : natural32 := 0;
    p : Multprec_Complex_Polynomials.Poly;
    q : Multprec_Floating_Polynomials.Poly;

  begin
    new_line;
    put_line
      ("Interactive testing of input/output of multprec complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    new_line;
    put_line("Your polynomial : "); put(p); new_line;
    put_line("Dropping the imaginary parts of the coefficients ...");
    q := Convert_Complex_to_Real(p);
    Multprec_Complex_Polynomials.Clear(p);
    put_line("Converting back to a complex poly for writing ...");
    p := Convert_Real_to_Complex(q);
    put_line("Your polynomial : "); put(p); new_line;
    Symbol_Table.Clear;
  end Test_Multprec_io;

  procedure Test_Standard_Eval is

  -- DESCRIPTION :
  --   Tests the evaluation of a polynomial in several variables
  --   with real coefficients.

    n : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    q : Standard_Floating_Polynomials.Poly;

    use Standard_Complex_Numbers;

  begin
    new_line;
    put_line
      ("Interactive testing of evaluation of standard complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    new_line;
    put_line("Your polynomial : "); put(p); new_line;
    put_line("Dropping the imaginary parts of the coefficients ...");
    q := Convert_Complex_to_Real(p);
    declare
      eq : constant Standard_Floating_Poly_Functions.Eval_Poly
         := Standard_Floating_Poly_Functions.Create(q);
      rx : Standard_Floating_Vectors.Vector(1..integer32(n));
      cx : Standard_Complex_Vectors.Vector(1..integer32(n));
      y1,y2 : double_float;
      z : Complex_Number;
    begin
      put("Give "); put(n,1); put_line(" floating numbers for x :");
      for i in rx'range loop
        get(rx(i));
        cx(i) := Create(rx(i));
      end loop;
      y1 := Standard_Floating_Poly_Functions.Eval(q,rx);
      put("The value at x : "); put(y1); new_line;
      y2 := Standard_Floating_Poly_Functions.Eval(eq,rx);
      put("   with Horner : "); put(y2); new_line;
      z := Standard_Complex_Poly_Functions.Eval(p,cx);
      put("  complex eval : "); put(z); new_line;
    end;
  end Test_Standard_Eval;

  procedure Test_Multprec_Eval is

  -- DESCRIPTION :
  --   Tests the evaluation of a polynomial in several variables
  --   with multprecision real coefficients.

    n : natural32 := 0;
    p : Multprec_Complex_Polynomials.Poly;
    q : Multprec_Floating_Polynomials.Poly;

    use Multprec_Complex_Numbers;

  begin
    new_line;
    put_line
      ("Interactive testing of evaluation of multprec complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    new_line;
    put_line("Your polynomial : "); put(p); new_line;
    put_line("Dropping the imaginary parts of the coefficients ...");
    q := Convert_Complex_to_Real(p);
    declare
      eq : Multprec_Floating_Poly_Functions.Eval_Poly
         := Multprec_Floating_Poly_Functions.Create(q);
      rx : Multprec_Floating_Vectors.Vector(1..integer32(n));
      cx : Multprec_Complex_Vectors.Vector(1..integer32(n));
      y1,y2 : Floating_Number;
      z : Complex_Number;
    begin
      put("Give "); put(n,1); put_line(" floating numbers for x :");
      for i in rx'range loop
        get(rx(i));
        cx(i) := Create(rx(i));
      end loop;
      y1 := Multprec_Floating_Poly_Functions.Eval(q,rx);
      put("The value at x : "); put(y1); new_line;
      y2 := Multprec_Floating_Poly_Functions.Eval(eq,rx);
      put("   with Horner : "); put(y2); new_line;
      z := Multprec_Complex_Poly_Functions.Eval(p,cx);
      put("  complex eval : "); put(z); new_line;
      Multprec_Floating_Poly_Functions.Clear(eq);
    end;
  end Test_Multprec_Eval;

  procedure Test_Standard_Jaco_Eval is

    clp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlp : Standard_Floating_Poly_Systems.Link_to_Poly_Sys;
    n : integer32;

  begin
    new_line;
    put_line("Creation and Evaluation of standard Jacobi matrices.");
    new_line;
    get(rlp);
    clp := new Standard_Complex_Poly_Systems.Poly_Sys'
                 (Convert_Real_to_Complex(rlp.all));
    n := integer32(Standard_Floating_Polynomials.Number_of_Unknowns
                    (rlp(rlp'first)));
    declare
      rjm : constant Standard_Floating_Jaco_Matrices.Jaco_Mat(rlp'range,1..n)
          := Standard_Floating_Jaco_Matrices.Create(rlp.all);
      erjm : Standard_Floating_Jaco_Matrices.Eval_Jaco_Mat(rlp'range,1..n)
           := Standard_Floating_Jaco_Matrices.Create(rjm);
      cjm : Standard_Complex_Jaco_Matrices.Jaco_Mat(clp'range,1..n)
          := Standard_Complex_Jaco_Matrices.Create(clp.all);
      rx : Standard_Floating_Vectors.Vector(1..n);
      cx : Standard_Complex_Vectors.Vector(1..n);
      ry1,ry2 : Standard_Floating_Matrices.Matrix(rjm'range,1..n);
      cy : Standard_Complex_Matrices.Matrix(rjm'range,1..n);
    begin
      put("Give "); put(rx'last,1); put_line(" floats for x :");
      for i in rx'range loop
        get(rx(i));
        cx(i) := Standard_Complex_Numbers.Create(rx(i));
      end loop;
      ry1 := Standard_Floating_Jaco_Matrices.Eval(rjm,rx);
      ry2 := Standard_Floating_Jaco_Matrices.Eval(erjm,rx);
      cy := Standard_Complex_Jaco_Matrices.Eval(cjm,cx);
      put_line("The real Jacobi matrix at x :"); put(ry1);
      put_line(" evaluated with Horner at x :"); put(ry2);
      put_line("The complex Jacobi matrix at x :"); put(cy);
    end;
  end Test_Standard_Jaco_Eval;

  procedure Test_Multprec_Jaco_Eval is

    clp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlp : Multprec_Floating_Poly_Systems.Link_to_Poly_Sys;
    n : integer32 := 0;

  begin
    new_line;
    put_line("Creation and Evaluation of multiprecision Jacobi matrices.");
    new_line;
    get(rlp);
    clp := new Multprec_Complex_Poly_Systems.Poly_Sys'
                 (Convert_Real_to_Complex(rlp.all));
    n := integer32(Multprec_Floating_Polynomials.Number_of_Unknowns
                     (rlp(rlp'first)));
    declare
      rjm : Multprec_Floating_Jaco_Matrices.Jaco_Mat(rlp'range,1..n)
          := Multprec_Floating_Jaco_Matrices.Create(rlp.all);
      erjm : Multprec_Floating_Jaco_Matrices.Eval_Jaco_Mat(rlp'range,1..n)
           := Multprec_Floating_Jaco_Matrices.Create(rjm);
      cjm : Multprec_Complex_Jaco_Matrices.Jaco_Mat(clp'range,1..n)
          := Multprec_Complex_Jaco_Matrices.Create(clp.all);
      rx : Multprec_Floating_Vectors.Vector(1..n);
      cx : Multprec_Complex_Vectors.Vector(1..n);
      ry1,ry2 : Multprec_Floating_Matrices.Matrix(rjm'range,1..n);
      cy : Multprec_Complex_Matrices.Matrix(rjm'range,1..n);
    begin
      put("Give "); put(rx'last,1); put_line(" floats for x :");
      for i in rx'range loop
        get(rx(i));
        cx(i) := Multprec_Complex_Numbers.Create(rx(i));
      end loop;
      ry1 := Multprec_Floating_Jaco_Matrices.Eval(rjm,rx);
      ry2 := Multprec_Floating_Jaco_Matrices.Eval(erjm,rx);
      cy := Multprec_Complex_Jaco_Matrices.Eval(cjm,cx);
      put_line("The real Jacobi matrix at x :"); put(ry1);
      put_line(" evaluated with Horner at x :"); put(ry2);
      put_line("The complex Jacobi matrix at x :"); put(cy);
    end;
  end Test_Multprec_Jaco_Eval;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of the operations on real polynomials.");
    new_line;
    put_line("MENU with testing operations :");
    put_line("  1. i/o of polynomials with standard coefficients;");
    put_line("  2. i/o of polynomials with multiprecision coefficients;");
    put_line("  3. polynomial evaluation with standard numbers;");
    put_line("  4. polynomial evaluation with multiprecision numbers;");
    put_line("  5. evaluation of Jacobi matrices in standard precision;");
    put_line("  6. evaluation of Jacobi matrices in multiprecision.");
    put("Type 1, 2, 3, 4, 5, or 6 to choose a test : ");
    Ask_Alternative(ans,"123456");
    case ans is
      when '1' => Test_Standard_io;
      when '2' => Test_Multprec_io;
      when '3' => Test_Standard_Eval;
      when '4' => Test_Multprec_Eval;
      when '5' => Test_Standard_Jaco_Eval;
      when '6' => Test_Multprec_Jaco_Eval;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_rpol;
