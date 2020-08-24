with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Random_Numbers;            use QuadDobl_Random_Numbers;
with QuadDobl_Mathematical_Functions;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;

procedure ts_tridbl is

-- DESCRIPTION :
--   This procedure develops triple double arithmetic,
--   interpolating between double double and quad double arithmetic,
--   with generated code from the CAMPARY library.

  function to_triple_double ( x : quad_double ) return triple_double is

  -- DESCRIPTION :
  --   Copies the first three words of x into the result.

    res : constant triple_double
        := create(hihi_part(x),lohi_part(x),hilo_part(x));

  begin
    return res;
  end to_triple_double;

  function create ( x : triple_double ) return quad_double is

  -- DESCRIPTION :
  --   Copies the words of x into the first three words of the result.

    res : constant quad_double
        := create(hi_part(x),mi_part(x),lo_part(x),0.0);

  begin
    return res;
  end create;

  procedure Test_Basic_Arithmetic is

  -- DESCRIPTION :
  --   Generates random numbers and compares with quad double arithmetic.
  --   The sin and cos functions are applied to have four nonzero words.

    qdx : constant quad_double
        := QuadDobl_Mathematical_Functions.sin(Random);
    qdy : constant quad_double
        := QuadDobl_Mathematical_Functions.cos(Random);
    qdz : constant quad_double := qdx + qdy;
    tdx : constant triple_double := to_triple_double(qdx);
    tdy : constant triple_double := to_triple_double(qdy);
    tdz : constant triple_double := tdx + tdy;
    qtx : constant quad_double := create(tdx);
    qty : constant quad_double := create(tdy);
    qtz : constant quad_double := create(tdz);
    df1 : constant quad_double := qdz - qtz;
    tdv : constant triple_double := tdz - tdy;
    qtv : constant quad_double := create(tdv);
    df2 : constant quad_double := qdx - qtv;
    qdp : constant quad_double := qdx * qdy;
    tdp : constant triple_double := tdx * tdy;
    qtp : constant quad_double := create(tdp);
    df3 : constant quad_double := qdp - qtp;
    dy : constant double_float := hihi_part(qdy);
    qdxdy : constant quad_double := qdx*dy;
    tdxdy : constant triple_double := tdx*dy;
    qtxdy : constant quad_double := create(tdxdy);
    df4 : constant quad_double := qdxdy - qtxdy;
    qdq : constant quad_double := qdx / qdy;
    tdq : constant triple_double := tdx / tdy;
    qtq : constant quad_double := create(tdq);
    df5 : constant quad_double := qdq - qtq;

  begin
    put_line("The sum of two random quad doubles :");
    put("    x : "); put(qdx); new_line;
    put("    y : "); put(qdy); new_line;
    put("x + y : "); put(qdz); new_line;
    put_line("The same sum with triple doubles :");
    put("    x : "); put(qtx); new_line;
    put("    y : "); put(qty); new_line;
    put("x + y : "); put(qtz); new_line;
    put_line("The difference :"); put(df1); new_line;
    put("(x+y)-y : "); put(qtv); new_line;
    put("      x : "); put(qtx); new_line;
    put_line("The difference :"); put(df2); new_line;
    put("qd x*y : "); put(qdp); new_line;
    put("td x*y : "); put(qtp); new_line;
    put_line("The difference :"); put(df3); new_line;
    put("qd x*dy : "); put(qdxdy); new_line;
    put("td x*dy : "); put(qtxdy); new_line;
    put_line("The difference :"); put(df4); new_line;
    put("qd x/y : "); put(qdq); new_line;
    put("td x/y : "); put(qtq); new_line;
    put_line("The difference :"); put(df5); new_line;
  end Test_Basic_Arithmetic;

  procedure Test_Read is

  -- DESCRIPTION :
  --   Reads a 50-digit approximation A for sqrt(2) from a string
  --   and shows the result of A*A - 2.

  --   >>> from sympy import evalf, sqrt
  --   >>> s2 = sqrt(2).evalf(50)
  --   >>> s2
  --   1.4142135623730950488016887242096980785696718753769
  --   >>> s2*s2 - 2
  --   -2.6727647100921956461405364671514818788151968801050e-51
    
    sqrt2 : constant string
          := "1.4142135623730950488016887242096980785696718753769";
    x,r : triple_double;
    two : constant triple_double := create(2.0);
    fail : boolean;
    qdx,qdr : quad_double;

  begin
    read(sqrt2,x,fail);
    new_line;
    if fail then
      put_line("The read procedure reports failure!");
    else
      put_line("The read procedure reports no failure.");
      qdx := create(x);
      put("qd x : "); put(qdx); new_line;
      r := x*x - 2.0;
      qdr := create(r);
      put_line("x*x - 2.0 : "); put(qdr); new_line;
      r := x*x - two;
      qdr := create(r);
      put_line("x*x - two : "); put(qdr); new_line;
    end if;
  end Test_Read;

  procedure Test_io is

  -- DESCRIPTION :
  --   Prompts the user for a number, reads and writes a triple double.

    x,y : triple_double;
    ans : character;
 
  begin
    new_line;
    loop
      put("Give x : "); get(x);
      put(" --> x : "); put(x); new_line;
      put("Give pair x y : "); get(x,y); 
      put(" --> x : "); put(x); new_line;
      put(" --> y : "); put(y); new_line;
      put("More tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_io;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to select a test
  --   and then runs the test.

    ans : character;

  begin
    new_line;
    put_line("Testing triple double arithmetic ...");
    put_line("  1. basic arithmetic");
    put_line("  2. test reading from string");
    put_line("  3. input and output");
    put("Type 1, 2, or 3 to select a test : "); Ask_Alternative(ans,"123");
    case ans is
      when '1' => Test_Basic_Arithmetic;
      when '2' => Test_Read;
      when '3' => Test_io;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_tridbl;
