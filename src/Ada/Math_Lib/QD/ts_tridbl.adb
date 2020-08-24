with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Random_Numbers;            use QuadDobl_Random_Numbers;
with QuadDobl_Mathematical_Functions;
with Triple_Double_Numbers;              use Triple_Double_Numbers;

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

  procedure Main is
  begin
    new_line;
    put_line("Testing triple double arithmetic ...");
    Test_Basic_Arithmetic;
  end Main;

begin
  Main;
end ts_tridbl;
