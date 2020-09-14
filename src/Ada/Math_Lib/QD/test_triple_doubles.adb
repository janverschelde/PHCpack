with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Random_Numbers;            use QuadDobl_Random_Numbers;
with QuadDobl_Mathematical_Functions;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with Triple_Double_Constants;
with TripDobl_Random_Numbers;

package body Test_Triple_Doubles is

  function to_triple_double ( x : quad_double ) return triple_double is

    res : constant triple_double
        := create(hihi_part(x),lohi_part(x),hilo_part(x));

  begin
    return res;
  end to_triple_double;

  function create ( x : triple_double ) return quad_double is

    res : constant quad_double
        := create(hi_part(x),mi_part(x),lo_part(x),0.0);

  begin
    return res;
  end create;

  procedure Test_Basic_Arithmetic is

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

  procedure Test_sqrt2 is

    n,x,y,z,e,a : triple_double;
    max_steps : constant natural32 := 7;
    sqrt2 : constant string
          := "1.4142135623730950488016887242096980785696718753769";
    fail : boolean;

  begin
    n := Create(2.0); Copy(n,x);
    new_line;
    put_line("running Newton's method for sqrt(2) ...");
    read(sqrt2,y,fail);
    if fail
     then put_line("reading value for sqrt2 from string failed!");
    end if;
    put(" sqrt2: "); put(y); new_line;
    put("step 0: "); put(x); new_line;
    for i in 1..max_steps loop
      z := x*x;
      z := z+n;
      z := z/x;
      z := 0.5*z;
      put("step "); put(i,1); put(": "); put(z); new_line;
      copy(z,x);
      e := x - y;
      a := abs(e);
      put("  error : "); put(a,3); new_line;
    end loop;
  end Test_sqrt2;

  procedure Test_td_eps is

    one : constant triple_double := create(1.0);
    eps : constant double_float := Triple_Double_Constants.td_eps/2.0;
    one_plus_td_eps : constant triple_double := one + eps;
    inc : constant double_float := (eps/2.0);
    one_plus_td_eps_half : constant triple_double := one + inc;

  begin
    new_line;
    put("    td_eps   :"); put(eps); new_line;
    put_line("1 + td_eps   : "); put(one_plus_td_eps,47); new_line;
    put_line("1 + td_eps/2 : "); put(one_plus_td_eps_half,47); new_line;
  end Test_td_eps;

  procedure Log_exp_of_Pi is

    td_pi : constant triple_double
          := create( 3.141592653589793116e+00, 1.224646799147353207e-16,
                    -2.994769809718339666e-33);
    exp_of_pi,log_of_exp_of_pi : triple_double;
    log_of_pi,exp_of_log_of_pi : triple_double;
    ans : character;
    x,expx,logx,logexpx,explogx : triple_double;

  begin
    new_line;
    put_line("testing log(exp(pi)) and exp(log(pi)) ...");
    exp_of_pi := exp(td_pi);
    put("     exp(pi) : "); put(exp_of_pi); new_line;
    log_of_exp_of_pi := log(exp_of_pi);
    put("log(exp(pi)) : "); put(log_of_exp_of_pi); new_line;
    put("        pi   : "); put(td_pi); new_line;
    put("  difference : "); put(log_of_exp_of_pi - td_pi,3); new_line;
    log_of_pi := log(td_pi);
    put("     log(pi) : "); put(log_of_pi); new_line;
    exp_of_log_of_pi := exp(log_of_pi);
    put("exp(log(pi)) : "); put(exp_of_log_of_pi); new_line;
    put("        pi   : "); put(td_pi); new_line;
    put("  difference : "); put(exp_of_log_of_pi - td_pi,3); new_line;
    loop
      put("Test a random number ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      x := abs(TripDobl_Random_Numbers.Random);
      expx := exp(x); logexpx := log(expx);
      put("log(exp(x)) : "); put(logexpx); new_line;
      put("        x   : "); put(x); new_line;
      put(" difference : "); put(logexpx - x,3); new_line;
      logx := log(x); explogx := exp(logx);
      put("exp(log(x)) : "); put(explogx); new_line;
      put("        x   : "); put(x); new_line;
      put(" difference : "); put(explogx - x,3); new_line;
    end loop;
  end Log_exp_of_Pi;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing triple double arithmetic ...");
    put_line("  1. basic arithmetic");
    put_line("  2. test reading from string");
    put_line("  3. input and output");
    put_line("  4. Newton's method for sqrt(2)");
    put_line("  5. test the value of td_eps");
    put_line("  6. test log(exp(pi)) = pi = exp(log(pi))");
    put("Type 1, 2, 3, 4, 5, or 6 to select a test : ");
    Ask_Alternative(ans,"123456");
    case ans is
      when '1' => Test_Basic_Arithmetic;
      when '2' => Test_Read;
      when '3' => Test_io;
      when '4' => Test_sqrt2;
      when '5' => Test_td_eps;
      when '6' => Log_exp_of_Pi;
      when others => null;
    end case;
  end Main;

end Test_Triple_Doubles;
