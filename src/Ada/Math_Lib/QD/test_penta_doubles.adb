with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Penta_Double_Numbers_io;            use Penta_Double_Numbers_io;
with Penta_Double_Constants;

package body Test_Penta_Doubles is

  procedure write ( x : penta_double ) is
  begin
    put("  thumb  : "); put(thumb_part(x),2,17,3); new_line;
    put("  index  : "); put(index_part(x),2,17,3); new_line;
    put("  middle : "); put(middle_part(x),2,17,3); new_line;
    put("  ring   : "); put(ring_part(x),2,17,3); new_line;
    put("  pink   : "); put(pink_part(x),2,17,3); new_line;
  end Write;

  function random return penta_double is

    res : penta_double;
    first : constant double_float := Standard_Random_Numbers.Random; 
    second : double_float := Standard_Random_Numbers.Random; 
    eps : constant double_float := 2.0**(-52);
    multiplier : double_float := eps;

  begin
    res := create(first);
    res := res + eps*second;
    for k in 3..5 loop
      multiplier := eps*multiplier;
      second := Standard_Random_Numbers.Random;
      res := res + multiplier*second;
    end loop;
    return res;
  end random;

  procedure Test_Addition_and_Subtraction is

    x : constant penta_double := random;
    y : constant penta_double := random;
    z : constant penta_double := x + y;
    v : constant penta_double := z - y;

  begin
   put_line("All parts of a random octo double x :"); Write(x);
   put_line("All parts of a random octo double y :"); Write(y);
   put_line("All parts of x + y :"); Write(z);
   put_line("All parts of (x + y) - y :"); Write(v);
  end Test_Addition_and_Subtraction;

  procedure Test_Multiplication_and_Division is

    x : constant penta_double := random;
    y : constant penta_double := random;
    z : constant penta_double := x*y;
    v : constant penta_double := z/y;

  begin
   put_line("All parts of a random penta double x :"); Write(x);
   put_line("All parts of a random penta double y :"); Write(y);
   put_line("All parts of x * y :"); Write(z);
   put_line("All parts of (x * y) / y :"); Write(v);
  end Test_Multiplication_and_Division;

  procedure Test_Read is

    sqrt2 : constant string
          := "1.4142135623730950488016887242096980785696718753769480"
           & "731766797379907324784621070";
    x,r : penta_double;
    two : constant penta_double := create(2.0);
    fail : boolean;

  begin
    read(sqrt2,x,fail);
    new_line;
    if fail then
      put_line("The read procedure reports failure!");
    else
      put_line("All parts of the sqrt(2) read :"); Write(x);
      r := x*x - 2.0;
      put_line("All parts of x*x - 2.0 : "); Write(r);
      r := x*x - two;
      put_line("All parts of x*x - two : "); Write(r);
    end if;
  end Test_Read;

  procedure Test_io is

    x,y : penta_double;
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

  --   >>> from sympy import evalf, sqrt
  --   >>> s2 = sqrt(2).evalf(128)
  --   >>> s2
  --   1.4142135623730950488016887242096980785696718753769480
  --   731766797379907324784621070

    n,x,y,z,e,a : penta_double;
    max_steps : constant natural32 := 8;
    sqrt2 : constant string
          := "1.4142135623730950488016887242096980785696718753769480"
           & "731766797379907324784621070";
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

  procedure Test_pd_eps is

    one : constant penta_double := create(1.0);
    eps : constant double_float := Penta_Double_Constants.pd_eps;
    one_plus_pd_eps : constant penta_double := one + eps;
    inc : constant double_float := (eps/2.0);
    one_plus_pd_eps_half : constant penta_double := one + inc;

  begin
    new_line;
    put("    pd_eps   :"); put(eps); new_line;
    put_line("1 + pd_eps   : "); put(one_plus_pd_eps,79); new_line;
    put_line("1 + pd_eps/2 : "); put(one_plus_pd_eps_half,79); new_line;
  end Test_pd_eps;

  procedure Log_exp_of_Pi is

    pd_pi : constant penta_double
          := create( 3.14159265358979312E+00, 1.22464679914735321E-16,
                    -2.99476980971833967E-33, 1.11245422086336528E-49,
                     5.67223197964031574E-66);
    exp_of_pi,log_of_exp_of_pi : penta_double;
    log_of_pi,exp_of_log_of_pi : penta_double;
    ans : character;
    x,expx,logx,logexpx,explogx : penta_double;

  begin
    new_line;
    put_line("testing log(exp(pi)) and exp(log(pi)) ...");
    exp_of_pi := exp(pd_pi);
    put("     exp(pi) : "); put(exp_of_pi); new_line;
    log_of_exp_of_pi := log(exp_of_pi);
    put("log(exp(pi)) : "); put(log_of_exp_of_pi); new_line;
    put("        pi   : "); put(pd_pi); new_line;
    put("  difference : "); put(log_of_exp_of_pi - pd_pi,3); new_line;
    log_of_pi := log(pd_pi);
    put("     log(pi) : "); put(log_of_pi); new_line;
    exp_of_log_of_pi := exp(log_of_pi);
    put("exp(log(pi)) : "); put(exp_of_log_of_pi); new_line;
    put("        pi   : "); put(pd_pi); new_line;
    put("  difference : "); put(exp_of_log_of_pi - pd_pi,3); new_line;
    loop
      put("Test a random number ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      x := abs(random);
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
    put_line("Testing penta double arithmetic ...");
    put_line("  1. test addition and subtraction");
    put_line("  2. test multiplication and division");
    put_line("  3. test reading from a string");
    put_line("  4. input and output");
    put_line("  5. Newton's method for sqrt(2)");
    put_line("  6. test the value of pd_eps");
    put_line("  7. test log(exp(pi)) = pi = exp(log(pi))");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select a test : ");
    Ask_Alternative(ans,"1234567");
    case ans is
      when '1' => Test_Addition_and_Subtraction;
      when '2' => Test_Multiplication_and_Division;
      when '3' => Test_Read;
      when '4' => Test_io;
      when '5' => Test_sqrt2;
      when '6' => Test_pd_eps;
      when '7' => Log_exp_of_Pi;
      when others => null;
    end case;
  end Main;

end Test_Penta_Doubles;
