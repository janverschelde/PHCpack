with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Double_Double_Constants;
with DoblDobl_Random_Numbers;

package body Test_Double_Doubles is

  procedure character_arithmetic is

    i,hc : integer32 := 0;
    c,ans,ch : character;

  begin
    new_line;
    loop
      put("give an integer : "); get(i);
      ch := character'val(i);
      put("character'val("); put(i,1); put(") = ");
      put(ch); new_line;
      hc := character'pos(ch);
      put("character('"); put(ch); put("') = ");
      put(hc,1); new_line;
      put("give a character : "); get(c);
      hc := character'pos(c);
      put("character'pos('"); put(c); put("') = ");
      put(hc,1); new_line;
      ch := character'val(hc);
      put("character'val("); put(hc,1); put(") = ");
      put(ch); new_line;
      put("more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end character_arithmetic;

  procedure Write ( d : in double_double ) is
  begin
    put("d hi = "); put(hi_part(d));
    put("  d lo = "); put(lo_part(d)); new_line;
  end Write;

  procedure Basic_Test is

    i,n : integer32 := 0;
    d,dn : double_double;

  begin
    new_line;
    put("Give an integer : "); get(i);
    d := Create(i);
    put_line("The integer as double double :"); Write(d);
    put("Give exponent n of 2 : "); get(n);
    dn := ldexp(d,integer(n));
    put_line("ldexp(d,n) :"); Write(dn);
  end Basic_Test;

  procedure Test_io is

    x,y : double_double;
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

  procedure Add_Sub_of_Pi_e is

    pi_lo : constant double_float := 3.141592653589793116e+00; -- pi hi word
    pi_hi : constant double_float := 1.224646799147353207e-16; -- pi lo word
    e_lo : constant double_float := 2.718281828459045091e+00;  -- e hi word
    e_hi : constant double_float := 1.445646891729250158e-16;  -- e lo word
    dd_pi,dd_e,pi_and_e,pi_and_e_minus_pi : double_double;

  begin
    new_line;
    put_line("testing pi + e - pi ...");
    dd_pi := Create(pi_lo,pi_hi);
    dd_e := Create(e_lo,e_hi);
    put_line("double double representation of pi :"); Write(dd_pi);
    put_line("double double representation of e :"); Write(dd_e);
    pi_and_e := dd_pi + dd_e;
    put_line("pi + e :"); Write(pi_and_e);
    pi_and_e_minus_pi := pi_and_e - dd_pi;
    put_line("pi + e - pi :"); Write(pi_and_e_minus_pi);
  end Add_Sub_of_Pi_e;

  procedure Div_sqr_of_Pi is

    pi_lo : constant double_float := 3.141592653589793116e+00; -- pi hi word
    pi_hi : constant double_float := 1.224646799147353207e-16; -- pi lo word
    dd_pi,sqr_of_pi,div_of_sqr_of_pi : double_double;

  begin
    new_line;
    put_line("testing sqr(pi)/pi ...");
    dd_pi := Create(pi_lo,pi_hi);
    put_line("double double representation of pi :"); Write(dd_pi);
    sqr_of_pi := sqr(dd_pi);
    put_line("pi^2 : "); Write(sqr_of_pi);
    div_of_sqr_of_pi := sqr_of_pi/dd_pi;
    put_line("pi^2/pi : "); Write(div_of_sqr_of_pi);
  end Div_sqr_of_Pi;

  procedure Log_exp_of_Pi is

    pi_lo : constant double_float := 3.141592653589793116e+00; -- pi hi word
    pi_hi : constant double_float := 1.224646799147353207e-16; -- pi lo word
    dd_pi,exp_of_pi,log_of_exp_of_pi : double_double;

  begin
    new_line;
    put_line("testing log(exp(pi)) ...");
    dd_pi := Create(pi_lo,pi_hi);
    put_line("double double representation of pi :"); Write(dd_pi);
    exp_of_pi := exp(dd_pi);
    put_line("exp(pi) :"); Write(exp_of_pi);
    log_of_exp_of_pi := log(exp_of_pi);
    put_line("log(exp(pi)) :"); Write(log_of_exp_of_pi);
  end Log_exp_of_Pi;

  procedure my_sqrt is

    n,x,y,z,e,a : double_double;
    max_steps : constant natural32 := 9;
    sqrt2 : constant string := "1.4142135623730950488016887242097\0";
    fail : boolean;

  begin
    n := Create(2.0); Copy(n,x);
    new_line;
    put_line("running Newton's method for sqrt(2) ...");
    read(sqrt2,y,fail);
    if fail
     then put_line("reading value for sqrt2 from string failed!");
    end if;
    Write(y);
    put("step 0: "); write(x,32); new_line;
    for i in 1..max_steps loop
      z := x*x;
      Add(z,n);
      Div(z,x);
      Mul(z,0.5);
      put("step "); put(i,1); put(": "); write(z,32);
      copy(z,x);
      e := x - y;
      a := abs(e);
      put("  error : "); write(a,3); new_line;
    end loop;
  end my_sqrt;

  procedure Test_Random is

    m : natural32 := 0;
    r : double_double;

  begin
    put("Give the magnitude : "); get(m);
    if m = 1
     then r := DoblDobl_Random_Numbers.Random;
     else r := DoblDobl_Random_Numbers.Random_Magnitude(m);
    end if;
    put("the random number : "); put(r); new_line;
  end Test_Random;

  procedure Test_dd_eps is

    one : constant double_double := create(1.0);
    eps : constant double_float := Double_Double_Constants.dd_eps;
    one_plus_dd_eps : constant double_double := one + eps;
    inc : constant double_float := (eps/2.0);
    one_plus_dd_eps_half : constant double_double := one + inc;

  begin
    new_line;
    put("    dd_eps   :"); put(eps); new_line;
    put("1 + dd_eps   : "); put(one_plus_dd_eps,31); new_line;
    put("1 + dd_eps/2 : "); put(one_plus_dd_eps_half,31); new_line;
  end Test_dd_eps;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test operations on double double numbers :");
    put_line("  0. test character arithmetic");
    put_line("  1. do a basic interactive test on ldexp");
    put_line("  2. test input/output of double double numbers");
    put_line("  3. run add/sub, sqr/div, log/exp, and a sqrt");
    put_line("  4. generate a random double double");
    put_line("  5. test the value of dd_eps");
    put("Type 0, 1, 2, 3, 4, or 5 to choose : ");
    Ask_Alternative(ans,"012345");
    case ans is
      when '0' => character_arithmetic;
      when '1' => Basic_Test;
      when '2' => Test_io;
      when '3' => Add_Sub_of_Pi_e;
                  Div_sqr_of_Pi;
                  Log_exp_of_Pi;
                  my_sqrt;
      when '4' => Test_Random;
      when '5' => Test_dd_eps;
      when others => null;
    end case;
  end Main;

end Test_Double_Doubles;
