with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Quad_Double_Constants;
with QuadDobl_Random_Numbers;

package body Test_Quad_Doubles is

  procedure Write ( q : in quad_double ) is
  begin
    put("  q hihi = "); put(hihi_part(q));
    put("  q lohi = "); put(lohi_part(q)); new_line;
    put("  q hilo = "); put(hilo_part(q));
    put("  q lolo = "); put(lolo_part(q)); new_line;
  end Write;

  procedure Basic_Test is

    i : integer32 := 0;
    q : quad_double;

  begin
    put("Give an integer : "); get(i);
    q := Create(i);
    put_line("The integer as quad double :"); Write(q);
  end Basic_Test;

  procedure Test_io is

    x,y : quad_double;
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

    pi_0 : constant double_float :=  3.141592653589793116e+00; -- qd_pi[0]
    pi_1 : constant double_float :=  1.224646799147353207e-16; -- qd_pi[1]
    pi_2 : constant double_float := -2.994769809718339666e-33; -- qd_pi[2]
    pi_3 : constant double_float :=  1.112454220863365282e-49; -- qd_pi[3] 
    qd_pi : constant quad_double := Create(pi_0,pi_1,pi_2,pi_3);
    e_0 : constant double_float :=  2.718281828459045091e+00;  -- qd_e[0] 
    e_1 : constant double_float :=  1.445646891729250158e-16;  -- qd_e[1]
    e_2 : constant double_float := -2.127717108038176765e-33;  -- qd_e[2] 
    e_3 : constant double_float :=  1.515630159841218954e-49;  -- qd_e[3]
    qd_e : constant quad_double := Create(e_0,e_1,e_2,e_3);
    pi_and_e,pi_and_e_minus_pi : quad_double;

  begin
    new_line;
    put_line("testing pi + e - pi ...");
    put_line("quad double representation of pi: "); Write(qd_pi);
    put_line("quad double representation of e: "); Write(qd_e);
    pi_and_e := qd_pi + qd_e;
    put_line("Pi + e :"); Write(pi_and_e);
    pi_and_e_minus_pi := pi_and_e - qd_pi;
    put_line("Pi + e - Pi :"); Write(pi_and_e_minus_pi);
  end Add_Sub_of_Pi_e;

  procedure Div_sqr_of_Pi is

    pi_0 : constant double_float :=  3.141592653589793116e+00; -- qd_pi[0]
    pi_1 : constant double_float :=  1.224646799147353207e-16; -- qd_pi[1]
    pi_2 : constant double_float := -2.994769809718339666e-33; -- qd_pi[2]
    pi_3 : constant double_float :=  1.112454220863365282e-49; -- qd_pi[3] 
    qd_pi : constant quad_double := Create(pi_0,pi_1,pi_2,pi_3);
    sqr_of_pi,div_of_sqr_of_pi : quad_double;

  begin
    new_line;
    put_line("testing sqr(pi)/pi ...");
    put_line("quad double representation of pi :"); Write(qd_pi);
    sqr_of_pi := sqr(qd_pi);
    put_line("pi^2 : "); Write(sqr_of_pi);
    div_of_sqr_of_pi := sqr_of_pi/qd_pi;
    put_line("pi^2/pi : "); Write(div_of_sqr_of_pi);
  end Div_sqr_of_Pi;

  procedure Log_exp_of_Pi is

    pi_0 : constant double_float :=  3.141592653589793116e+00; -- qd_pi[0]
    pi_1 : constant double_float :=  1.224646799147353207e-16; -- qd_pi[1]
    pi_2 : constant double_float := -2.994769809718339666e-33; -- qd_pi[2]
    pi_3 : constant double_float :=  1.112454220863365282e-49; -- qd_pi[3] 
    qd_pi : constant quad_double := Create(pi_0,pi_1,pi_2,pi_3);
    exp_of_pi,log_of_exp_of_pi : quad_double;

  begin
    new_line;
    put_line("testing log(exp(pi)) ...");
    put_line("quad double representation of pi :"); Write(qd_pi);
    exp_of_pi := exp(qd_pi);
    put_line("exp(pi) :"); Write(exp_of_pi);
    log_of_exp_of_pi := log(exp_of_pi);
    put_line("log(exp(pi)) :"); Write(log_of_exp_of_pi);
  end Log_exp_of_Pi;

  procedure my_sqrt is

    n,x,y,z,e,a : quad_double;
    max_steps : constant natural32 := 9;
    sqrt2 : constant string
      := "1.414213562373095048801688724209698078569671875376948073176679738\0";
    fail : boolean;

  begin
    n := Create(2.0); Copy(n,x);
    new_line;
    put_line("running Newton's method for sqrt(2) ...");
    read(sqrt2,y,fail);
    if fail
     then put_line("reading value for sqrt2 from string failed!");
    end if;
    put(" sqrt2: "); write(y,64); new_line;
    Write(y);
    put("step 0: "); write(x,64); new_line;
    for i in 1..max_steps loop
      z := x*x;
      Add(z,n);
      Div(z,x);
      Mul(z,0.5);
      put("step "); put(i,1); put(": "); write(z,64);
      copy(z,x);
      e := x - y;
      a := abs(e);
      put("  error : "); write(a,3); new_line;
    end loop;
  end my_sqrt;

  procedure Test_Random is

    m : natural32 := 0;
    r : quad_double;

  begin
    put("Give the magnitude : "); get(m);
    if m = 1
     then r := QuadDobl_Random_Numbers.Random;
     else r := QuadDobl_Random_Numbers.Random_Magnitude(m);
    end if;
    put("the random number : "); put(r); new_line;
  end Test_Random;

  procedure Test_qd_eps is

    one : constant quad_double := create(1.0);
    eps : constant double_float := Quad_Double_Constants.qd_eps;
    one_plus_qd_eps : constant quad_double := one + eps;
    inc : constant double_float := (eps/2.0);
    one_plus_qd_eps_half : constant quad_double := one + inc;

  begin
    new_line;
    put("    qd_eps   :"); put(eps); new_line;
    put_line("1 + qd_eps   : "); put(one_plus_qd_eps,63); new_line;
    put_line("1 + qd_eps/2 : "); put(one_plus_qd_eps_half,63); new_line;
  end Test_qd_eps;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test operations on quad double arithmetic :");
    put_line("  0. do a basic test");
    put_line("  1. test input/output of quad double numbers");
    put_line("  2. run test on add/sub, div/sqr, log/exp, and a sqrt");
    put_line("  3. generate a random quad double");
    put_line("  4. test the value of qd_eps");
    put("Type 0, 1, 2, 3, or 4 to choose : "); Ask_Alternative(ans,"01234");
    case ans is
      when '0' => Basic_Test;
      when '1' => Test_io;
      when '2' => Add_Sub_of_Pi_e;
                  Div_sqr_of_Pi;
                  Log_exp_of_Pi;
                  my_sqrt;
      when '3' => Test_Random;
      when '4' => Test_qd_eps;
      when others => null;
    end case;
  end Main;

end Test_Quad_Doubles;
