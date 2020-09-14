with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with Octo_Double_Constants;

package body Test_Octo_Doubles is

  function random return octo_double is

    res : octo_double;
    first : constant double_float := Standard_Random_Numbers.Random; 
    second : double_float := Standard_Random_Numbers.Random; 
    eps : constant double_float := 2.0**(-52);
    multiplier : double_float := eps;

  begin
    res := create(first);
    res := res + eps*second;
    for k in 3..8 loop
      multiplier := eps*multiplier;
      second := Standard_Random_Numbers.Random;
      res := res + multiplier*second;
    end loop;
    return res;
  end random;

  procedure Write ( x : in octo_double ) is
  begin
    put("  hihihi : "); put(hihihi_part(x),2,17,3); new_line;
    put("  lohihi : "); put(lohihi_part(x),2,17,3); new_line;
    put("  hilohi : "); put(hilohi_part(x),2,17,3); new_line;
    put("  lolohi : "); put(lolohi_part(x),2,17,3); new_line;
    put("  hihilo : "); put(hihilo_part(x),2,17,3); new_line;
    put("  lohilo : "); put(lohilo_part(x),2,17,3); new_line;
    put("  hilolo : "); put(hilolo_part(x),2,17,3); new_line;
    put("  lololo : "); put(lololo_part(x),2,17,3); new_line;
  end Write;

  procedure Test_Add_and_Subtract is

    x : constant octo_double := random;
    y : constant octo_double := random;
    z : constant octo_double := x + y;
    v : constant octo_double := z - y;

  begin
   put_line("All parts of a random octo double x :"); Write(x);
   put_line("All parts of a random octo double y :"); Write(y);
   put_line("All parts of x + y :"); Write(z);
   put_line("All parts of (x + y) - y :"); Write(v);
  end Test_Add_and_Subtract;

  procedure Test_Multiplication_and_Division is

    x : constant octo_double := random;
    y : constant octo_double := random;
    z : constant octo_double := x*y;
    v : constant octo_double := z/y;

  begin
   put_line("All parts of a random octo double x :"); Write(x);
   put_line("All parts of a random octo double y :"); Write(y);
   put_line("All parts of x * y :"); Write(z);
   put_line("All parts of (x * y) / y :"); Write(v);
  end Test_Multiplication_and_Division;

  procedure Test_Read is

  --   >>> from sympy import evalf, sqrt
  --   >>> s2 = sqrt(2).evalf(128)
  --   >>> s2
  --   1.414213562373095048801688724209698078569671875376948073176679737
  --   9907324784621070388503875343276415727350138462309122970249248361
    
    sqrt2 : constant string
      := "1.414213562373095048801688724209698078569671875376948073176679737"
       & "9907324784621070388503875343276415727350138462309122970249248361";
    x,r : octo_double;
    two : constant octo_double := create(2.0);
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

  procedure Log10log2exp1_doubles is

  -- >>> from sympy import log, evalf
  -- >>> x = log(10)
  -- >>> x.evalf(128)
  -- 2.3025850929940456840179914546843642076011014886287729760333279009
  -- 675726096773524802359972050895982983419677840422862486334095255
  -- >>> from sympy import exp, evalf
  -- >>> x = exp(1)
  -- >>> x.evalf(128)
  -- 2.7182818284590452353602874713526624977572470936999595749669676277
  -- 240766303535475945713821785251664274274663919320030599218174136
  -- >>> from sympy import log, evalf
  -- >>> x = log(2)
  -- >>> x.evalf(128)
  -- 0.6931471805599453094172321214581765680755001343602552541206800094
  -- 9339362196969471560586332699641868754200148102057068573368552024

    log10 : constant string
      := "2.3025850929940456840179914546843642076011014886287729760333279009"
       & "675726096773524802359972050895982983419677840422862486334095255";
    log2 : constant string
      := "0.6931471805599453094172321214581765680755001343602552541206800094"
       & "9339362196969471560586332699641868754200148102057068573368552024";
    exp1 : constant string
      := "2.7182818284590452353602874713526624977572470936999595749669676277"
       & "240766303535475945713821785251664274274663919320030599218174136";

    x : octo_double;
    fail : boolean;

  begin
    read(log10,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the log(10) read :"); Write(x);
    end if;
    read(log2,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the log(2) read :"); Write(x);
    end if;
    read(exp1,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the exp(1) read :"); Write(x);
    end if;
  end Log10log2exp1_doubles;

  procedure inverse_factorials is

  -- >>> from sympy import factorial, evalf
  -- >>> f = [1/factorial(k) for k in range(3,15+3)]
  -- >>> for x in f: print(x.evalf(128))

    f0 : constant string
       := "0.1666666666666666666666666666666666666666666666666666666666666666"
        & "6666666666666666666666666666666666666666666666666666666666666667";
    f1 : constant string
       := "0.0416666666666666666666666666666666666666666666666666666666666666"
        & "66666666666666666666666666666666666666666666666666666666666666667";
    f2 : constant string
       := "0.0083333333333333333333333333333333333333333333333333333333333333"
        & "33333333333333333333333333333333333333333333333333333333333333333"
        & "3";
    f3 : constant string
       := "0.0013888888888888888888888888888888888888888888888888888888888888"
        & "88888888888888888888888888888888888888888888888888888888888888888"
        & "9";
    f4 : constant string
       := "0.0001984126984126984126984126984126984126984126984126984126984126"
        & "98412698412698412698412698412698412698412698412698412698412698412"
        & "70";
    f5 : constant string
       := "0.0000248015873015873015873015873015873015873015873015873015873015"
        & "87301587301587301587301587301587301587301587301587301587301587301"
        & "587";
    f6 : constant string
       := "0.0000027557319223985890652557319223985890652557319223985890652557"
        & "31922398589065255731922398589065255731922398589065255731922398589"
        & "0653";
    f7 : constant string
       := "0.0000002755731922398589065255731922398589065255731922398589065255"
        & "73192239858906525573192239858906525573192239858906525573192239858"
        & "90653";
    f8 : constant string
       := "0.0000000250521083854417187750521083854417187750521083854417187750"
        & "52108385441718775052108385441718775052108385441718775052108385441"
        & "718775";
    f9 : constant string
       := "0.0000000020876756987868098979210090321201432312543423654534765645"
        & "87675698786809897921009032120143231254342365453476564587675698786"
        & "8098979";
    f10 : constant string
        := "0.0000000001605904383682161459939237717015494793272571050348828126"
         & "60590438368216145993923771701549479327257105034882812660590438368"
         & "21614599";
    f11 : constant string
        := "0.0000000000114707455977297247138516979786821056662326503596344866"
         & "18613602740586867570994555121539248523375507502491629475756459883"
         & "444010428";
    f12 : constant string
        := "0.0000000000007647163731819816475901131985788070444155100239756324"
         & "41240906849372457838066303674769283234891700500166108631717097325"
         & "56293402854";
    f13 : constant string
        := "0.0000000000000477947733238738529743820749111754402759693764984770"
         & "27577556678085778614879143979673080202180731281260381789482318582"
         & "847683376784";
    f14 : constant string
        := "0.0000000000000028114572543455207631989455830103200162334927352045"
         & "31033973922240339918522302587039592953069454781250610693498959916"
         & "6380990221638";
    x : octo_double;
    fail : boolean;
    ans : character;

  begin
    read(f0,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(0) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f1,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(1) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f2,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(2) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f3,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(3) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f4,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(4) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f5,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(5) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f6,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(6) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f7,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(7) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f8,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(8) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f9,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(9) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f10,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(10) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f11,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(11) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f12,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(12) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f13,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(13) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f14,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(14) read :"); Write(x);
    end if;
  end inverse_factorials;

  procedure Test_io is

    x,y : octo_double;
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

    n,x,y,z,e,a : octo_double;
    max_steps : constant natural32 := 8;
    sqrt2 : constant string
      := "1.414213562373095048801688724209698078569671875376948073176679737"
       & "9907324784621070388503875343276415727350138462309122970249248361";
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

  procedure Test_od_eps is

    one : constant octo_double := create(1.0);
    eps : constant double_float := Octo_Double_Constants.od_eps;
    one_plus_od_eps : constant octo_double := one + eps;
    inc : constant double_float := (eps/2.0);
    one_plus_od_eps_half : constant octo_double := one + inc;

  begin
    new_line;
    put("    od_eps   :"); put(eps); new_line;
    put_line("1 + od_eps   : "); put(one_plus_od_eps,127); new_line;
    put_line("1 + od_eps/2 : "); put(one_plus_od_eps_half,127); new_line;
  end Test_od_eps;

  procedure Log_exp_of_Pi is

    od_pi : constant octo_double
          := create( 3.14159265358979312E+00, 1.22464679914735321E-16,
                    -2.99476980971833967E-33, 1.11245422086336528E-49,
                     5.67223197964031574E-66, 1.74498621613524860E-83,
                     6.02937273224953984E-100, 1.91012354687998999E-116);
    exp_of_pi,log_of_exp_of_pi : octo_double;
    log_of_pi,exp_of_log_of_pi : octo_double;
    ans : character;
    x,expx,logx,logexpx,explogx : octo_double;

  begin
    new_line;
    put_line("testing log(exp(pi)) and exp(log(pi)) ...");
    exp_of_pi := exp(od_pi);
    put("     exp(pi) : "); put(exp_of_pi); new_line;
    log_of_exp_of_pi := log(exp_of_pi);
    put("log(exp(pi)) : "); put(log_of_exp_of_pi); new_line;
    put("        pi   : "); put(od_pi); new_line;
    put("  difference : "); put(log_of_exp_of_pi - od_pi,3); new_line;
    log_of_pi := log(od_pi);
    put("     log(pi) : "); put(log_of_pi); new_line;
    exp_of_log_of_pi := exp(log_of_pi);
    put("exp(log(pi)) : "); put(exp_of_log_of_pi); new_line;
    put("        pi   : "); put(od_pi); new_line;
    put("  difference : "); put(exp_of_log_of_pi - od_pi,3); new_line;
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
    put_line("Testing octo double arithmetic ...");
    put_line("  1. test addition and subtraction");
    put_line("  2. test multiplication and division");
    put_line("  3. test reading from a string");
    put_line("  4. write 8 leading doubles of log(10), log(2), exp(1)");
    put_line("  5. write 8 leading double for inverse factorials");
    put_line("  6. input and output");
    put_line("  7. Newton's method for sqrt(2)");
    put_line("  8. test the value of od_eps");
    put_line("  9. test log(exp(pi)) = pi = exp(log(pi))");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, or 9 to select a test : ");
    Ask_Alternative(ans,"123456789");
    case ans is
      when '1' => Test_Add_and_Subtract;
      when '2' => Test_Multiplication_and_Division;
      when '3' => Test_Read;
      when '4' => Log10log2exp1_doubles;
      when '5' => inverse_factorials;
      when '6' => Test_io;
      when '7' => Test_sqrt2;
      when '8' => Test_od_eps;
      when '9' => Log_exp_of_Pi;
      when others => null;
    end case;
  end Main;

end Test_Octo_Doubles;
