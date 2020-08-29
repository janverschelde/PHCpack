with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;

package body Test_Mathematical_Functions is

  function C_COS ( x : double_float ) return double_float is

    function cos ( x : double_float ) return double_float;
    pragma import(C,cos); -- pragma interface(C,cos);

  begin
    return cos(x);
  end C_COS;

  function C_SIN ( x : double_float ) return double_float is

    function sin ( x : double_float ) return double_float;
    pragma import(C,sin); -- pragma interface(C,sin);

  begin
    return sin(x);
  end C_SIN;

  procedure Read ( f : in out Floating_Number; name : in string ) is

    n : natural32 := 0;

  begin
    put("Give " & name & " : "); get(f);
    put("Current size is "); put(Size_Fraction(f),1);
    put(".  Give expansion factor : "); get(n);
    if n > 0
     then Expand(f,n);
    end if;
  end Read;

  procedure Test_Standard_COS_and_SIN is

    x,cx,sx : double_float := 0.0;
    ans : character;
    
  begin
    new_line;
    put_line("Testing whether cos^2(x) + sin^2(x) = 1 holds.");
    loop
      new_line;
      put("Give x : "); get(x);
     -- cx := COS(x);  cx := cx*cx;
     -- sx := SIN(x);  sx := sx*sx;
      cx := C_COS(x);  cx := cx*cx;
      sx := C_SIN(x);  sx := sx*sx;
      put("cos^2(x) + sin^2(x) = "); put(cx+sx); new_line;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Standard_COS_and_SIN;

  procedure Test_Standard_SQRT is

    x,y,diff : double_float := 0.0;

  begin
    put("Give x : "); get(x);
    y := SQRT(x);
    put("SQRT(x) : "); put(y); new_line;
    y := y*y;
    put("(SQRT(x))**2 : "); put(y); new_line;
    diff := x - y;
    put("x - (SQRT(x))**2 : "); put(diff); new_line;
  end Test_Standard_SQRT;

  procedure Test_Multprec_SQRT
              ( x : in Floating_Number; output : in boolean ) is

    y,y2,diff : Floating_Number;

  begin
    y := SQRT(x);
    if output
     then put("SQRT(x) : "); put(y); new_line;
    end if;
    y2 := y*y;
    if output
     then put("SQRT(x))**2 : "); put(y2); new_line;
    end if;
    diff := x - y2;
    put("x - (SQRT(x))**2 : "); put(diff); new_line;
    Clear(y); Clear(y2); Clear(diff);
  end Test_Multprec_SQRT;

  procedure Interactive_Test_Multprec_SQRT is

    x : Floating_Number;
    ans : character;

  begin
    loop
      new_line;
      Read(x,"x");
      Test_Multprec_SQRT(x,true);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Interactive_Test_Multprec_SQRT;

  procedure Random_Test_Multprec_SQRT is

    sz,nb : natural32 := 0;
    x : Floating_Number;
    timer : timing_widget;

  begin
    new_line;
    put("Give size of the numbers : "); get(sz);
    put("Give the number of samples : "); get(nb);
    tstart(timer);
    for i in 1..nb loop
      x := Random(sz);
      Test_Multprec_SQRT(x,false);
      Clear(x);
    end loop;
    tstop(timer);
    put("tested "); put(nb,1); put_line(" SQRT computations");
    new_line;
    print_times(Standard_Output,timer,"testing random SQRT");
  end Random_Test_Multprec_SQRT;

  procedure ln_one_plus_x ( x,eps : in Floating_Number;
                            nit : out natural32; max : in natural32;
                            y : out Floating_Number; fail : out boolean ) is

    sum,acc,quot : Floating_Number;   
    nk : natural32 := 0;
    fk : Floating_Number := Create(natural32(1));
    plus_sign : boolean := false;

  begin
    Copy(x,sum);
    Copy(x,acc);
    for i in 1..max loop
      nk := nk + 1;
      Mul(acc,x);
      Add(fk,1.0);
      quot := acc/fk;
      if plus_sign
       then Add(sum,quot);
       else Sub(sum,quot);
      end if;
      if quot > 0.0
       then fail := (quot > eps);
       else fail := (quot < eps);
      end if;
      Clear(quot);
      exit when not fail;
      plus_sign := not plus_sign;
    end loop;
    Clear(acc);
    nit := nk;
    y := sum;
  end ln_one_plus_x;

  procedure Test_Series is

    x,y,eps : Floating_Number;
    dbl_eps : double_float;
    steps,n,size : natural32 := 0;
    nb_digits : natural32 := 0;
    fail : boolean;

  begin
    new_line;
    put_line("Testing the convergence of the series for ln(1+x)...");
    new_line;
    put("Give a floating-point number : "); get(x);
    put("How many digits in the answer should be correct ? "); get(nb_digits);
    put("What is the maximal number of steps ? "); get(n);
    size := Decimal_to_Size(nb_digits);
    Set_Size(x,size);
    dbl_eps := 10.0**(-integer(nb_digits));
    eps := Create(dbl_eps);
    ln_one_plus_x(x,eps,steps,n,y,fail);
    put("Computed value is "); put(y); new_line;
    if fail then
      put("failed to reach the desired accuracy");
      put(eps,3); put(" in "); put(n,1); put_line(" steps.");
    else
      put("reached the desired accuracy");
      put(eps,3); put(" in "); put(steps,1); put_line(" steps.");
    end if;
  end Test_Series;

  function Read_Base return character is

    ans : character;

  begin
    put_line("The following bases are available :");
    put_line("  e : natural logarithm;");
    put_line("  2 : binary logarithm;");
    put_line("  A : decimal logarithm.");
    put("Type e, 2, or A to select the base : ");
    Ask_Alternative(ans,"e2A");
    return ans;
  end Read_Base;

  procedure Test_Multprec_Logarithm is

    x,y : Floating_Number;
    dbl_x,dbl_y,check : double_float;
    nb_digits,size : natural32 := 0;
    ans,base : character;

  begin
    new_line;
    put_line("Testing the multi-precision logarithm for bases e, 2, and 10.");
    new_line;
    loop
      put("Give a floating-point number : "); get(x);
      put("How many digits in the answer should be correct ? ");
      get(nb_digits);
      size := Decimal_to_Size(nb_digits);
      Set_Size(x,size);
      dbl_x := Round(x);
      base := Read_Base;
      case base is
        when 'e' => y := LN(x);    dbl_y := LN(dbl_x);
                                   check := EXP(dbl_y);
        when '2' => y := LOG2(x);  dbl_y := LOG2(dbl_x);
                                   check := 2.0**dbl_y;
        when 'A' => y := LOG10(x); dbl_y := LOG10(dbl_x);
                                   check := 10.0**dbl_y;
        when others => null;
      end case;
      put(" multi-precision logarithm : "); put(y); new_line;
      put("double precision logarithm : "); put(dbl_y); new_line;
      put("sanity check, with inverse : "); put(check); new_line;
      Clear(x); Clear(y);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Multprec_Logarithm;

  procedure Test_Multprec_Exponential is

    x,y : Floating_Number;
    dbl_x,dbl_y,check : double_float;
    nb_digits,size : natural32 := 0;
    ans,base : character;

  begin
    new_line;
    put_line("Testing multi-precision exponentiation for bases e, 2, and 10.");
    new_line;
    loop
      put("Give a floating-point number : "); get(x);
      put("How many digits in the answer should be correct ? ");
      get(nb_digits);
      size := Decimal_to_Size(nb_digits);
      Set_Size(x,size);
      dbl_x := Round(x);
      base := Read_Base;
      case base is
        when '2' => y := 2.0**x;  dbl_y := 2.0**dbl_x;  check := LOG2(dbl_y);
        when 'e' => y := EXP(x);  dbl_y := EXP(dbl_x);  check := LN(dbl_y);
        when 'A' => y := 10.0**x; dbl_y := 10.0**dbl_x; check := LOG10(dbl_y);
        when others => null;
      end case;
      put(" multi-precision exponentiation : "); put(y); new_line;
      put("double precision exponentiation : "); put(dbl_y); new_line;
      put(" sanity check, with the inverse : "); put(check); new_line;
      Clear(x); Clear(y);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Multprec_Exponential;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing some mathematical functions.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. Exit this program.");
      put_line("  1. COS/SIN for standard floating-point numbers.");
      put_line("  2. SQRT for standard floating-point numbers.");
      put_line("  3. SQRT of given multi-precision floating numbers.");
      put_line("  4. SQRT of random multi-precision floating numbers.");
      put_line("  5. test convergence of series for ln(1+x).");
      put_line("  6. multi-precision logarithm for bases e, 2, and 10.");
      put_line("  7. multi-precision exponentiation for bases e, 2, and 10.");
      put("Type 0,1,2,3,4,5,6 or 7 to make your selection : ");
      Ask_Alternative(ans,"01234567");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Standard_COS_and_SIN;
        when '2' => Test_Standard_SQRT;
        when '3' => Interactive_Test_Multprec_SQRT;
        when '4' => Random_Test_Multprec_SQRT;
        when '5' => Test_Series;
        when '6' => Test_Multprec_Logarithm;
        when '7' => Test_Multprec_Exponential;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Mathematical_Functions;
