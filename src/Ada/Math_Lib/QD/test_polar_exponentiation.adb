with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;
with Standard_Mathematical_Functions;
with Standard_Random_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Floating_Constants;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Double_Double_Constants;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Quad_Double_Constants;
with Multprec_DoblDobl_Convertors;       use Multprec_DoblDobl_Convertors;
with Multprec_QuadDobl_Convertors;       use Multprec_QuadDobl_Convertors;
with Standard_Complex_Exponentiation;    use Standard_Complex_Exponentiation;

package body Test_Polar_Exponentiation is

  procedure Test_Standard_DivModTwoPi is

    twopi : constant double_float := 2.0*Standard_Mathematical_Functions.Pi;
    x : double_float := 0.0;
    q : integer32;
    r,y : double_float;

  begin
    new_line;
    put("Give a double float x : "); get(x);
    put("-> your double float x : "); put(x); new_line;
    DivModTwoPi(x,q,r);
    put("-> quotient  div 2*Pi, q = "); put(q,1); new_line;
    put("-> remainder mod 2*Pi, r = "); put(r); new_line;
    y := twopi*double_float(q) + r;
    put("-> 2*Pi*q + r = "); put(y); new_line;
    put("->          x = "); put(x); new_line;
  end Test_Standard_DivModTwoPi;

  procedure Test_DoblDobl_DivModTwoPi is

    twopi : constant double_double := Double_Double_Constants.twopi;
    x : double_double := Create(0.0);
    q : Integer_Number;
    r,y,dd_q : double_double;

  begin
    new_line;
    put("Give a double double x : "); get(x);
    put("-> your double double x : "); put(x); new_line;
    DivModTwoPi(x,q,r);
    put("-> quotient  div 2*Pi, q = "); put(q); new_line;
    put("-> remainder mod 2*Pi, r = "); put(r); new_line;
    dd_q := to_double_double(q);
    y := twopi*dd_q + r;
    put("-> 2*Pi*q + r = "); put(y); new_line;
    put("->          x = "); put(x); new_line;
  end Test_DoblDobl_DivModTwoPi;

  procedure Test_QuadDobl_DivModTwoPi is

    twopi : constant quad_double := Quad_Double_Constants.twopi;
    x : quad_double := Create(0.0);
    q : Integer_Number;
    r,y,qd_q : quad_double;

  begin
    new_line;
    put("2*Pi = "); put(twopi); new_line;
    new_line;
    put("Give a quad double x : "); get(x);
    put("-> your quad double x : "); put(x); new_line;
    DivModTwoPi(x,q,r);
    put("-> quotient  div 2*Pi, q = "); put(q); new_line;
    put("-> remainder mod 2*Pi, r = "); put(r); new_line;
    qd_q := to_quad_double(q);
    y := twopi*qd_q + r;
    put("-> 2*Pi*q + r = "); put(y); new_line;
    put("->          x = "); put(x); new_line;
  end Test_QuadDobl_DivModTwoPi;

  procedure Test_Multprec_DivModTwoPi is

    x,r,y,fq,twopi : Floating_Number;
    dp : natural32;
    q : Integer_Number;
    ans : character;

  begin
    new_line;
    put("Give a floating number x : "); get(x);
    put("-> your floating number x : "); put(x); new_line;
    dp := Multprec_Floating_Numbers_io.Character_Size(x);
    put("-> number of decimal places : "); put(dp,1); new_line;
    put("Change number of decimal places ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put("Give new number of decimal places : "); get(dp);
    end if;
    DivModTwoPi(x,dp,q,r);
    put("-> quotient  div 2*Pi, q = "); put(q); new_line;
    put("-> remainder mod 2*Pi, r = "); put(r); new_line;
    fq := Create(q);
    twopi := Multprec_Floating_Constants.TwoPi(dp);
    y := twopi*fq + r;
    put("-> 2*Pi*q + r = "); put(y); new_line;
    put("->          x = "); put(x); new_line;
  end Test_Multprec_DivModTwoPi;

  function Binary_Exponentiation 
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    res,acc,tmp : Complex_Number;
    abse : natural32;

  begin
    if e = 0 then
      res := Create(1.0);
    else
      if e > 0
       then abse := natural32(e);
       else abse := natural32(-e);
      end if;
      res := x; acc := Create(1.0);
      if abse > 1 then          -- use binary exponentiation
        while abse > 0 loop
          if abse mod 2 = 1
           then tmp := res; Mul(acc,tmp); -- to avoid compiler warning
          end if;
          abse := abse/2;
          if abse > 0
           then tmp := res; Mul(res,tmp); -- to avoid compiler warning
          end if;
        end loop;
      else
        acc := res;
      end if;
      if e < 0
       then res := 1.0/acc;          -- compute reciprocal
       else res := acc;
      end if;
    end if;
    return res;
  end Binary_Exponentiation;

  procedure Test_Standard_Exponentiation is

    ans : character;
    x,y,z,w : Complex_Number;
    v : double_float;
    e : integer32 := 0;
    u : boolean;

  begin
    new_line;
    put("Exponentiation of a random number ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      x := Standard_Random_Numbers.Random1;
      put("-> a random number x = "); 
    else
      put("Give a complex number x : "); get(x);
      put("-> your number x = "); 
    end if;
    u := (ans = 'y');
    put(x); new_line;
    v := Standard_Complex_Numbers_Polar.Radius(x);
    put("  |x| = "); put(v); new_line;
    put("Give an exponent e : "); get(e);
    y := x**integer(e);
    v := Standard_Complex_Numbers_Polar.Radius(y);
    put("-> x^"); put(e,1); put(" = "); put(y); new_line;
    put("  |x^"); put(e,1); put("| = "); put(v); new_line;
    if u
     then z := Standard_Complex_Numbers_Polar.Polar_Exponentiation_of_Unit(x,e);
     else z := Standard_Complex_Numbers_Polar.Polar_Exponentiation(x,e);
    end if;
    v := Standard_Complex_Numbers_Polar.Radius(z);
    put_line("using the polar representation of x :");
    put("-> x^"); put(e,1); put(" = "); put(z); new_line;
    put("  |x^"); put(e,1); put("| = "); put(v); new_line;
    if u 
     then z := Polar_Exponentiation_ModTwoPi_of_Unit(x,e);
     else z := Polar_Exponentiation_ModTwoPi(x,e);
    end if;
    v := Standard_Complex_Numbers_Polar.Radius(z);
    put_line("using the polar representation of x (mod 2*Pi) :");
    put("-> x^"); put(e,1); put(" = "); put(z); new_line;
    put("  |x^"); put(e,1); put("| = "); put(v); new_line;
    w := Binary_Exponentiation(x,e);
    v := Standard_Complex_Numbers_Polar.Radius(w);
    put_line("using binary exponentiation :");
    put("-> x^"); put(e,1); put(" = "); put(w); new_line;
    put("  |x^"); put(e,1); put("| = "); put(v); new_line;
  end Test_Standard_Exponentiation;

  procedure Test_Multprec_Exponentiation is

    x : constant Complex_Number := Standard_Random_Numbers.Random1;
    v : double_float;
    e : Integer_Number;
    y : Complex_Number;
    e32 : integer32;

  begin
    new_line;
    put("-> x = "); put(x); new_line;
    v := Standard_Complex_Numbers_Polar.Radius(x);
    put("  |x| = "); put(v); new_line;
    put("Give an exponent e : "); get(e);
    y := Polar_Exponentiation_ModTwoPi_of_Unit(x,e);
    put("-> x^"); put(e); put(" = "); put(y); new_line;
    v := Standard_Complex_Numbers_Polar.Radius(y);
    put("  |x^"); put(e); put("| = "); put(v); new_line;
    if Decimal_Places(e) > 9 then
      put_line("Exponent is too large to compare to 32-bit exponentiation.");
    else
      e32 := Create(e);
      put_line("Comparing to 32-bit exponentiation...");
     -- y := Polar_Exponentiation_ModTwoPi_of_Unit(x,e32);
     -- y := Binary_Exponentiation(x,e32);
      y := Standard_Complex_Numbers_Polar.Polar_Exponentiation_of_Unit(x,e32);
      put("-> x^"); put(e32,1); put(" = "); put(y); new_line;
      v := Standard_Complex_Numbers_Polar.Radius(y);
      put("  |x^"); put(e32,1); put("| = "); put(v); new_line;
    end if;
  end Test_Multprec_Exponentiation;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to text exponentiation operations :");
    put_line("  1. test division modulo two times Pi with standard doubles;");
    put_line("  2. test division modulo two times Pi with double doubles;");
    put_line("  3. test division modulo two times Pi with quad doubles;");
    put_line("  4. test multiprecision division modulo two times Pi;");
    put_line("  5. exponentiation of a complex number with standard exponent;");
    put_line("  6. exponentiation of a complex with multiprecision exponent.");
    put("Type 1, 2, 3, 4, 5, or 6 to choose : ");
    Ask_Alternative(ans,"123456");
    case ans is
      when '1' => Test_Standard_DivModTwoPi;
      when '2' => Test_DoblDobl_DivModTwoPi;
      when '3' => Test_QuadDobl_DivModTwoPi;
      when '4' => Test_Multprec_DivModTwoPi;
      when '5' => Test_Standard_Exponentiation;
      when '6' => Test_Multprec_Exponentiation;
      when others => null;
    end case;
  end Main;

end Test_Polar_Exponentiation;
