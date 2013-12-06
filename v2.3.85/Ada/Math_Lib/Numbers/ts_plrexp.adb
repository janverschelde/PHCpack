with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;
with Standard_Random_Numbers;

procedure ts_plrexp is

-- DESCRIPTION :
--   Test on exponentiation of complex numbers
--   via polar representation.

  twopi : constant double_float := 2.0*Standard_Mathematical_Functions.Pi;

  procedure DivModTwoPi ( x : in double_float;
                          q : out integer32; r : out double_float ) is

  -- DESCRIPTION :
  --   Computes quotient q and remainder r after division of 2*Pi,
  --   so on return we have x = 2*Pi*q + r.
  --   If x < 2*Pi, then q = 0 and r = x on return.

    fq : double_float;
    iq : integer32;

  begin
    if x < twopi then
      q := 0; r := x;
    else
      fq := x/twopi;
      iq := integer32(fq);
      r := x - double_float(iq)*twopi;
      if r < 0.0
       then r := r + twopi; iq := iq - 1;
      end if;
      q := iq;
    end if;
  end DivModTwoPi;

  procedure Test_DivModTwoPi is

  -- DESCRIPTION :
  --   Interactive test to compute division and remainder
  --   of any integer number by 2 times Pi.

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
  end Test_DivModTwoPi;

  function Polar_Exponentiation_ModTwoPi
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns x^e via polar coordinates of x, computing the remainder
  --   of the exponent modulo 2*Pi.

    r : constant double_float := Standard_Complex_Numbers_Polar.Radius(x);
    a : constant double_float := Standard_Complex_Numbers_Polar.Angle(x);
    s : constant double_float := r**integer(e);
    f : constant double_float := double_float(e)*a;
    q : integer32;
    b,re,im : double_float;
    res : Complex_Number;

  begin
    DivModTwoPi(f,q,b);
    re := s*COS(b);
    im := s*SIN(b);
    res := Create(re,im);
    return res;
  end Polar_Exponentiation_ModTwoPi;

  function Binary_Exponentiation 
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

  -- DESCRIPTION :
  --   Uses binary exponentiation to compute x^e.

    res,acc : Complex_Number;
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
           then Mul(acc,res);
          end if;
          abse := abse/2;
          if abse > 0
           then Mul(res,res);
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

  procedure Test_Exponentiation is

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
    z := Polar_Exponentiation_ModTwoPi(x,e);
    v := Standard_Complex_Numbers_Polar.Radius(z);
    put_line("using the polar representation of x (mod 2*Pi) :");
    put("-> x^"); put(e,1); put(" = "); put(z); new_line;
    put("  |x^"); put(e,1); put("| = "); put(v); new_line;
    w := Binary_Exponentiation(x,e);
    v := Standard_Complex_Numbers_Polar.Radius(w);
    put_line("using binary exponentiation :");
    put("-> x^"); put(e,1); put(" = "); put(w); new_line;
    put("  |x^"); put(e,1); put("| = "); put(v); new_line;
  end Test_Exponentiation;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to text exponentiation operations :");
    put_line("  1. test division modulo two times Pi;");
    put_line("  2. test exponentiation of a complex number.");
    put("Type 1 or 2 to choose : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then Test_DivModTwoPi;
     else Test_Exponentiation;
    end if;
  end Main;

begin
  Main;
end ts_plrexp;
