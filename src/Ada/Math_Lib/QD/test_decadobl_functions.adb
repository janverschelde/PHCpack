with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Deca_Double_Numbers_io;             use Deca_Double_Numbers_io;
with DecaDobl_Random_Numbers;
with DecaDobl_Mathematical_Functions;    use DecaDobl_Mathematical_Functions;

package body Test_DecaDobl_Functions is

  procedure Test_SQRT is

    x,y,z,err : deca_double;
    ans : character;

  begin
    new_line;
    put_line("Testing the square root function ...");
    loop
      put("Generate a random number ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans /= 'y' then
        put("Give x : "); get(x);
      else
        x := DecaDobl_Random_Numbers.Random;
        if is_negative(x)
         then Min(x);
        end if;
      end if;
      put("        x : "); put(x); new_line;
      y := SQRT(x);
      put("  sqrt(x) : "); put(y); new_line;
      z := y*y;
      put("sqrt(x)^2 : "); put(z); new_line;
      err := x - z;
      put("x - sqrt(x)^2 : "); put(err,3); new_line;
      put("More tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_SQRT;

  procedure Test_Radius is

  -- (6, 8, 10) is a Pythagorean triple, 6^2 + 8^2 = 10^2,
  -- which allows for a quick visual check on the accuracy
  -- of the radius of 6*x and 8*x for any random number x.

    x,a,b,c : deca_double;
    ans : character;

  begin
    new_line;
    put_line("Testing the radius with a (6, 8, 10) ...");
    loop
      x := DecaDobl_Random_Numbers.Random;
      if is_negative(x)
       then Min(x);
      end if;
      put("   x : "); put(x); new_line;
      put("10*x : "); put(10.0*x); new_line;
      a := 6.0*x;
      b := 8.0*x;
      c := Radius(a,b);
      put_line("See if R = Radius(6*x,8*x) = 10*x ...");
      put("   R : "); put(c); new_line;
      put("More tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Radius;

  procedure Test_SINCOS is

    x,s,c,y,t,a,inv_s,inv_c,err : deca_double;
    one : constant deca_double := create(1.0);
    ans : character;

  begin
    new_line;
    put_line("testing sin and cos function ...");
    loop
      put("generate random number ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then x := DecaDobl_Random_Numbers.Random;
       else put("give x : "); get(x);
      end if;
                   put("     x : "); put(x); new_line;
      s := SIN(x); put("sin(x) : "); put(s); new_line;
      c := COS(x); put("cos(x) : "); put(c); new_line;
      y := s*s + c*c;
      put("sin(x)^2 + cos(x)^2 : "); put(y); new_line;
      err := y - one;
      put("sin(x)^2 + cos(x)^2 - 1 : "); put(err,3); new_line;
      put("sin(x)/cos(x) : "); put(s/c); new_line;
      t := tan(x);
      put("       tan(x) : "); put(t); new_line;
      a := arctan(t);
      put("    arctan(tan(x)) : "); put(a); new_line;
      err := x - a;
      put("x - arctan(tan(x)) : "); put(err,3); new_line;
      inv_s := arcsin(s);
      put("arcsin(sin(x)) : "); put(inv_s); new_line;
      err := x - inv_s;
      put("x - arcsin(sin(x)) : "); put(err,3); new_line;
      inv_c := arccos(c);
      put("arccos(cos(x)) : "); put(inv_c); new_line;
      err := x - inv_c;
      put("x - arccos(cos(x)) : "); put(err,3); new_line;
      put("Test more ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_SINCOS;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test deca double mathematical functions :");
    put_line("  1. test square root of deca doubles");
    put_line("  2. test radius for deca doubles");
    put_line("  3. test cos, sin, tan, and their inverses");
    put("Type 1, 2, or 3 to select a test : "); Ask_Alternative(ans,"123");
    case ans is
      when '1' => Test_SQRT;
      when '2' => Test_Radius;
      when '3' => Test_SINCOS;
      when others => null;
    end case;
  end Main;

end Test_DecaDobl_Functions;
