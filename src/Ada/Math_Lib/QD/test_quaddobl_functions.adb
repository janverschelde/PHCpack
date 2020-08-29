with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Random_Numbers;
with QuadDobl_Mathematical_Functions;    use QuadDobl_Mathematical_Functions;

package body Test_QuadDobl_Functions is

  procedure Test_SQRT is

    x,y,z,err : quad_double;
    ans : character;

  begin
    new_line;
    put_line("testing the square root function ...");
    loop
      put("Generate a random number ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans /= 'y' then
        put("Give x : "); get(x);
      else
        x := QuadDobl_Random_Numbers.Random;
        if is_negative(x)
         then Min(x);
        end if;
      end if;
      put("        x : "); put(x); new_line;
      y := sqrt(x);
      put("  sqrt(x) : "); put(y); new_line;
      z := sqr(y);
      put("sqrt(x)^2 : "); put(z); new_line;
      err := x - z;
      put("x - sqrt(x)^2 : "); put(err,3); new_line;
      put("More tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_SQRT;

  procedure Test_SINCOS is

    x,s,c,y,t,a,inv_s,inv_c,err : quad_double;
    one : constant quad_double := create(1.0);
    ans : character;

  begin
    new_line;
    put_line("testing sin and cos function ...");
    loop
      put("generate random number ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then x := QuadDobl_Random_Numbers.Random;
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
    put_line("MENU for testing quad double mathematical functions :");
    put_line("  1. test square root of quad doubles");
    put_line("  2. test cos, sin, tan and their inverses for quad doubles");
    put("Type 1 or 2 : "); Ask_Alternative(ans,"12");
    if ans = '1'
     then Test_SQRT;
     else Test_SINCOS;
    end if;
  end Main;

end Test_QuadDobl_Functions;
