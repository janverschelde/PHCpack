with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with TripDobl_Random_Numbers;
with TripDobl_Mathematical_Functions;    use TripDobl_Mathematical_Functions;

package body Test_TripDobl_Functions is

  procedure Test_SQRT is

    x,y,z,err : triple_double;
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
        x := TripDobl_Random_Numbers.Random;
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

    x,a,b,c : triple_double;
    ans : character;

  begin
    new_line;
    put_line("Testing the radius with a (6, 8, 10) ...");
    loop
      x := TripDobl_Random_Numbers.Random;
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

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test triple double mathematical functions :");
    put_line("  1. test square root of triple doubles");
    put_line("  2. test radius for triple doubles");
    put("Type 1 or 2 to select a test : "); Ask_Alternative(ans,"12");
    if ans = '1'
     then Test_SQRT;
     else Test_Radius;
    end if;
  end Main;

end Test_TripDobl_Functions;
