with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard64_Common_Divisors;         use Standard64_Common_Divisors;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Multprec_Common_Divisors;           use Multprec_Common_Divisors;

package body Test_Greatest_Common_Divisors is

  procedure Interactive_Test_Standard_GCD is

    a,b,k,l,d : integer64 := 0;
    ans : character;

  begin
    loop
      put("Give a : "); get(a);
      put("Give b : "); get(b);
      d := lcm(a,b);
      put("lcm("); put(a,1); put(','); put(b,1); put(") = "); put(d,1);
      new_line;
      gcd(a,b,k,l,d);
      put("gcd("); put(a,1); put(','); put(b,1); put(") = "); put(d,1);
      new_line;
      put("  k = "); put(k,1); new_line;
      put("  l = "); put(l,1); new_line;
      put(" k*a + l*b = "); put(k*a + l*b,1); new_line;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /='y';
    end loop;
  end Interactive_Test_Standard_GCD;

  procedure Test_Standard_GCD ( low,upp : in integer64 ) is

    a,b,k,l,d1,d2 : integer64;

    procedure Report_Bug is
    begin
      put(" a : "); put(a,1); new_line;
      put(" b : "); put(b,1); new_line;
    end Report_Bug;

  begin
    a := Random(low,upp);
    b := Random(low,upp);
    d1 := gcd(a,b);
    gcd(a,b,k,l,d2);
    put(k,1); put("*"); put(a,1); put("+");
    put(l,1); put("*"); put(b,1); put("="); put(d2,1);
    if d1 = d2 and k*a + l*b = d2
     then put("  okay"); new_line;
     else put("  Bug!"); Report_Bug;
    end if;
  end Test_Standard_GCD;

  procedure Random_Test_Standard_GCD is

    nb : natural32 := 0;
    low,upp : integer64 := 0;

  begin
    put("Give the number of tests : "); get(nb);
    put("Give a lower bound for the numbers : "); get(low);
    put("Give an upper bound for the numbers : "); get(upp);
    for i in 1..nb loop
      Test_Standard_GCD(low,upp);
    end loop;
  end Random_Test_Standard_GCD;

  procedure Interactive_Test_Multprec_GCD is

    a,b,k,l,d : Integer_Number;
    ans : character;

  begin
    loop
      put("Give a : "); get(a);
      put("Give b : "); get(b);
      d := lcm(a,b);
      put("lcm("); put(a); put(','); put(b); put(") = "); put(d);
      new_line;
      gcd(a,b,k,l,d);
      put("gcd("); put(a); put(','); put(b); put(") = "); put(d);
      new_line;
      put("  k = "); put(k); new_line;
      put("  l = "); put(l); new_line;
      put(" k*a + l*b = "); put(k*a + l*b); new_line;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /='y';
    end loop;
  end Interactive_Test_Multprec_GCD;

  procedure Test_Multprec_GCD ( sz1,sz2 : in natural32 ) is

    a,b,k,l,d1,d2,acc1,acc2 : Integer_Number;

    procedure Report_Bug is
    begin
      put(" a : "); put(a); new_line;
      put(" b : "); put(b); new_line;
    end Report_Bug;

  begin
    a := Random(sz1);
    b := Random(sz2);
    d1 := gcd(a,b);
    gcd(a,b,k,l,d2);
    put(k); put("*"); put(a); put("+");
    put(l); put("*"); put(b); put("="); put(d2);
    acc1 := k*a;
    acc2 := l*b;
    Add(acc1,acc2);
    if Equal(d1,d2) and Equal(acc1,d2)
     then put("  okay"); new_line;
     else put("  Bug!"); Report_Bug;
    end if;
    Clear(a); Clear(b); Clear(d1);
    Clear(k); Clear(l); Clear(d2);
    Clear(acc1); Clear(acc2);
  end Test_Multprec_GCD;

  procedure Random_Test_Multprec_GCD is

    nb,sz1,sz2 : natural32 := 0;

  begin
    put("Give the number of tests : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    for i in 1..nb loop
      Test_Multprec_GCD(sz1,sz2);
    end loop;
  end Random_Test_Multprec_GCD;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of gcd-computations.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program.");
      put_line("  1. interactive gcd of standard integer numbers.");
      put_line("  2. gcd of randomly generated standard integer numbers.");
      put_line("  3. interactive gcd of multi-precision integer numbers.");
      put_line("  4. gcd of random multi-precision integer numbers.");
      put("Type 0,1,2,3, or 4 to select your choice : "); get(ans);
      exit when ans = '0';
      case ans is
        when '1' => Interactive_Test_Standard_GCD;
        when '2' => Random_Test_Standard_GCD;
        when '3' => Interactive_Test_Multprec_GCD;
        when '4' => Random_Test_Multprec_GCD;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Greatest_Common_Divisors;
