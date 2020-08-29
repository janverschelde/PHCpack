with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;
with Multprec_Integer64_Numbers;
with Multprec_Integer64_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;

package body Test_Random_Numbers is

  function Modulus ( c : Standard_Complex_Numbers.Complex_Number ) 
                   return double_float is

    x : constant double_float := Standard_Complex_Numbers.REAL_PART(c);
    y : constant double_float := Standard_Complex_Numbers.IMAG_PART(c);

  begin
    return SQRT(x*x + y*y);
  end Modulus;

  procedure Random_Standard_Integer is

    n : natural32 := 0;
    low,upp,rndi : integer64 := 0;

  begin
    new_line;
    put_line("Testing the random generation of standard integer numbers.");
    new_line;
    put("Give lower bound : "); get(low);
    put("Give upper bound : "); get(upp);
    put("Give number of random integers : "); get(n);
    put_line("Generating random numbers : ");
    for i in 1..n loop
      rndi := Random(low,upp);
      put(rndi); new_line;
    end loop;
  end Random_Standard_Integer;

  procedure Random_Standard_Complex is

    n : natural32 := 0;
    m : natural32 := 1;
    use Standard_Complex_Numbers;
    rndc : Complex_Number;
    absv : double_float;

  begin
    new_line;
    put_line("Testing the random generation of standard complex numbers.");
    new_line;
    put("Give the number of randoms to be generated : "); get(n);
    put("Give the magnitude (1 for on unit circle) : "); get(m);
    if m = 1 then
      for k in 1..n loop
        rndc := Random1;
        put(" x = "); put(rndc); 
        absv := Modulus(rndc);
        put("   |x| = "); put(absv,3,3,3); new_line;
      end loop;
    else
      for k in 1..n loop
        rndc := Random_Magnitude(m);
        put(" x = "); put(rndc); 
        absv := Modulus(rndc);
        put("   |x| = "); put(absv,3,3,3); new_line;
      end loop;
    end if;
  end Random_Standard_Complex;

  procedure Random_Multprec_Natural is

    n,sz : natural32 := 0;
    rnd : Natural_Number;

  begin
    new_line;
    put_line("Testing the random generation of multi-precision naturals.");
    new_line;
    put("Give the size of the numbers : "); get(sz);
    put("Give number of random naturals : "); get(n);
    put_line("Generating random numbers : ");
    for i in 1..n loop
      rnd := Random(sz); put(rnd); new_line;
    end loop;
  end Random_Multprec_Natural;

  procedure Random_Multprec_Integer is

    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;

    n,sz : natural32 := 0;
    rnd : Integer_Number;

  begin
    new_line;
    put_line("Testing the random generation of multi-precision integers.");
    new_line;
    put("Give the size of the numbers : "); get(sz);
    put("Give number of random integers : "); get(n);
    put_line("Generating random numbers : ");
    for i in 1..n loop
      rnd := Random(sz); put(rnd); new_line;
    end loop;
  end Random_Multprec_Integer;

  procedure Random_Multprec_Integer64 is

    use Multprec_Integer64_Numbers,Multprec_Integer64_Numbers_io;

    n,sz : natural32 := 0;
    rnd : Integer_Number;

  begin
    new_line;
    put_line("Testing the random generation of multi-precision integers.");
    new_line;
    put("Give the size of the numbers : "); get(sz);
    put("Give number of random integers : "); get(n);
    put_line("Generating random numbers : ");
    for i in 1..n loop
      rnd := Random(sz); put(rnd); new_line;
    end loop;
  end Random_Multprec_Integer64;

  procedure Random_Multprec_Floating is

    n,sz : natural32 := 0;
    rnd : Floating_Number;

  begin
    new_line;
    put_line("Testing the random generation of multi-precision floats.");
    new_line;
    put("Give the size of the numbers : "); get(sz);
    put("Give number of random floating numbers : "); get(n);
    put_line("Generating random numbers : ");
    for i in 1..n loop
      rnd := Random(sz); put(rnd); new_line;
    end loop;
  end Random_Multprec_Floating;

  procedure Random_Multprec_Complex is

    n,sz : natural32 := 0;
    use Multprec_Complex_Numbers;
    rnd : Complex_Number;

  begin
    new_line;
    put_line("Testing the random generation of multi-precision complex");
    new_line;
    put("Give the size of the numbers : "); get(sz);
    put("Give number of random complex numbers : "); get(n);
    put_line("Generating random numbers : ");
    for i in 1..n loop
      rnd := Random(sz); put(rnd); new_line;
    end loop;
  end Random_Multprec_Complex;

  procedure Main is

    ans,lng : character;

  begin
    new_line;
    put_line("Testing the Random Number Generators");
    loop
      new_line;
      put_line("Choose one of the following :");
      put_line("  0. exit the program.");
      put_line("  1. standard integer numbers.");
      put_line("  2. standard complex numbers.");
      put_line("  3. multi-precision natural numbers.");
      put_line("  4. multi-precision integer numbers.");
      put_line("  5. multi-precision floating numbers.");
      put_line("  6. multi-precision complex numbers.");
      put("Make your choice (0,1,2,3,4,5 or 6) : ");
      Ask_Alternative(ans,"0123456");
      exit when (ans = '0');
      if ans = '4' then
        new_line;
        put("Use 64-bit arithmetic ? (y/n) ");
        Ask_Yes_or_No(lng);
      end if;
      if lng = 'y' then
        case ans is
          when '4' => Random_Multprec_Integer64;
          when others => null;
        end case;
      else
        case ans is
          when '1' => Random_Standard_Integer;
          when '2' => Random_Standard_Complex;
          when '3' => Random_Multprec_Natural;
          when '4' => Random_Multprec_Integer;
          when '5' => Random_Multprec_Floating;
          when '6' => Random_Multprec_Complex;
          when others => null;
        end case;
      end if;
    end loop;
  end Main;

end Test_Random_Numbers;
