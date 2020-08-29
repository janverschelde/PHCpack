with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Multprec_Natural_Coefficients;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Multprec_Floating_Numbers_io;
with Multprec_Floating64_Numbers_io;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;

package body Test_Floating_Numbers is

  procedure Read ( f : in out Multprec_Floating_Numbers.Floating_Number;
                   name : in string ) is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;
    n : natural32 := 0;

  begin
    put("Give " & name & " : "); get(f);
    put("Current size is "); put(Size_Fraction(f),1);
    put(".  Give expansion factor : "); get(n);
    if n > 0
     then Expand(f,n);
    end if;
  end Read;

  procedure Read ( f : in out Multprec_Floating64_Numbers.Floating_Number;
                   name : in string ) is

    use Multprec_Floating64_Numbers,Multprec_Floating64_Numbers_io;
    n : natural32 := 0;

  begin
    put("Give " & name & " : "); get(f);
    put("Current size is "); put(Size_Fraction(f),1);
    put(".  Give expansion factor : "); get(n);
    if n > 0
     then Expand(f,n);
    end if;
  end Read;

  procedure Formatted_Output
              ( f : in Multprec_Floating_Numbers.Floating_Number ) is

    use Multprec_Floating_Numbers_io;
    fore,aft,exp : natural32 := 0;
    
  begin
    put("Give the number of places before the decimal point : ");
    get(fore);
    put("Give the number of places after the decimal point : ");
    get(aft);
    put("Give the number of places of the exponent : ");
    get(exp);
    put("-> formatted : "); put(f,fore,aft,exp); new_line;
  end Formatted_Output;

  procedure Formatted64_Output
              ( f : in Multprec_Floating64_Numbers.Floating_Number ) is

    use Multprec_Floating64_Numbers_io;
    fore,aft,exp : natural32 := 0;
    
  begin
    put("Give the number of places before the decimal point : ");
    get(fore);
    put("Give the number of places after the decimal point : ");  
    get(aft);
    put("Give the number of places of the exponent : ");
    get(exp);
    put("-> formatted : "); put(f,fore,aft,exp); new_line;
  end Formatted64_Output;

  procedure Test_io is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    f,abf : Floating_Number;
    ans,c : character;
    nc : natural32;

  begin
    put_line("Testing input/output for multiprecision floating numbers.");
    loop
      c := ' ';
      put("Give a floating number : "); get(f,c);
      put("-> your floating number : "); put(f); new_line;
      put_line("   was terminated by the character " & c);
      abf := AbsVal(f);
      put("-> its absolute value : "); put(abf); new_line;
      put("-> its fraction : "); put(Fraction(f)); new_line;
      put("-> #decimal places in fraction : ");
      put(Decimal_Places_Fraction(f),1); new_line;
      put("-> its exponent : "); put(Exponent(f)); new_line;
      put("-> #decimal places in exponent : ");
      put(Decimal_Places_Exponent(f),1); new_line;
      nc := Character_Size(f);
      put("#characters in string representation : ");
      put(nc,1); new_line;
      declare
        sf : string(1..integer(nc));
      begin
        put(sf,f);
        put("the floating-point number : "); put(f); new_line;
        put("its string representation : "); put(sf); new_line;
      end;
      put("Do you want formatted output ? (y/n) "); get(ans);
      if ans = 'y'
       then Formatted_Output(f);
      end if;
      Clear(f); Clear(abf);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_io;

  procedure Test64_io is

    use Multprec_Floating64_Numbers,Multprec_Floating64_Numbers_io;

    f,abf : Floating_Number;
    ans,c : character;
    nc : natural32;

  begin
    put_line("Testing input/output for multiprecision floating numbers.");
    loop
      c := ' ';
      put("Give a floating number : "); get(f,c);
      put("-> your floating number : "); put(f); new_line;
      put_line("   was terminated by the character " & c);
      abf := AbsVal(f);
      put("-> its absolute value : "); put(abf); new_line;
      put("-> #decimal places in fraction : ");
      put(Decimal_Places_Fraction(f),1); new_line;
      put("-> #decimal places in exponent : ");
      put(Decimal_Places_Exponent(f),1); new_line;
      nc := Character_Size(f);
      put("#characters in string representation : ");
      put(nc,1); new_line;
      declare
        sf : string(1..integer(nc));
      begin
        put(sf,f);
        put("the floating-point number : "); put(f); new_line;
        put("its string representation : "); put(sf); new_line;
      end;
      put("Do you want formatted output ? (y/n) "); get(ans);
      if ans = 'y'
       then Formatted64_Output(f);
      end if;
      Clear(f); Clear(abf);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test64_io;
   
  procedure Test_Creation is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    f : Floating_Number;
    d,fd : double_float := 0.0;
    i : integer32 := 0;
    ans : character;

  begin
    put_line("Testing the creation of multiprecision floating numbers.");
    loop
      put("Give an integer : "); get(i);
      put("-> your integer : "); put(i,1); new_line;
      f := Create(i);
      put("-> as floating number : "); put(f); new_line;
      put("Give a standard float : "); get(d);
      put("-> your float : "); put(d); new_line;
      f := Create(d);
      put("-> as floating number : "); put(f); new_line;
      fd := Round(f);
      put("-> rounded as standard float : "); put(fd); new_line;
      if d = fd
       then put_line("Creation/Rounding test is successful.");
       else put_line("Difference up to working precision ?");
            put("d - Round(Create(d)) : "); put(f-fd); new_line;
      end if;
      put("Give a floating number : "); get(f);
      put("-> your floating number : "); put(f); new_line;
      d := Round(f);
      put("-> rounded as float     :"); put(d); new_line;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Creation;

  procedure Test64_Creation is

    use Multprec_Floating64_Numbers,Multprec_Floating64_Numbers_io;

    f : Floating_Number;
    d,fd : double_float := 0.0;
    i : integer64 := 0;
    ans : character;

  begin
    put_line("Testing the creation of multiprecision floating numbers.");
    loop
      put("Give an integer : "); get(i);
      put("-> your integer : "); put(i,1); new_line;
      f := Create(i);
      put("-> as floating number : "); put(f); new_line;
      put("Give a standard float : "); get(d);
      put("-> your float : "); put(d); new_line;
      f := Create(d);
      put("-> as floating number : "); put(f); new_line;
      fd := Round(f);
      put("-> rounded as standard float : "); put(fd); new_line;
      if d = fd
       then put_line("Creation/Rounding test is successful.");
       else put_line("Difference up to working precision ?");
            put("d - Round(Create(d)) : "); put(f-fd); new_line;
      end if;
      put("Give a floating number : "); get(f);
      put("-> your floating number : "); put(f); new_line;
      d := Round(f);
      put("-> rounded as float     :"); put(d); new_line;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test64_Creation;

--  procedure Test_Compare ( f1 : in Floating_Number;
--                           f2 : in double_float ) is
--  begin
--    if Equal(f1,f2)
--     then put_line("The numbers are equal.");
--     else put_line("The numbers are different.");
--    end if;
--    if f1 < f2
--     then put_line("First number is less than second number.");
--     else put_line("First number is not less than second number.");
--    end if;
--    if f1 > f2
--     then put_line("First number is greater than second number.");
--     else put_line("First number is not greater than second number.");
--    end if;
--  end Test_Compare;

  procedure Test_Compare
              ( f1,f2 : in Multprec_Floating_Numbers.Floating_Number ) is

    use Multprec_Floating_Numbers;

  begin
    if Equal(f1,f2)
     then put_line("The numbers are equal.");
     else put_line("The numbers are different.");
    end if;
    if f1 < f2
     then put_line("First number is less than second number.");
     else put_line("First number is not less than second number.");
    end if;
    if f1 > f2
     then put_line("First number is greater than second number.");
     else put_line("First number is not greater than second number.");
    end if;
  end Test_Compare;

  procedure Test64_Compare
              ( f1,f2 : in Multprec_Floating64_Numbers.Floating_Number ) is

    use Multprec_Floating64_Numbers;

  begin
    if Equal(f1,f2)
     then put_line("The numbers are equal.");
     else put_line("The numbers are different.");
    end if;
    if f1 < f2
     then put_line("First number is less than second number.");
     else put_line("First number is not less than second number.");
    end if;
    if f1 > f2
     then put_line("First number is greater than second number.");
     else put_line("First number is not greater than second number.");
    end if;
  end Test64_Compare;

  procedure Zero_Test ( f : in Multprec_Floating_Numbers.Floating_Number ) is
  begin
    if Multprec_Floating_Numbers.Equal(f,0.0)
     then put_line(" equals zero.");
     else put_line(" is different from zero.");
    end if;
  end Zero_Test;

  procedure Zero_Test ( f : in Multprec_Floating64_Numbers.Floating_Number ) is
  begin
    if Multprec_Floating64_Numbers.Equal(f,0.0)
     then put_line(" equals zero.");
     else put_line(" is different from zero.");
    end if;
  end Zero_Test;

  procedure Test_Comparison is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    f1,f2 : Floating_Number;
   -- f2 : double_float;
    ans : character;

  begin
    put_line("Testing comparison/copying for multiprecision floats.");
    loop
      put("Give 1st number f1 : "); get(f1);
      put(" f1 : "); put(f1); 
      Zero_Test(f1);
      put("Give 2nd number f2 : "); get(f2);
      put(" f2 : "); put(f2);
      Zero_Test(f2);
      Test_Compare(f1,f2);
     -- Copy(f1,f2);
     -- put_line("After copy :");
     -- Test_Compare(f1,f2);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Comparison;

  procedure Test64_Comparison is

    use Multprec_Floating64_Numbers,Multprec_Floating64_Numbers_io;

    f1,f2 : Floating_Number;
   -- f2 : double_float;
    ans : character;

  begin
    put_line("Testing comparison/copying for multiprecision floats.");
    loop
      put("Give 1st number f1 : "); get(f1);
      put(" f1 : "); put(f1); 
      Zero_Test(f1);
      put("Give 2nd number f2 : "); get(f2);
      put(" f2 : "); put(f2);
      Zero_Test(f2);
      Test64_Compare(f1,f2);
     -- Copy(f1,f2);
     -- put_line("After copy :");
     -- Test_Compare(f1,f2);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Comparison;

  procedure Test_Size is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    f : Floating_Number;
    ans : character;
    factor : integer32 := 0;
    rnd : boolean;

  begin
    put_line("Testing trunc/round/expand for multiprecision floats");
    loop
      put("Give a floating number : "); get(f);
      put("-> your floating : "); put(f); new_line;
      put("The size of the fraction : ");
      put(Size_Fraction(f),1); new_line;
      loop
        put("Give size modificator : "); get(factor);
        if factor <= 0 then
          put("Do you want to truncate or to round ? (t/r) "); get(ans);
          rnd := (ans = 'r');
        end if;
        if factor > 0 then
          Expand(f,natural32(factor));
          put("expanded : "); put(f);
          new_line;
        elsif factor < 0 then
          if rnd
           then Round(f,natural32(-factor)); put("rounded : ");
           else Trunc(f,natural32(-factor)); put("truncated : ");
          end if;
          put(f); -- put(mf);
          new_line;
        else
          if rnd
           then Round(f,natural32(factor)); put("rounded : ");
           else Trunc(f,natural32(factor)); put("truncated : ");
          end if;
          put(f);
          new_line;
          Expand(f,natural32(factor));
          put("expanded : "); put(f);
          new_line;
        end if;
        put("Do you want other size modificators ? (y/n) "); get(ans);
        exit when (ans /= 'y');
      end loop;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Size;

  procedure Test_Addition is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    ans : character;
    f1,f2,sum1,sum2 : Floating_Number;

  begin
    put_line("Testing the addition operations.");
    loop
      Read(f1,"f1");
     -- put("Give 1st number f1 : "); get(f1);
      put("-> f1 : "); put(f1); new_line;
      Read(f2,"f2");
     -- put("Give 2nd number f2 : "); get(f2);
      put("-> f2 : "); put(f2); new_line;
      sum1 := f1+f2;
      put("f1+f2 : "); put(sum1); new_line;
      sum2 := f2+f1;
      put("f2+f1 : "); put(sum2); new_line;
      if Equal(sum1,sum2)
       then put_line("Test on commutativity is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Addition;

  procedure Test64_Addition is

    use Multprec_Floating64_Numbers,Multprec_Floating64_Numbers_io;

    ans : character;
    f1,f2,sum1,sum2 : Floating_Number;

  begin
    put_line("Testing the addition operations.");
    loop
      Read(f1,"f1");
     -- put("Give 1st number f1 : "); get(f1);
      put("-> f1 : "); put(f1); new_line;
      Read(f2,"f2");
     -- put("Give 2nd number f2 : "); get(f2);
      put("-> f2 : "); put(f2); new_line;
      sum1 := f1+f2;
      put("f1+f2 : "); put(sum1); new_line;
      sum2 := f2+f1;
      put("f2+f1 : "); put(sum2); new_line;
      if Equal(sum1,sum2)
       then put_line("Test on commutativity is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Addition;

  procedure Test_Subtraction is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    ans : character;
    f1,f2,diff : Floating_Number;

  begin
    put_line("Testing the subtraction operations.");
    loop
      Read(f1,"f1");
     -- put("Give 1st number f1 : "); get(f1);
      put("-> f1 : "); put(f1); new_line;
      Read(f2,"f2");
     -- put("Give 2nd number f2 : "); get(f2);
      put("-> f2 : "); put(f2); new_line;
      diff := f1-f2;
      put("f1 - f2 : "); put(diff); new_line;
      Add(diff,f2);
      put("(f1-f2)+f2 : "); put(diff); new_line;
      if Equal(diff,f1)
       then put_line("Test of subtraction is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Subtraction;

  procedure Test64_Subtraction is

    use Multprec_Floating64_Numbers,Multprec_Floating64_Numbers_io;

    ans : character;
    f1,f2,diff : Floating_Number;

  begin
    put_line("Testing the subtraction operations.");
    loop
      Read(f1,"f1");
     -- put("Give 1st number f1 : "); get(f1);
      put("-> f1 : "); put(f1); new_line;
      Read(f2,"f2");
     -- put("Give 2nd number f2 : "); get(f2);
      put("-> f2 : "); put(f2); new_line;
      diff := f1-f2;
      put("f1 - f2 : "); put(diff); new_line;
      Add(diff,f2);
      put("(f1-f2)+f2 : "); put(diff); new_line;
      if Equal(diff,f1)
       then put_line("Test of subtraction is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Subtraction;

  procedure Test_Multiplication is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    ans : character;
    f1,f2,prod1,prod2 : Floating_Number;
    szf1,szf2,szprod1,szprod2 : natural32;

  begin
    put_line("Testing the multiplication operations.");
    loop
      Read(f1,"f1");
     -- put("Give 1st number : "); get(f1);
      szf1 := Size_Fraction(f1);
      put("Size of fraction of f1 : "); put(szf1,1); new_line;
      put("-> f1 : "); put(f1); new_line;
      Read(f2,"f2");
     -- put("Give 2nd number : "); get(f2);
      szf2 := Size_Fraction(f2);
      put("Size of fraction of f2 : "); put(szf2,1); new_line;
      put("-> f2 : "); put(f2); new_line;
      prod1 := f1*f2;
      szprod1 := Size_Fraction(prod1);
      put("Size of fraction of product : ");
      put(szprod1,1); new_line;
      put("Product f1*f2 : "); put(prod1); new_line;
      prod2 := f2*f1;
      szprod2 := Size_Fraction(prod2);
      put("Size of fraction of product : ");
      put(szprod2,1); new_line;
      put("Product f2*f1 : "); put(prod2); new_line;
      if Equal(prod1,prod2)
       then put_line("Test on commutativity is successful.");
       else put_line("Failure, product not commutative: bug!");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Multiplication;

  procedure Test_Exponentiation is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    ans : character;
    e1,e2 : Integer_Number;
    f,exp1,exp2,prod,expo : Floating_Number;

  begin
    put_line("Testing the exponentiation operations.");
    loop
      Read(f,"f");
     -- put("Give a number : "); get(f);
      put("-> your number f : "); put(f); new_line;
      put("Give 1st exponent : "); get(e1);
      put("-> your 1st exponent e1 : "); put(e1); new_line;
      exp1 := f**e1;
      put("f**e1 : "); put(exp1); new_line;
      put("Give 2nd exponent : "); get(e2);
      put("-> your 2nd exponent e2 : "); put(e2); new_line;
      exp2 := f**e2;
      put("f**e2 : "); put(exp2); new_line;
      prod := exp1*exp2;
      put("(f**e1)*(f**e2) : "); put(prod); new_line;
      expo := f**(e1+e2);
      put("f**(e1+e2)      : "); put(expo); new_line;
      if Equal(prod,expo)
       then put_line("Test of exponentiation is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Exponentiation;

  procedure Test_Division is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    ans : character;
    f1,f2,quot,prod,diff : Floating_Number;

  begin
    put_line("Testing the division operations.");
    loop
      Read(f1,"f1");
     -- put("Give 1st number f1 : "); get(f1);
      put("-> f1 : "); put(f1); new_line;
      Read(f2,"f2");
     -- put("Give 2nd number f2 : "); get(f2);
      put("-> f2 : "); put(f2); new_line;
      prod := f1*f2;
      put("f1*f2 : "); put(prod); new_line;
      quot := prod/f2;
      put("(f1*f2)/f2 : "); put(quot); new_line;
      if Equal(quot,f1)
       then put_line("Test of division is successful.");
       else put("Failure, bug detected?");
            put_line("  Difference up to working precision?");
            diff := quot - f1;
            put("(f1*f2)/f2 - f1 : "); put(diff); new_line;
      end if;
      quot := f1/f2;
      put("f1/f2 : "); put(quot); new_line; Clear(quot);
      Copy(f1,quot);
      Div(quot,f2);    put("f1/f2 : "); put(quot); new_line;
      prod := quot*f2; put("(f1/f2)*f2 : "); put(prod); new_line;
                       put(" f1        : "); put(f1); new_line;
      if Equal(prod,f1) then
        put_line("Test of division/remainder computation is successful.");
      else
        put("Failure, bug detected?");
        put_line("  Difference up to working precision?");
        if prod > f1
         then diff := prod - f1;
         else diff := f1 - prod;
        end if;
        put("(f1/f2)*f2 - f1 : "); put(diff); new_line;
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Division;

  function Min ( a,b : natural32 ) return natural32 is
  begin
    if a <= b
     then return a;
     else return b;
    end if;
  end Min;

  function Max ( a,b : natural32 ) return natural32 is
  begin
    if a >= b
     then return a;
     else return b;
    end if;
  end Max;

  function Acceptable_Tolerance
             ( sz1,sz2 : natural32 )
             return Multprec_Floating_Numbers.Floating_Number is

    use Multprec_Floating_Numbers;

    res : Floating_Number;
    the_radix : constant natural32 := Multprec_Natural_Coefficients.Radix;
    the_expo : constant natural32 := Multprec_Natural_Coefficients.Exponent;
    sz : constant natural32 := Min(sz1,sz2) + 1;
    exp : constant integer32 := -integer32(sz*the_expo+2);

  begin
    res := Create(integer32(the_radix),exp);
    return res;
  end Acceptable_Tolerance;

  procedure Report_Bug ( n1,n2 : Multprec_Floating_Numbers.Floating_Number ) is

    use Multprec_Floating_Numbers_io;

  begin
    new_line;
    put("  n1 : "); put(n1); new_line;
    put("  n2 : "); put(n2); new_line;
  end Report_Bug;

  procedure Compare_with_Tolerance
              ( n1,n2,resi : in Multprec_Floating_Numbers.Floating_Number; 
                acctol : in Multprec_Floating_Numbers.Floating_Number ) is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    absresi : Floating_Number := AbsVal(resi);

  begin
    if absresi < acctol then
      put(" < "); put(acctol,3); put_line(" okay");
    else
      put(" > "); put(acctol,3); put_line(" Bug!");
      Report_Bug(n1,n2);
    end if;
    Clear(absresi);
  end Compare_with_Tolerance;

  procedure Random_Addition_and_Subtraction
              ( sz1,sz2 : in natural32; low,upp : in integer32;
                acctol : in Multprec_Floating_Numbers.Floating_Number ) is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    n1,n2,sum1,sum2,tmp : Floating_Number;
    szsum1,szsum2 : natural32;

  begin
    n1 := Random(sz1,low,upp);
    n2 := Random(sz2,low,upp);
    sum1 := n1+n2;
    szsum1 := Size_Fraction(sum1);
    sum2 := sum1-n2;
    szsum2 := Size_Fraction(sum2);
    if szsum1 <= max(sz1,sz2) and szsum2 <= max(sz1,sz2)
     then put("size okay  ");
     else put("size bug!  "); Report_Bug(n1,n2);
    end if;
    if Equal(sum2,n1) then
      put("n1+n2-n2 okay");
    else
      put("n1+n2-n2 Bug?");
      put("  diff : "); tmp := sum2-n1; put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Add(sum2,n2);
    szsum2 := Size_Fraction(sum2);
    if szsum2 <= max(sz1,sz2)
     then put("  size okay");
     else put("  size Bug!");
    end if;
    if Equal(sum2,sum1) then
      put("  Add okay");
    else
      put("  Add Bug?");
      put("  diff : "); tmp := sum2-sum1; put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Sub(sum2,n1);
    szsum2 := Size_Fraction(sum2);
    if szsum2 <= max(sz1,sz2)
     then put("  size okay");
     else put("  size Bug!");
    end if;
    if Equal(sum2,n2) then
      put("  Sub okay"); new_line;
    else
      put("  Sub Bug?");
      put("  diff : "); tmp := sum2-n2; put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Clear(n1); Clear(n2);
    Clear(sum1); Clear(sum2);
  exception
    when others => put_line("input caused exception:"); 
                   Report_Bug(n1,n2); raise;
  end Random_Addition_and_Subtraction;

  procedure Additions_and_Subtractions_on_Randoms is

    use Multprec_Floating_Numbers;

    nb,sz1,sz2 : natural32 := 0;
    low,upp : integer32 := 0;
    timer : Timing_Widget;
    acctol : Floating_Number;

  begin
    put("Give the number of tests : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    put("Give lower bound on exponent : "); get(low);
    put("Give upper bound on exponent : "); get(upp);
    acctol := Acceptable_Tolerance(sz1,sz2);
    tstart(timer);
    for i in 1..nb loop
      Random_Addition_and_Subtraction(sz1,sz2,low,upp,acctol);
    end loop;
    tstop(timer);
    Clear(acctol);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Additions_and_Subtractions_on_Randoms;

  procedure Random_Multiplication_and_Division
               ( sz1,sz2 : in natural32; low,upp : in integer32;
                 acctol : in Multprec_Floating_Numbers.Floating_Number ) is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;

    n1,n2,prod,quot,tmp : Floating_Number;
    szprod,szquot : natural32;

  begin
    n1 := Random(sz1,low,upp);
    n2 := Random(sz2,low,upp);
    prod := n1*n2;
    szprod := Size_Fraction(prod);
    quot := prod/n2;
    szquot := Size_Fraction(quot);
    if szprod <= max(sz1,sz2) and szquot <= max(sz1,sz2)
     then put("size okay  ");
     else put("size bug!  ");
    end if;
    if Equal(quot,n1) then
      put("n1*n2/n2 okay");
    else
      put("n1*n2/n2 Bug?");
      put("  Diff : "); tmp := quot-n1; put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Mul(quot,n2);
    szquot := Size_Fraction(quot);
    if szquot <= max(sz1,sz2)
     then put("  size okay");
     else put("  size Bug!");
    end if;
    if Equal(prod,quot) then
      put("  Mul okay");
    else
      put("  Mul Bug?");
      put("  Diff : "); tmp := quot-prod; put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Div(prod,n1);
    szprod := Size_Fraction(prod);
    if szprod <= max(sz1,sz2)
     then put("  size okay");
     else put("  size Bug!");
    end if;
    if Equal(prod,n2) then
      put("  Div okay"); new_line;
    else
      put("  Div Bug?");
      put("  Diff : "); tmp := prod-n2; put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Clear(n1); Clear(n2);
    Clear(prod); Clear(quot);
  exception
    when others => put_line("input caused exception :");
                   Report_Bug(n1,n2); raise;
  end Random_Multiplication_and_Division;

  procedure Multiplications_and_Divisions_on_Randoms is

    use Multprec_Floating_Numbers;

    nb,sz1,sz2 : natural32 := 0;
    low,upp : integer32 := 0;
    acctol : Floating_Number;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    put("Give lower bound on exponent : "); get(low);
    put("Give upper bound on exponent : "); get(upp);
    acctol := Acceptable_Tolerance(sz1,sz2);
    tstart(timer);
    for i in 1..nb loop
      Random_Multiplication_and_Division(sz1,sz2,low,upp,acctol);
    end loop;
    tstop(timer);
    Clear(acctol);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Multiplications_and_Divisions_on_Randoms;

  procedure Multiplication_Timer is

    use Multprec_Floating_Numbers;

    nb,sz1,sz2 : natural32 := 0;
    low,upp : integer32 := 0;
    timer : Timing_Widget;
    n1,n2,prod : Floating_Number;

  begin
    put("Give the number of cycles : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    put("Give lower bound on the exponent : "); get(low);
    put("Give upper bound on the exponent : "); get(upp);
    tstart(timer);
    for i in 1..nb loop
      n1 := Random(sz1,low,upp);
      n2 := Random(sz2,low,upp);
      prod := n1*n2;
      Clear(n1); Clear(n2); Clear(prod);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Multiplication_Timer;

  procedure Multiplication64_Timer is

    use Multprec_Floating64_Numbers;

    nb,sz1,sz2 : natural32 := 0;
    low,upp : integer64 := 0;
    timer : Timing_Widget;
    n1,n2,prod : Floating_Number;

  begin
    put("Give the number of cycles : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    put("Give lower bound on the exponent : "); get(low);
    put("Give upper bound on the exponent : "); get(upp);
    tstart(timer);
    for i in 1..nb loop
      n1 := Random(sz1,low,upp);
      n2 := Random(sz2,low,upp);
      prod := n1*n2;
      Clear(n1); Clear(n2); Clear(prod);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Multiplication64_Timer;

  procedure Division_Timer is

    use Multprec_Floating_Numbers;

    nb,sz1,sz2 : natural32 := 0;
    low,upp : integer32 := 0;
    timer : Timing_Widget;
    n1,n2,quot : Floating_Number;

  begin
    put("Give the number of cycles : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    put("Give lower bound on the exponent : "); get(low);
    put("Give upper bound on the exponent : "); get(upp);
    tstart(timer);
    for i in 1..nb loop
      n1 := Random(sz1,low,upp);
      n2 := Random(sz2,low,upp);
      quot := n1/n2;
      Clear(quot);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Division_Timer;

  procedure Division64_Timer is

    use Multprec_Floating64_Numbers;

    nb,sz1,sz2 : natural32 := 0;
    low,upp : integer64 := 0;
    timer : Timing_Widget;
    n1,n2,quot : Floating_Number;

  begin
    put("Give the number of cycles : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    put("Give lower bound on the exponent : "); get(low);
    put("Give upper bound on the exponent : "); get(upp);
    tstart(timer);
    for i in 1..nb loop
      n1 := Random(sz1,low,upp);
      n2 := Random(sz2,low,upp);
      quot := n1/n2;
      Clear(quot);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Division64_Timer;

  procedure Nearest_Integer is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;
    f : Floating_Number;
    i : Integer_Number;

  begin
    put("Give your float : "); get(f);
    put("-> your float is "); put(f); new_line;
    put("Fraction : "); put(Fraction(f)); new_line;
    put("Exponent : "); put(Exponent(f)); new_line;
    i := Truncate_to_Nearest_Integer(f);
    put("Nearest integer : "); put(i); new_line;
  end Nearest_Integer;

  procedure Nearest_Integer64 is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;
    f : Floating_Number;
    i : Integer_Number;

  begin
    put("Give your float : "); get(f);
    put("-> your float is "); put(f); new_line;
    put("Fraction : "); put(Fraction(f)); new_line;
    put("Exponent : "); put(Exponent(f)); new_line;
    i := Truncate_to_Nearest_Integer(f);
    put("Nearest integer : "); put(i); new_line;
  end Nearest_Integer64;

  procedure Main is

    ans,lng : character := ' ';

  begin
    new_line;
    put_line("Interactive testing of multiprecision floating numbers.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit program      1. Input/Output     2. Creation      ");
      put_line("  3. Comparison/Copy   4. Addition         5. Subtraction   ");
      put_line("  6. Multiplication    7. Exponentiation   8. Division      ");
      put_line("  9. Truncate/Round/Expand                                  ");
      put_line("  A. Addition/subtraction on randomly generated numbers.    ");
      put_line("  B. Multiplication/division/remainder on random numbers.   ");
      put_line("  C. Timer on multiplication cycle on random numbers.       ");
      put_line("  D. Timer on division cycle on random numbers.             ");
      put_line("  E. Nearest integer of a multiprecision floating number.   ");
      put("Type in your choice (0,1,2,3,4,5,6,7,8,9,A,B,C,D, or E) : ");
      Ask_Alternative(ans,"0123456789ABCDE");
      exit when (ans = '0');
      new_line;
      if ans = '1' or ans = '2' or ans = '3' or ans = '4' or ans = '5'
        or ans = 'C' or ans = 'D' or ans = 'E'
       then put("Use 64-bit arithmetic ? (y/n) "); Ask_Yes_or_No(lng);
            new_line;
      end if;
      if lng = 'y' then
        case ans is
          when '1' => Test64_io;
          when '2' => Test64_Creation;
          when '3' => Test64_Comparison;
          when '4' => Test64_Addition;
          when '5' => Test64_Subtraction;
          when 'C' => Multiplication64_Timer;
          when 'D' => Division64_Timer;
          when 'E' => Nearest_Integer64;
          when others => null;
        end case;
        lng := 'n'; -- for the next test ...
      else
        case ans is
          when '1' => Test_io;
          when '2' => Test_Creation;
          when '3' => Test_Comparison;
          when '4' => Test_Addition;
          when '5' => Test_Subtraction;
          when '6' => Test_Multiplication;
          when '7' => Test_Exponentiation;
          when '8' => Test_Division;
          when '9' => Test_Size;
          when 'A' => Additions_and_Subtractions_on_Randoms;
          when 'B' => Multiplications_and_Divisions_on_Randoms;
          when 'C' => Multiplication_Timer;
          when 'D' => Division_Timer;
          when 'E' => Nearest_Integer;
          when others => null;
        end case;
      end if;
    end loop;
  end Main;

end Test_Floating_Numbers;
