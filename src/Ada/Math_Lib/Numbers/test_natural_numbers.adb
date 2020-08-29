with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Multprec_Natural_Coefficients;
with Multprec_Natural_Numbers_io;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;
with Multprec_Natural64_Coefficients;
with Multprec_Natural64_Numbers;
with Multprec_Natural64_Numbers_io;

package body Test_Natural_Numbers is

  procedure Test_Creation is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;

    n : natural32 := 0;
    sc : natural32;
    nn : Natural_Number;
    ans : character;

  begin
    put_line("Testing the creation of a natural number.");
    loop
      put("Give a standard natural number : "); get(n);
      sc := Multprec_Natural_Coefficients.Size_of_Coefficient(n);
      put("-> number of digits : "); put(sc,1);
      if sc > Multprec_Natural_Numbers.Exponent
       then put_line(" not a coefficient");
       else put_line(" is a coefficient");
      end if;
      nn := Create(n);
      put("-> as natural number : "); put(nn); new_line;
      Clear(nn);
      put("Do you want more tests ? (y/n) "); Ask_Yes_or_No(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Creation;

  procedure Test64_Creation is

    n : natural64 := 0;
    sc : natural32;
    ans : character;
    nn : Multprec_Natural64_Numbers.Natural_Number;

  begin
    put_line("Testing the creation of a 64-bit natural number.");
    loop
      put("Give a standard natural number : "); get(n);
      sc := Multprec_Natural64_Coefficients.Size_of_Coefficient(n);
      put("-> number of digits : "); put(sc,1);
      if sc > Multprec_Natural64_Numbers.Exponent
       then put_line(" not a coefficient");
       else put_line(" is a coefficient");
      end if;
      put("-> your number : "); put(n); new_line;
      nn := Multprec_Natural64_Numbers.Create(n);
      put("as a multiprecision natural : ");
      Multprec_Natural64_Numbers_io.put(nn); new_line;
      put("Do you want more tests ? (y/n) "); Ask_Yes_or_No(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Creation;

  procedure Test_io is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;
    use Multprec_Natural_Coefficients;

    ans : character;
    n : Natural_Number;

  begin
    put_line("Testing the input/output operations.");
    loop
      put("Give a number : "); get(n);
      put("-> your number : "); put(n); new_line;
      declare
        c : constant Array_of_Naturals := Coefficients(n);
        m : constant natural32
          := (Size(n)+1)*Multprec_Natural_Coefficients.Exponent;
        sc : string(1..natural(m));
        d : constant natural32 := Decimal_Places(n);
        sn : string(1..natural(d));
        n2 : Natural_Number;
      begin
        put("#decimal places : "); put(d,1);
        put("  Number of characters : "); put(m,1); new_line;
        put("The coefficients of the number : "); put(c); new_line;
        put(sc,c);
        put("String representation of coeff : "); put(sc); new_line;
        put(sn,n);
        put("String representation of number : "); put(sn); new_line;
        get(sn,n2);
        put("Number from parsing the string  : "); put(n2); new_line;
      end;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_io;

  procedure Test64_io is

    use Multprec_Natural64_Numbers,Multprec_Natural64_Numbers_io;
    use Multprec_Natural64_Coefficients;

    ans : character;
    n : Natural_Number;

  begin
    put_line("Testing the input/output operations.");
    loop
      put("Give a number : "); get(n);
      put("-> your number : "); put(n); new_line;
      declare
        c : constant Array_of_Naturals := Coefficients(n);
        m : constant natural32
          := (Size(n)+1)*Multprec_Natural64_Coefficients.Exponent;
        sc : string(1..natural(m));
        d : constant natural32 := Decimal_Places(n);
        sn : string(1..natural(d));
        n2 : Natural_Number;
      begin
        put("#decimal places : "); put(d,1);
        put("  Number of characters : "); put(m,1); new_line;
        put("The coefficients of the number : "); put(c); new_line;
        put(sc,c);
        put("String representation of coeff : "); put(sc); new_line;
        put(sn,n);
        put("String representation of number : "); put(sn); new_line;
        get(sn,n2);
        put("Number from parsing the string  : "); put(n2); new_line;
      end;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_io;

 -- procedure Test_Compare ( n1 : in Multprec_Natural_Numbers.Natural_Number;
 --                          n2 : in natural ) is
 --
 --   use Multprec_Natural_Numbers;
 -- 
 -- begin
 --   if Equal(n1,n2)
 --    then put_line("The numbers are equal.");
 --    else put_line("The numbers are different.");
 --   end if;
 --   if n1 < n2
 --    then put_line("First is less than second.");
 --    else put_line("First not less than second.");
 --   end if;
 --   if n1 > n2
 --    then put_line("First is greater than second.");
 --    else put_line("First is not greater than second.");
 --   end if;
 -- end Test_Compare;

 -- procedure Test_Compare
 --              ( n1 : in natural;
 --                n2 : in Multprec_Natural_Numbers.Natural_Number ) is
 --
 --   use Multprec_Natural_Numbers;
 --
 -- begin
 --   if Equal(n2,n1)
 --    then put_line("Both numbers are equal.");
 --    else put_line("The numbers are different.");
 --   end if;
 --   if n1 < n2
 --    then put_line("First is less than second.");
 --    else put_line("First not less than second.");
 --   end if;
 --   if n1 > n2
 --    then put_line("First is greater than second.");
 --    else put_line("First is not greater than second.");
 --   end if;
 -- end Test_Compare;

  procedure Test_Compare
              ( n1,n2 : in Multprec_Natural_Numbers.Natural_Number ) is

    use Multprec_Natural_Numbers;

  begin
    if Equal(n1,n2)
     then put_line("Both numbers are equal.");
     else put_line("The numbers are different.");
    end if;
    if n1 < n2
     then put_line("First is less than second.");
     else put_line("First not less than second.");
    end if;
    if n1 > n2
     then put_line("First is greater than second.");
     else put_line("First is not greater than second.");
    end if;
  end Test_Compare;

  procedure Test_Comparison is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;

    ans : character;
    n1,n2 : Natural_Number;
   -- n1 : natural;

  begin 
    put_line("Testing the comparison operations.");
    loop
      put("Give 1st number n1 : "); get(n1);
      put("-> n1 : "); put(n1); new_line;
      put("Give 2nd number n2 : "); get(n2);
      put("-> n2 : "); put(n2); new_line;
      Test_Compare(n1,n2);
      Copy(n1,n2);
      put_line("Tests after copying : ");
      Test_Compare(n1,n2);
      Div(n1,10);
      put_line("After dividing n1 by 10 :"); 
      put(" n1 : "); put(n1); new_line;
      put(" n2 : "); put(n2); new_line;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Comparison;

  procedure Test_Addition is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;
 
    ans : character;
    n1,n2,sum1,sum2 : Natural_Number;

  begin
    put_line("Testing the addition operations.");
    loop
      put("Give 1st number : "); get(n1);
      put("-> your 1st number n1 : "); put(n1); new_line;
      put("Give 2nd number : "); get(n2);
      put("-> your 2nd number n2 : "); put(n2); new_line;
      sum1 := n1+n2;
      put("n1+n2 : "); put(sum1); new_line;
      sum2 := n2+n1;
      put("n2+n1 : "); put(sum2); new_line;
      if Equal(sum1,sum2)
       then put_line("Test on commutativity is successful.");
       else put_line("Failure, sum not commutative: bug!");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Addition;

  procedure Test64_Addition is

    use Multprec_Natural64_Numbers,Multprec_Natural64_Numbers_io;
 
    ans : character;
    n1,n2,sum1,sum2 : Natural_Number;

  begin
    put_line("Testing the addition operations.");
    loop
      put("Give 1st number : "); get(n1);
      put("-> your 1st number n1 : "); put(n1); new_line;
      put("Give 2nd number : "); get(n2);
      put("-> your 2nd number n2 : "); put(n2); new_line;
      sum1 := n1+n2;
      put("n1+n2 : "); put(sum1); new_line;
      sum2 := n2+n1;
      put("n2+n1 : "); put(sum2); new_line;
      if Equal(sum1,sum2)
       then put_line("Test on commutativity is successful.");
       else put_line("Failure, sum not commutative: bug!");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Addition;

 -- function Mult_by_Add ( n1 : Multprec_Natural_Numbers.Natural_Number;
 --                        n2 : natural ) 
 --                      return Multprec_Natural_Numbers.Natural_Number is
 --
  -- DESCRIPTION :
  --   Does the multiplication by adding up n1 to itself as many times
  --   as the number n2.  Only to be used as test of course.
 --
 --   use Multprec_Natural_Numbers;
 --
 --   res : Natural_Number;
 --
 -- begin
 --   if n2 = 0 then
 --     return res;
 --   else
 --     Copy(n1,res);
 --     for i in 1..n2-1 loop
 --       Add(res,n1);
 --     end loop;
 --     return res;
 --   end if;
 -- end Mult_by_Add;

  function Mult_by_Add ( n1,n2 : Multprec_Natural_Numbers.Natural_Number )
                       return Multprec_Natural_Numbers.Natural_Number is

    use Multprec_Natural_Numbers;
    res,cnt : Natural_Number;

  begin
    if Equal(n2,0) then
      return res;
    else
      Copy(n1,res);
      cnt := Create(natural32(1));
      while not Equal(cnt,n2) loop
        Add(res,n1);
        Add(cnt,1);
      end loop;
      Clear(cnt);
      return res;
    end if;
  end Mult_by_Add;

  procedure Test_Multiplication is

  -- NOTE : to test n1*n2 with n2 : natural, change the declaration of n2.

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;
 
    ans : character;
    n1,n2,prod1,prod2,prod3 : Natural_Number;
   -- n2 : natural := 0;

  begin
    put_line("Testing the multiplication operations.");
    loop
      put("Give 1st number : "); get(n1);
      put("-> your 1st number : "); put(n1); new_line;
      put("Give 2nd number : "); get(n2);
      put("-> your 2nd number : "); put(n2); new_line;
      prod1 := n1*n2;
      put("Product n1*n2 : "); put(prod1); new_line;
      prod2 := n2*n1;
      put("Product n2*n1 : "); put(prod2); new_line;
      if Equal(prod1,prod2)
       then put_line("Test on commutativity is successful.");
       else put_line("Failure, product not commutative: bug!");
      end if;
      put("Do you want multiplication by addition ? (y/n) "); get(ans);
      if ans = 'y' then
        put_line("Testing the multiplication by addition.  Be patient...");
        prod3 := Mult_by_Add(n1,n2);
        put("After adding "); put(n2); put(" times : "); put(prod3);
        new_line;
        if Equal(prod1,prod3)
         then put_line("Test of multiplication is successful.");
         else put_line("Failure, bug detected.");
        end if;
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Multiplication;

  procedure Test_Exponentiation is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;

    ans : character;
    n,e1,e2,exp1,exp2,prod,expo : Natural_Number;

  begin
    put_line("Testing the exponentiation operations.");
    loop
      put("Give a number : "); get(n);
      put("-> your number n : "); put(n); new_line;
      put("Give 1st exponent : "); get(e1);
      put("-> your 1st exponent e1 : "); put(e1); new_line;
      exp1 := n**e1;
      put("n**e1 : "); put(exp1); new_line;
      put("Give 2nd exponent : "); get(e2);
      put("-> your 2nd exponent e2 : "); put(e2); new_line;
      exp2 := n**e2;
      put("n**e2 : "); put(exp2); new_line;
      prod := exp1*exp2;
      put("(n**e1)*(n**e2) : "); put(prod); new_line;
      expo := n**(e1+e2);
      put("n**(e1+e2)      : "); put(expo); new_line;
      if Equal(prod,expo)
       then put_line("Test of exponentiation is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Exponentiation;

  procedure Test_Subtraction is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;

    ans : character;
    n1,n2,diff : Natural_Number;

  begin
    put_line("Testing the subtraction operations.");
    loop
      put("Give 1st number : "); get(n1);
      put("-> your 1st number n1 : "); put(n1); new_line;
      put("Give 2nd number : "); get(n2);
      put("-> your 2nd number n2 : "); put(n2); new_line;
      diff := n1-n2;
      put("n1 - n2 : "); put(diff); new_line;
      Add(diff,n2);
      put("(n1-n2)+n2 : "); put(diff); new_line;
      if Equal(diff,n1)
       then put_line("Test of subtraction is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Subtraction;

  procedure Divide10 ( n : in Multprec_Natural_Numbers.Natural_Number ) is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;

    quot,prod : Natural_Number;
    rest : natural32;

  begin
    put("n : "); put(n); new_line;
    for i in 1..10 loop
      rest := Rmd(n,natural32(i));
      quot := n/natural32(i);
      if rest = 0
       then put("Divisible by "); put(natural32(i),1);
       else put("Not divisible by "); put(natural32(i),1);
      end if;
      put("  rest : "); put(rest,1); new_line;
      put("quotient : "); put(quot); new_line;
      prod := quot*natural32(i) + rest;
      if Equal(prod,n)
       then put_line("Test on Remainder/Division is successful.");
       else put_line("Failure, bug detected.");
      end if;
    end loop;
  end Divide10;

  procedure Test_Short_Division is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;
    n1,quot,prod : Natural_Number;
    n2,rest : natural32 := 0;

  begin
    put("Give 1st number : "); get(n1);
    put("-> your 1st number n1 : "); put(n1); new_line;
    put("Give 2nd number : "); get(n2);
    put("-> your 2nd number n2 : "); put(n2);
    if n2 <= Multprec_Natural_Numbers.Base then
      put(" <= "); put(Multprec_Natural_Numbers.Base); new_line;
    else
      put(" > "); put(Multprec_Natural_Numbers.Base);
      put_line(" please try again!"); return;
    end if;
    prod := n1*n2;
    put("n1*n2 : "); put(prod); new_line;
    quot := prod/n2; rest := Rmd(prod,n2);
    put("(n1*n2)/n2 : "); put(quot); new_line;
    put("Remainder : "); put(rest); new_line;
    if Equal(quot,n1) and (rest = 0)
     then put_line("Test of division is successful.");
     else put_line("Failure, bug detected.");
    end if;
    Div(n1,n2,quot,rest);
    put("n1/n2 : "); put(quot); new_line;
    put("rest : "); put(rest); new_line;
    prod := quot*n2 + rest;
    if Equal(prod,n1)
     then put_line("Test of division/remainder computation is successful.");
     else put_line("Failure, bug detected.");
    end if;
    if n2 <= 10
     then Divide10(n1);
    end if;
  end Test_Short_Division;

  procedure Test_Long_Division is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;
    n1,n2,quot,prod,rest : Natural_Number;

  begin
    put("Give 1st number : "); get(n1);
    put("-> your 1st number n1 : "); put(n1); new_line;
    put("Give 2nd number : "); get(n2);
    put("-> your 2nd number n2 : "); put(n2); new_line;
    prod := n1*n2;
    put("n1*n2 : "); put(prod); new_line;
   -- quot := prod/n2; rest := Rmd(prod,n2);
    Div(prod,n2,quot,rest);
    put("(n1*n2)/n2 : "); put(quot); new_line;
    put("Remainder : "); put(rest); new_line;
    if Equal(quot,n1) and Equal(rest,0)
     then put_line("Test of division is successful.");
     else put_line("Failure, bug detected.");
    end if;
    Div(n1,n2,quot,rest);
    put("n1/n2 : "); put(quot); new_line;
    put("rest : "); put(rest); new_line;
    prod := quot*n2 + rest;
    if Equal(prod,n1)
     then put_line("Test of division/remainder computation is successful.");
     else put_line("Failure, bug detected.");
          put("q*n2 + r : "); put(prod); new_line;
          put("      n1 : "); put(n1); new_line;
    end if;
  end Test_Long_Division;

  procedure Test_Division ( short : in boolean ) is

    ans : character;

  begin
    put_line("Testing the division operations.");
    loop
      if short
       then Test_Short_Division;
       else Test_Long_Division;
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Division;

  procedure Test64_Division is

    use Multprec_Natural64_Numbers,Multprec_Natural64_Numbers_io;

    ans : character;
    n1,n2,quot,prod,rest : Natural_Number;
   -- n2,rest : natural;

  begin
    put_line("Testing the division operations.");
    loop
      put("Give 1st number : "); get(n1);
      put("-> your 1st number n1 : "); put(n1); new_line;
      put("Give 2nd number : "); get(n2);
      put("-> your 2nd number n2 : "); put(n2); new_line;
      prod := n1*n2;
      put("n1*n2 : "); put(prod); new_line;
      quot := prod/n2; rest := Rmd(prod,n2);
      put("(n1*n2)/n2 : "); put(quot); new_line;
      put("Remainder : "); put(rest); new_line;
      if Equal(quot,n1) and Equal(rest,0)
       then put_line("Test of division is successful.");
       else put_line("Failure, bug detected.");
      end if;
      Div(n1,n2,quot,rest);
      put("n1/n2 : "); put(quot); new_line;
      put("rest : "); put(rest); new_line;
      prod := quot*n2 + rest;
      if Equal(prod,n1)
       then put_line("Test of division/remainder computation is successful.");
       else put_line("Failure, bug detected.");
      end if;
     -- if n2 <= 10
     --  then Divide10(n1);
     -- end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Division;

  procedure Random_Addition_and_Subtraction ( sz1,sz2 : in natural32 ) is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;

    n1,n2,sum1,sum2 : Natural_Number;

    procedure Report_Bug is
    begin
      new_line;
      put("  n1 : "); put(n1); new_line;
      put("  n2 : "); put(n2); new_line;
    end Report_Bug;

  begin
    n1 := Random(sz1);
    n2 := Random(sz2);
    sum1 := n1+n2;
    sum2 := sum1-n2;
    if Equal(sum2,n1)
     then put("n1+n2-n2 okay");
     else put("n1+n2-n2 Bug!"); Report_Bug;
    end if;
    Add(sum2,n2);
    if Equal(sum2,sum1)
     then put("  Add okay");
     else put("  Add Bug!"); Report_Bug;
    end if;
    Sub(sum2,n1);
    if Equal(sum2,n2)
     then put("  Sub okay"); new_line;
     else put("  Sub Bug!"); Report_Bug;
    end if;
    Clear(n1); Clear(n2);
    Clear(sum1); Clear(sum2);
  end Random_Addition_and_Subtraction;

  procedure Additions_and_Subtractions_on_Randoms is

    nb,sz1,sz2 : natural32 := 0;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    tstart(timer);
    for i in 1..nb loop
      Random_Addition_and_Subtraction(sz1,sz2);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Additions_and_Subtractions_on_Randoms;

  procedure Random_Multiplication_and_Division ( sz1,sz2 : in natural32 ) is

    use Multprec_Natural_Numbers,Multprec_Natural_Numbers_io;

    n1,n2,prod1,prod2,quot1,quot2,quot3,rest1,rest2 : Natural_Number;

    procedure Report_Bug is
    begin
      new_line;
      put("  n1 : "); put(n1); new_line;
      put("  n2 : "); put(n2); new_line;
    end Report_Bug;

  begin
    n1 := Random(sz1);
    n2 := Random(sz2);
    prod1 := n1*n2;
    quot1 := prod1/n2;
    if Equal(quot1,n1)
     then put("n1*n2/n2 okay");
     else put("n1*n2/n2 Bug!"); Report_Bug;
    end if;
    Mul(quot1,n2);
    if Equal(prod1,quot1)
     then put("  Mul okay");
     else put("  Mul Bug!"); Report_Bug;
    end if;
    Div(prod1,n1);
    if Equal(prod1,n2)
     then put("  Div okay");
     else put("  Div Bug!"); Report_Bug;
    end if;
    rest1 := Rmd(n1,n2);
    quot2 := n1/n2;
    prod2 := quot2*n2;
    Add(prod2,rest1);
    if Equal(prod2,n1)
     then put("  Rmd okay");
     else put("  Rmd Bug!"); Report_Bug;
    end if;
    Div(n1,n2,quot3,rest2);
    Mul(quot3,n2);
    Add(quot3,rest2);
    if Equal(quot3,n1)
     then put("  Div/Rmd okay"); new_line;
     else put("  Div/Rmd Bug!"); Report_Bug;
    end if;
    Clear(n1);    Clear(n2);
    Clear(prod1); Clear(quot1);
    Clear(prod2); Clear(quot2);
    Clear(quot3); Clear(rest1); Clear(rest2);
  end Random_Multiplication_and_Division;

  procedure Random64_Multiplication_and_Division ( sz1,sz2 : in natural32 ) is

    use Multprec_Natural64_Numbers,Multprec_Natural64_Numbers_io;

    n1,n2,prod1,prod2,quot1,quot2,quot3,rest1,rest2 : Natural_Number;

    procedure Report_Bug is
    begin
      new_line;
      put("  n1 : "); put(n1); new_line;
      put("  n2 : "); put(n2); new_line;
    end Report_Bug;

  begin
    n1 := Random(sz1);
    n2 := Random(sz2);
    prod1 := n1*n2;
    quot1 := prod1/n2;
    if Equal(quot1,n1)
     then put("n1*n2/n2 okay");
     else put("n1*n2/n2 Bug!"); Report_Bug;
    end if;
    Mul(quot1,n2);
    if Equal(prod1,quot1)
     then put("  Mul okay");
     else put("  Mul Bug!"); Report_Bug;
    end if;
    Div(prod1,n1);
    if Equal(prod1,n2)
     then put("  Div okay");
     else put("  Div Bug!"); Report_Bug;
    end if;
    rest1 := Rmd(n1,n2);
    quot2 := n1/n2;
    prod2 := quot2*n2;
    Add(prod2,rest1);
    if Equal(prod2,n1)
     then put("  Rmd okay");
     else put("  Rmd Bug!"); Report_Bug;
    end if;
    Div(n1,n2,quot3,rest2);
    Mul(quot3,n2);
    Add(quot3,rest2);
    if Equal(quot3,n1)
     then put("  Div/Rmd okay"); new_line;
     else put("  Div/Rmd Bug!"); Report_Bug;
    end if;
    Clear(n1);    Clear(n2);
    Clear(prod1); Clear(quot1);
    Clear(prod2); Clear(quot2);
    Clear(quot3); Clear(rest1); Clear(rest2);
  end Random64_Multiplication_and_Division;

  procedure Multiplications_and_Divisions_on_Randoms is

    nb,sz1,sz2 : natural32 := 0;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    tstart(timer);
    for i in 1..nb loop
      Random_Multiplication_and_Division(sz1,sz2);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Multiplications_and_Divisions_on_Randoms;

  procedure Multiplications64_and_Divisions_on_Randoms is

    nb,sz1,sz2 : natural32 := 0;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    tstart(timer);
    for i in 1..nb loop
      Random64_Multiplication_and_Division(sz1,sz2);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Multiplications64_and_Divisions_on_Randoms;

  procedure Multiplication_Timer is

    use Multprec_Natural_Numbers;

    nb,sz1,sz2 : natural32 := 0;
    timer : Timing_Widget;
    n1,n2,prod : Natural_Number;

  begin
    put("Give the number of cycles : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    tstart(timer);
    for i in 1..nb loop
      n1 := Random(sz1);
      n2 := Random(sz2);
      prod := n1*n2;
      Clear(n1); Clear(n2); Clear(prod);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Multiplication_Timer;

  procedure Multiplication64_Timer is

    use Multprec_Natural64_Numbers;

    nb,sz1,sz2 : natural32 := 0;
    timer : Timing_Widget;
    n1,n2,prod : Natural_Number;

  begin
    put("Give the number of cycles : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    tstart(timer);
    for i in 1..nb loop
      n1 := Random(sz1);
      n2 := Random(sz2);
      prod := n1*n2;
      Clear(n1); Clear(n2); Clear(prod);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Multiplication64_Timer;

  procedure Division_Remainder_Timer is

    use Multprec_Natural_Numbers;

    nb,sz1,sz2 : natural32 := 0;
    timer : Timing_Widget;
    n1,n2,quot,rest : Natural_Number;

  begin
    put("Give the number of cycles : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    tstart(timer);
    for i in 1..nb loop
      n1 := Random(sz1);
      n2 := Random(sz2);
      Div(n1,n2,quot,rest);
      Clear(n1); Clear(n2);
      Clear(quot); Clear(rest);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Division_Remainder_Timer;

  procedure Division64_Remainder_Timer is

    use Multprec_Natural64_Numbers;

    nb,sz1,sz2 : natural32 := 0;
    timer : Timing_Widget;
    n1,n2,quot,rest : Natural_Number;

  begin
    put("Give the number of cycles : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    tstart(timer);
    for i in 1..nb loop
      n1 := Random(sz1);
      n2 := Random(sz2);
      Div(n1,n2,quot,rest);
      Clear(n1); Clear(n2);
      Clear(quot); Clear(rest);
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Division64_Remainder_Timer;

  procedure Main is

    ans,lng,shrt : character;

  begin
    new_line;
    put_line("Interactive testing of multi-precision natural numbers.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit program      1. Input/Output     2. Creation      ");
      put_line("  3. Comparison/Copy   4. Addition         5. Subtraction   ");
      put_line("  6. Multiplication    7. Exponentiation   8. Division      ");
      put_line("  9. Addition/subtraction on randomly generated numbers.    ");
      put_line("  A. Multiplication/division/remainder on random numbers.   ");
      put_line("  B. Timer on cycle of multiplication on random numbers.    ");
      put_line("  C. Timer on cycle of division/remainder on random numbers.");
      put("Type in your choice (0,1,2,3,4,5,6,7,8,9,A,B or C) : ");
      Ask_Alternative(ans,"0123456789ABC");
      exit when (ans = '0');
      new_line;
      if ans = '1' or ans = '2' or ans = '4' or ans = '8' 
         or ans = 'A' or ans = 'B' or ans = 'C' then
        put("Use 64-bit arithmetic ? (y/n) ");
        Ask_Yes_or_No(lng);
        if ans = '8' then
          put("Test short division with small divisor ? (y/n) ");
          Ask_Yes_or_No(shrt);
        end if;
      end if;
      new_line;
      if lng = 'y' then
        case ans is
          when '1' => Test64_io;
          when '2' => Test64_Creation;
          when '4' => Test64_Addition;
          when '8' => Test64_Division;
          when 'A' => Multiplications64_and_Divisions_on_Randoms;
          when 'B' => Multiplication64_Timer;
          when 'C' => Division64_Remainder_Timer;
          when others => null;
        end case;
      else
        case ans is
          when '1' => Test_io;
          when '2' => Test_Creation;
          when '3' => Test_Comparison;
          when '4' => Test_Addition;
          when '5' => Test_Subtraction;
          when '6' => Test_Multiplication;
          when '7' => Test_Exponentiation;
          when '8' => Test_Division(shrt='y');
          when '9' => Additions_and_Subtractions_on_Randoms;
          when 'A' => Multiplications_and_Divisions_on_Randoms;
          when 'B' => Multiplication_Timer;
          when 'C' => Division_Remainder_Timer;
          when others => null;
        end case;
      end if;
    end loop;
  end Main;

end Test_Natural_Numbers;
