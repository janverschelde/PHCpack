with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Timing_Package;                     use Timing_Package;
with Multprec_Natural_Numbers;
with Multprec_Natural_Numbers_io;
with Multprec_Integer_Numbers_io;
with Multprec_Natural64_Numbers;
with Multprec_Natural64_Numbers_io;
with Multprec_Integer64_Numbers;
with Multprec_Integer64_Numbers_io;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;

package body Test_Integer_Numbers is

  procedure Test_Creation is

    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;

    i1 : integer32 := 0;
    i2 : Integer_Number;
    ans : character;

  begin
    put_line("Testing the creation of an integer number.");
    loop
      put("Give a standard integer number : "); get(i1);
      i2 := Create(i1);
      put("-> as integer number : "); put(i2); new_line;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Creation;

  procedure Test64_Creation is

    use Multprec_Integer64_Numbers,Multprec_Integer64_Numbers_io;

    i1 : integer64 := 0;
    i2 : Integer_Number;
    ans : character;

  begin
    put_line("Testing the creation of an integer number.");
    loop
      put("Give a standard integer number : "); get(i1);
      put("-> your given integer : "); put(i1); new_line;
      i2 := Create64(i1);
      put("--> as integer number : "); put(i2); new_line;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Creation;

  procedure Test_io is

    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;

    ans : character;
    i : Integer_Number;
    nd,nc : natural32 := 0;

  begin
    put_line("Testing the input/output operations.");
    loop
      put("Give a number : "); get(i);
      put("-> your number : "); put(i); new_line;
      nd := Decimal_Places(i);
      put("#decimal places : "); put(nd,1); new_line;
      nc := nd;
      if Negative(i)
       then nc := nc + 1;
      end if;
      declare
        s : string(1..integer(nc));
        i2 : Integer_Number;
      begin
        put(s,i);
        put("String representation : "); put(s); new_line;
        get(s,i2);
        put(" after parsing string : "); put(i2); new_line;
      end;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_io;

  procedure Test64_io is

    use Multprec_Integer64_Numbers,Multprec_Integer64_Numbers_io;

    ans : character;
    i : Integer_Number;
    nd,nc : natural32 := 0;

  begin
    put_line("Testing the input/output operations.");
    loop
      put("Give a number : "); get(i);
      put("-> your number : "); put(i); new_line;
      nd := Decimal_Places(i);
      put("#decimal places : "); put(nd,1); new_line;
      nc := nd;
      if Negative(i)
       then nc := nc + 1;
      end if;
      declare
        s : string(1..integer(nc));
        i2 : Integer_Number;
      begin
        put(s,i);
        put("String representation : "); put(s); new_line;
        get(s,i2);
        put(" after parsing string : "); put(i2); new_line;
      end;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_io;

  procedure Test_Sign ( i : in Multprec_Integer_Numbers.Integer_Number ) is

    use Multprec_Integer_Numbers;

  begin
    if Multprec_Integer_Numbers.Positive(i)
     then put("This number is positive,");
     else put("This number is not positive,");
    end if;
    if Negative(i)
     then put(" is negative ");
     else put(" is not negative ");
    end if;
    put("and its sign is ");
    if Sign(i) > 0 then
      put("+"); 
    elsif Sign(i) < 0 then
      put("-"); 
    else
      put("0");
    end if;
    put_line(".");
  end Test_Sign;

  procedure Test_Compare
              ( i1,i2 : in Multprec_Integer_Numbers.Integer_Number ) is

    use Multprec_Integer_Numbers;

  begin
    if Equal(i1,i2)
     then put_line("The numbers are equal.");
     else put_line("The numbers are different.");
    end if;
    if i1 < i2
     then put_line("First is less than second.");
     else put_line("First not less than second.");
    end if;
    if i1 > i2
     then put_line("First is greater than second.");
     else put_line("First is not greater than second.");
    end if;
  end Test_Compare;

  procedure Zero_Test ( i : Multprec_Integer_Numbers.Integer_Number ) is
  begin
    if Multprec_Integer_Numbers.Equal(i,0)
     then put_line(" equals zero");
     else put_line(" is different from zero");
    end if;
  end Zero_Test;

  procedure Test_Comparison is

    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;

    ans : character;
    i1,i2 : Integer_Number;

  begin 
    put_line("Testing the comparison operations.");
    loop
      put("Give 1st number i1 : "); get(i1);
      put("-> i1 : "); put(i1); 
      Zero_Test(i1);
      Test_Sign(i1);
      put("Give 2nd number i2 : "); get(i2);
      put("-> i2 : "); put(i2);
      Zero_Test(i2);
      Test_Sign(i2);
      Test_Compare(i1,i2);
      Copy(i1,i2);
      put_line("Tests after copying : ");
      Test_Compare(i1,i2);
      Div(i1,10);
      put_line("After dividing i1 by 10 :");
      put(" i1 : "); put(i1); new_line;
      put(" i2 : "); put(i2); new_line;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Comparison;

  procedure Test_Addition is

    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;
 
    ans : character;
    i1,i2,sum1,sum2 : Integer_Number;

  begin
    put_line("Testing the addition operations.");
    loop
      put("Give 1st number : "); get(i1);
      put("-> your 1st number i1 : "); put(i1); new_line;
      put("Give 2nd number : "); get(i2);
      put("-> your 2nd number i2 : "); put(i2); new_line;
      sum1 := i1+i2;
      put("i1+i2 : "); put(sum1); new_line;
      sum2 := i2+i1;
      put("i2+i1 : "); put(sum2); new_line;
      if Equal(sum1,sum2)
       then put_line("Test on commutativity is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Addition;

  procedure Test64_Addition is

    use Multprec_Integer64_Numbers,Multprec_Integer64_Numbers_io;
 
    ans : character;
    i1,i2,sum1,sum2 : Integer_Number;

  begin
    put_line("Testing the addition operations.");
    loop
      put("Give 1st number : "); get(i1);
      put("-> your 1st number i1 : "); put(i1); new_line;
      put("Give 2nd number : "); get(i2);
      put("-> your 2nd number i2 : "); put(i2); new_line;
      sum1 := i1+i2;
      put("i1+i2 : "); put(sum1); new_line;
      sum2 := i2+i1;
      put("i2+i1 : "); put(sum2); new_line;
      if Equal(sum1,sum2)
       then put_line("Test on commutativity is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Addition;

 -- function Mult_by_Add
 --            ( i1 : Multprec_Integer_Numbers.Integer_Number; i2 : integer ) 
 --            return Multprec_Integer_Numbers.Integer_Number is
 --
  -- DESCRIPTION :
  --   Does the multiplication by adding up i1 to itself as many times
  --   as the number i2.  Only to be used as test of course.
 --
 --   use Multprec_Integer_Numbers;
 --
 --   res : Integer_Number;
 --   n : natural := 0;
 --
 -- begin
 --   if i2 = 0 then
 --     return res;
 --   else
 --     Copy(i1,res);
 --     if i2 < 0
 --      then n := -i2;
 --      else n := i2;
 --     end if;
 --     for i in 1..n-1 loop
 --       Add(res,i1);
 --     end loop;
 --     if i2 < 0
 --      then Min(res);
 --     end if;
 --     return res;
 --   end if;
 -- end Mult_by_Add;

  function Mult_by_Add
             ( i1,i2 : Multprec_Integer_Numbers.Integer_Number )
             return Multprec_Integer_Numbers.Integer_Number is

  -- DESCRIPTION :
  --   Does the multiplication by adding up n1 to itself as many times
  --   as the number i2.  Only to be used as test of course.
  --   This can be quite time consuming as i2 gets large.

    use Multprec_Natural_Numbers,Multprec_Integer_Numbers;

    res : Integer_Number;
    cnt,tot : Natural_Number;

  begin
    if Equal(i2,0) then
      return res;
    else
      Copy(i1,res);
      cnt := Create(natural32(1));
      tot := Unsigned(i2);
      while not Equal(cnt,tot) loop
        Add(res,i1);
        Add(cnt,1);
      end loop;
      Clear(cnt);
      if Negative(i2)
       then Min(res);
      end if;
      return res;
    end if;
  end Mult_by_Add;

  procedure Test_Multiplication is

    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;
 
    ans : character;
    i1,i2,prod1,prod2,prod3 : Integer_Number;
   -- i2 : integer;

  begin
    put_line("Testing the multiplication operations.");
    loop
      put("Give 1st number : "); get(i1);
      put("-> your 1st number i1 : "); put(i1); new_line;
      put("Give 2nd number : "); get(i2);
      put("-> your 2nd number i2 : "); put(i2); new_line;
      prod1 := i1*i2;
      put("Product i1*i2 : "); put(prod1); new_line;
      prod2 := i2*i1;
      put("Product i2*i1 : "); put(prod2); new_line;
      if Equal(prod1,prod2)
       then put_line("Test on commutativity is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want multiplication by addition ? (y/n) "); get(ans);
      if ans = 'y' then
        put_line("Testing the multiplication by addition.  Be patient...");
        prod3 := Mult_by_Add(i1,i2);
        put("After adding "); put(i2); put(" times : "); put(prod3);
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
    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;

    ans : character;
    e1,e2 : Natural_Number;
    i,exp1,exp2,prod,expo : Integer_Number;

  begin
    put_line("Testing the exponentiation operations.");
    loop
      put("Give a number : "); get(i);
      put("-> your number i : "); put(i); new_line;
      put("Give 1st exponent : "); get(e1);
      put("-> your 1st exponent e1 : "); put(e1); new_line;
      exp1 := i**e1;
      put("i**e1 : "); put(exp1); new_line;
      put("Give 2nd exponent : "); get(e2);
      put("-> your 2nd exponent e2 : "); put(e2); new_line;
      exp2 := i**e2;
      put("i**e2 : "); put(exp2); new_line;
      prod := exp1*exp2;
      put("(i**e1)*(i**e2) : "); put(prod); new_line;
      expo := i**(e1+e2);
      put("i**(e1+e2)      : "); put(expo); new_line;
      if Equal(prod,expo)
       then put_line("Test of exponentiation is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Exponentiation;

  procedure Test64_Exponentiation is

    use Multprec_Natural64_Numbers,Multprec_Natural64_Numbers_io;
    use Multprec_Integer64_Numbers,Multprec_Integer64_Numbers_io;

    ans : character;
    e1,e2 : Natural_Number;
    i,exp1,exp2,prod,expo : Integer_Number;

  begin
    put_line("Testing the exponentiation operations.");
    loop
      put("Give a number : "); get(i);
      put("-> your number i : "); put(i); new_line;
      put("Give 1st exponent : "); get(e1);
      put("-> your 1st exponent e1 : "); put(e1); new_line;
      exp1 := i**e1;
      put("i**e1 : "); put(exp1); new_line;
      put("Give 2nd exponent : "); get(e2);
      put("-> your 2nd exponent e2 : "); put(e2); new_line;
      exp2 := i**e2;
      put("i**e2 : "); put(exp2); new_line;
      prod := exp1*exp2;
      put("(i**e1)*(i**e2) : "); put(prod); new_line;
      expo := i**(e1+e2);
      put("i**(e1+e2)      : "); put(expo); new_line;
      if Equal(prod,expo)
       then put_line("Test of exponentiation is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Exponentiation;

  procedure Test_Subtraction is

    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;

    ans : character;
    i1,i2,diff : Integer_Number;
   -- i2 : integer;

  begin
    put_line("Testing the subtraction operations.");
    loop
      put("Give 1st number : "); get(i1);
      put("-> your 1st number i1 : "); put(i1); new_line;
      put("Give 2nd number : "); get(i2);
      put("-> your 2nd number i2 : "); put(i2); new_line;
      diff := i1-i2;
      put("i1 - i2 : "); put(diff); new_line;
      Add(diff,i2);
      put("(i1-i2)+i2 : "); put(diff); new_line;
      if Equal(diff,i1)
       then put_line("Test of subtraction is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Subtraction;

  procedure Test64_Subtraction is

    use Multprec_Integer64_Numbers,Multprec_Integer64_Numbers_io;

    ans : character;
    i1,i2,diff : Integer_Number;
   -- i2 : integer;

  begin
    put_line("Testing the subtraction operations.");
    loop
      put("Give 1st number : "); get(i1);
      put("-> your 1st number i1 : "); put(i1); new_line;
      put("Give 2nd number : "); get(i2);
      put("-> your 2nd number i2 : "); put(i2); new_line;
      diff := i1-i2;
      put("i1 - i2 : "); put(diff); new_line;
      Add(diff,i2);
      put("(i1-i2)+i2 : "); put(diff); new_line;
      if Equal(diff,i1)
       then put_line("Test of subtraction is successful.");
       else put_line("Failure, bug detected.");
      end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Subtraction;

 -- procedure Divide10 ( i : in Multprec_Integer_Numbers.Integer_Number ) is
 --
 --   use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;
 --
  -- DESCRIPTION :
  --   Checks whether the number i is divisible by 1..10.
 --
 --   quot,prod : Integer_Number;
 --   rest : integer := 0;
 --
 -- begin
 --   put("i : "); put(i); new_line;
 --   for j in 1..10 loop
 --     rest := Rmd(i,j);
 --     quot := i/j;
 --     if rest = 0
 --      then put("Divisible by "); put(j,1);
 --      else put("Not divisible by "); put(j,1);
 --     end if;
 --     put("  rest : "); put(rest,1); new_line;
 --     put("quotient : "); put(quot); new_line;
 --     prod := quot*j + rest;
 --     if Equal(prod,i)
 --      then put_line("Test on Remainder/Division is successful.");
 --      else put_line("Failure, bug detected.");
 --     end if;
 --   end loop;
 -- end Divide10;

  procedure Test_Division is

    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;

    ans : character;
    i1,quot,prod : Integer_Number;
    i2,rest : Integer_Number;
   -- i2,rest : integer;

  begin
    put_line("Testing the division operations.");
    loop
      put("Give 1st number : "); get(i1);
      put("-> your 1st number i1 : "); put(i1); new_line;
      put("Give 2nd number : "); get(i2);
      put("-> your 2nd number i2 : "); put(i2); new_line;
      prod := i1*i2;
      put("i1*i2 : "); put(prod); new_line;
      quot := prod/i2; rest := Rmd(prod,i2);
      put("(i1*i2)/i2 : "); put(quot); new_line;
      put("Remainder : "); put(rest); new_line;
      if Equal(quot,i1) and Equal(rest,0) -- rest = 0
       then put_line("Test of division is successful.");
       else put_line("Failure, bug detected.");
      end if;
      Div(i1,i2,quot,rest);
      put("i1/i2 : "); put(quot); new_line;
      put("rest : "); put(rest); new_line;
      prod := quot*i2 + rest;
      if Equal(prod,i1)
       then put_line("Test of division/remainder computation is successful.");
       else put_line("Failure, bug detected.");
      end if;
     -- if i2 <= 10
     --  then Divide10(i1);
     -- end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Division;

  procedure Test64_Division is

    use Multprec_Integer64_Numbers,Multprec_Integer64_Numbers_io;

    ans : character;
    i1,quot,prod : Integer_Number;
    i2,rest : Integer_Number;
   -- i2,rest : integer;

  begin
    put_line("Testing the division operations.");
    loop
      put("Give 1st number : "); get(i1);
      put("-> your 1st number i1 : "); put(i1); new_line;
      put("Give 2nd number : "); get(i2);
      put("-> your 2nd number i2 : "); put(i2); new_line;
      prod := i1*i2;
      put("i1*i2 : "); put(prod); new_line;
      quot := prod/i2; rest := Rmd(prod,i2);
      put("(i1*i2)/i2 : "); put(quot); new_line;
      put("Remainder : "); put(rest); new_line;
      if Equal(quot,i1) and Equal(rest,0) -- rest = 0
       then put_line("Test of division is successful.");
       else put_line("Failure, bug detected.");
      end if;
      Div(i1,i2,quot,rest);
      put("i1/i2 : "); put(quot); new_line;
      put("rest : "); put(rest); new_line;
      prod := quot*i2 + rest;
      if Equal(prod,i1)
       then put_line("Test of division/remainder computation is successful.");
       else put_line("Failure, bug detected.");
      end if;
     -- if i2 <= 10
     --  then Divide10(i1);
     -- end if;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test64_Division;

  procedure Random_Addition_and_Subtraction ( sz1,sz2 : in natural32 ) is

    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;

  -- DESCRIPTION :
  --   Three tests are performed:
  --   1) n1+n2-n2 = n1, with "+" and "-".
  --   2) Add(n1,n2) is the same as n1 := n1+n2?
  --   3) Sub(n1+n2,n1) leads to n2?

    n1,n2,sum1,sum2 : Integer_Number;

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

  procedure Random64_Addition_and_Subtraction ( sz1,sz2 : in natural32 ) is

    use Multprec_Integer64_Numbers,Multprec_Integer64_Numbers_io;

    n1,n2,sum1,sum2 : Integer_Number;

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
  end Random64_Addition_and_Subtraction;

  procedure Additions_and_Subtractions_on_Randoms ( sixtyfour : in boolean ) is

    nb,sz1,sz2 : natural32 := 0;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    if sixtyfour then
      tstart(timer);
      for i in 1..nb loop
        Random64_Addition_and_Subtraction(sz1,sz2);
      end loop;
      tstop(timer);
    else
      tstart(timer);
      for i in 1..nb loop
        Random_Addition_and_Subtraction(sz1,sz2);
      end loop;
      tstop(timer);
    end if;
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Additions_and_Subtractions_on_Randoms;

  procedure Random_Multiplication_and_Division ( sz1,sz2 : in natural32 ) is

    use Multprec_Integer_Numbers,Multprec_Integer_Numbers_io;

    n1,n2,prod1,prod2,quot1,quot2,quot3,rest1,rest2 : Integer_Number;

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

    use Multprec_Integer64_Numbers,Multprec_Integer64_Numbers_io;

    n1,n2,prod1,prod2,quot1,quot2,quot3,rest1,rest2 : Integer_Number;

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

  procedure Multiplications_and_Divisions_on_Randoms ( six4 : in boolean ) is

    nb,sz1,sz2 : natural32 := 0;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    put("Give the size of the 1st number : "); get(sz1);
    put("Give the size of the 2nd number : "); get(sz2);
    if six4 then
      tstart(timer);
      for i in 1..nb loop
        Random64_Multiplication_and_Division(sz1,sz2);
      end loop;
      tstop(timer);
    else
      tstart(timer);
      for i in 1..nb loop
        Random_Multiplication_and_Division(sz1,sz2);
      end loop;
      tstop(timer);
    end if;
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Multiplications_and_Divisions_on_Randoms;

  procedure Main is

    ans,lng : character;

  begin
    new_line;
    put_line("Interactive testing of multi-precision integer numbers.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit program      1. Input/Output     2. Creation      ");
      put_line("  3. Comparison/Copy   4. Addition         5. Subtraction   ");
      put_line("  6. Multiplication    7. Exponentiation   8. Division      ");
      put_line("  9. Addition/subtraction on randomly generated numbers.    ");
      put_line("  A. Multiplication/division/remainder on random numbers.   ");
      put("Type in your choice (0,1,2,3,4,5,6,7,8,9 or A) : ");
      Ask_Alternative(ans,"0123456789A");
      exit when (ans = '0');
      if ans /= '3' and ans /= '6' then
        new_line;
        put("Use 64-bit arithmetic ? (y/n) "); Ask_Yes_or_No(lng);
      end if;
      new_line;
      if lng = 'y' then
        case ans is
          when '1' => Test64_io;
          when '2' => Test64_Creation;
          when '4' => Test64_Addition;
          when '5' => Test64_Subtraction;
          when '7' => Test64_Exponentiation;
          when '8' => Test64_Division;
          when '9' => Additions_and_Subtractions_on_Randoms(true);
          when 'A' => Multiplications_and_Divisions_on_Randoms(true);
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
          when '8' => Test_Division;
          when '9' => Additions_and_Subtractions_on_Randoms(false);
          when 'A' => Multiplications_and_Divisions_on_Randoms(false);
          when others => null;
        end case;
      end if;
    end loop;
  end Main;

end Test_Integer_Numbers;
