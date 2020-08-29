with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;     use Standard_Complex_Numbers_Polar;
with Multprec_Natural_Coefficients;
with Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers_io;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Multprec_Complex_Numbers_Polar;     use Multprec_Complex_Numbers_Polar;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;

package body Test_Complex_Numbers is

  procedure Test_Standard_io is

    use Standard_Complex_Numbers,Standard_Complex_Numbers_io;
    c : Complex_Number;

  begin
    new_line;
    put_line("Testing input/output for standard complex numbers.");
    new_line;
    put("Give a complex number c : "); get(c);
    put("-> c : "); put(c); new_line;
    put("-> 1/c : "); put(1.0/c); new_line;
    put("-> 1/c : "); put(Create(1.0)/c); new_line;
  end Test_Standard_io;

  procedure Test_Multprec_io is

    use Multprec_Complex_Numbers,Multprec_Complex_Numbers_io;
    c : Complex_Number;

  begin
    new_line;
    put_line("Testing input/output for multi-precision complex numbers.");
    new_line;
    put("Give a complex number c : "); get(c);
    put("-> c : "); put(c); new_line;
  end Test_Multprec_io;

  procedure Test_Standard_Roots is

    use Standard_Complex_Numbers,Standard_Complex_Numbers_io;

    d,k : natural32 := 0;
    a,c,prod : Complex_Number;
    ans : character;

  begin
    new_line;
    put_line("Solving x^d - c = 0, with c a standard complex number.");
    new_line;
    put("Give the degree d : "); get(d);
    put("Give the constant c : "); get(c);
    loop
      put("Which root do you want ? "); get(k);
      a := Root(c,d,k);
      put("The root is "); put(a); new_line;
      prod := a;
      for j in 2..d loop
        prod := prod*a;
      end loop;
      put("root^d  =   "); put(prod); new_line;
      if Equal(prod,c)
       then put_line("root^d = c, test is successful.");
       else put_line("root^d /= c, bug detected? ");
            put("Difference : "); put(prod-c); new_line;
      end if;
      put("Do you want other roots ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Standard_Roots;

  procedure Test_Multprec_Roots is

    use Standard_Complex_Numbers,Standard_Complex_Numbers_io;
    use Multprec_Complex_Numbers,Multprec_Complex_Numbers_io;

    d,k,size : natural32 := 0;
    a,c,prod,diff : Multprec_Complex_Numbers.Complex_Number;
    stc,sta : Standard_Complex_Numbers.Complex_Number;
    ans : character;

  begin
    new_line;
    put_line("Solving x^d - c = 0, with c a multi-precision complex number.");
    new_line;
    put("Give the degree d : "); get(d);
    put("Give the constant c : "); get(c);
    stc := Round(c);
    put("Give size of results : "); get(size);
    Set_Size(c,size);
    loop
      put("Which root do you want ? "); get(k);
      a := Root(c,d,k);
      sta := Root(stc,d,k);
      put("The root is "); put(a); new_line;
      put("-> Standard "); put(sta); new_line;
      Copy(a,prod);
      for j in 2..d loop
        Mul(prod,a);
      end loop;
      put("root^d  =   "); put(prod); new_line;
      if Equal(prod,c) then
        put_line("root^d = c, test is successful.");
      else
        put_line("root^d /= c, bug detected? ");
        diff := prod-c;
        put("Difference : "); put(diff); new_line;
        Clear(diff);
      end if;
      Clear(prod); Clear(a);
      put("Do you want other roots ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Test_Multprec_Roots;

  procedure Standard_Random_Addition_and_Subtraction is

    use Standard_Complex_Numbers,Standard_Complex_Numbers_io;

    n1,n2,sum1,sum2 : Complex_Number;

    procedure Report_Bug is
    begin
      new_line;
      put("  n1 : "); put(n1); new_line;
      put("  n2 : "); put(n2); new_line;
    end Report_Bug;

  begin
    n1 := Random;
    n2 := Random;
    sum1 := n1+n2;
    sum2 := sum1-n2;
    if Equal(sum2,n1) then
      put("n1+n2-n2 okay");
    else
      put("n1+n2-n2 Bug?"); Report_Bug;
      put("diff : "); put(sum2-n1); new_line;
    end if;
    Add(sum2,n2);
    if Equal(sum2,sum1) then
      put("  Add okay");
    else
      put("  Add Bug?"); Report_Bug;
      put("diff : "); put(sum2-sum1); new_line;
    end if;
    Sub(sum2,n1);
    if Equal(sum2,n2) then
      put("  Sub okay"); new_line;
    else
      put("  Sub Bug?"); Report_Bug;
      put("diff : "); put(sum2-n2); new_line;
    end if;
  exception
    when CONSTRAINT_ERROR => put_line("input caused exception:");
                             Report_Bug; raise;
  end Standard_Random_Addition_and_Subtraction;

  procedure Standard_Additions_and_Subtractions_on_Randoms is

    nb : natural32 := 0;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    tstart(timer);
    for i in 1..nb loop
      Standard_Random_Addition_and_Subtraction;
    end loop;
    tstop(timer);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Standard_Additions_and_Subtractions_on_Randoms;

  function Min ( a,b : natural32 ) return natural32 is

  begin
    if a <= b
     then return a;
     else return b;
    end if;
  end Min;

  function Acceptable_Tolerance
              ( sz1,sz2 : natural32 )
              return Multprec_Floating_Numbers.Floating_Number is

    use Multprec_Floating_Numbers;

    res : Floating_Number;
    the_radix : constant natural32 := Multprec_Natural_Coefficients.Radix;
    the_expo : constant natural32 := Multprec_Natural_Coefficients.Exponent;
    sz : constant natural32 := Min(sz1,sz2) + 1;
    exp : constant integer32 := (-integer32(sz))*integer32(the_expo+2);

  begin
    res := Create(integer32(the_radix),exp);
    return res;
  end Acceptable_Tolerance;

  procedure Report_Bug ( n1,n2 : Multprec_Complex_Numbers.Complex_Number ) is

    use Multprec_Complex_Numbers_io;

  begin
    new_line;
    put("  n1 : "); put(n1); new_line;
    put("  n2 : "); put(n2); new_line;
  end Report_Bug;

  procedure Compare_with_Tolerance 
              ( n1,n2,resi : in Multprec_Complex_Numbers.Complex_Number; 
                acctol : in Multprec_Floating_Numbers.Floating_Number ) is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;
    use Multprec_Complex_Numbers;
    absresi : Floating_Number := AbsVal(resi);

  begin
    if absresi < acctol
     then put(" < "); put(acctol,3); put_line(" okay");
     else put(" > "); put(acctol,3); put_line(" Bug!"); Report_Bug(n1,n2);
    end if;
    Clear(absresi);
  end Compare_with_Tolerance;

  procedure Multprec_Random_Addition_and_Subtraction
              ( sz1,sz2 : in natural32; low,upp : in integer32;
                acctol : in Multprec_Floating_Numbers.Floating_Number ) is

    use Multprec_Complex_Numbers,Multprec_Complex_Numbers_io;

    n1,n2,sum1,sum2,tmp : Multprec_Complex_Numbers.Complex_Number;

  begin
    n1 := Random(sz1,low,upp);
    n2 := Random(sz2,low,upp);
    sum1 := n1+n2;
    sum2 := sum1-n2;
    if Equal(sum2,n1) then
      put("n1+n2-n2 okay");
    else
      put("n1+n2-n2 Bug?"); 
      tmp := sum2-n1;
      put("  diff : "); put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Add(sum2,n2);
    if Equal(sum2,sum1) then
      put("  Add okay");
    else
      put("  Add Bug?");
      tmp := sum2-sum1;
      put("  diff : "); put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Sub(sum2,n1);
    if Equal(sum2,n2) then
      put("  Sub okay"); new_line;
    else
      put("  Sub Bug?");
      tmp := sum2-n2;
      put("diff : "); put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Clear(n1); Clear(n2);
    Clear(sum1); Clear(sum2);
  exception
    when CONSTRAINT_ERROR => put_line("input caused exception:");
                             Report_Bug(n1,n2); raise;
  end Multprec_Random_Addition_and_Subtraction;

  procedure Multprec_Additions_and_Subtractions_on_Randoms is

    nb,sz1,sz2 : natural32 := 0;
    low,upp : integer32 := 0;
    acctol : Multprec_Floating_Numbers.Floating_Number;
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
      Multprec_Random_Addition_and_Subtraction(sz1,sz2,low,upp,acctol);
    end loop;
    tstop(timer);
    Multprec_Floating_Numbers.Clear(acctol);
    put("Elapsed User CPU Time : ");
    print_hms(Standard_Output,Elapsed_User_Time(timer));
    new_line;
  end Multprec_Additions_and_Subtractions_on_Randoms;

  procedure Simul_Div ( x : in out Standard_Complex_Numbers.Complex_Number;
                        y : in Standard_Complex_Numbers.Complex_Number ) is

    acc,avyre,avyim,resre,resim,xre,xim,yre,yim : double_float;

    use Standard_Complex_Numbers;

  begin
    put_line("Simulated division on standard numbers");
    xre := REAL_PART(x); xim := IMAG_PART(x);
    yre := REAL_PART(y); yim := IMAG_PART(y);
    if Equal(yim,0.0) then
      put_line("imaginary part of divisor is zero"); 
      Div(xre,yre); Div(xim,yre);
      x := Create(xre,xim);
    elsif Equal(yre,0.0) then
      put_line("real part of divisor is zero");
      resre := xim/yim;
      resim := xre/yim; Min(resim);
      x := Create(resre,resim);
    else
      avyre := AbsVal(yre); avyim := AbsVal(yim);
      if avyre < avyim then
        put_line("imaginary part of divisor is larger than real part");
        acc := yre/yim;
        resre := xre*acc; Add(resre,xim);
        resim := xim*acc; Sub(resim,xre);
        Mul(acc,yre);
        Add(acc,yim);
        Div(resre,acc);
        Div(resim,acc);
      elsif avyre > avyim then
        put_line("real part of divisor is larger than imaginary");
        acc := yim/yre;
        put("yim/yre : "); put(acc); new_line;
        resre := xim*acc; Add(resre,xre);
        put("resre : "); put(resre); new_line;
        resim := xre*acc; Sub(resim,xim); Min(resim);
        put("resim : "); put(resim); new_line;
        Mul(acc,yim); Add(acc,yre);
        put("acc : "); put(acc); new_line;
        Div(resre,acc);
        put("resre : "); put(resre); new_line;
        Div(resim,acc);
        put("resim : "); put(resim); new_line;
      elsif Equal(yre,yim) then
        put_line("real part of divisor equals imaginary");
        acc := 2.0*yre;
        resre := xre + xim; Div(resre,acc);
        resim := xim - xre; Div(resim,acc);
      else -- y.RE = -y.IM then
        put_line("real part is negative of imaginar");
        acc := 2.0*yre;
        resre := xre - xim; Div(resre,acc);
        resim := xim + xre; Div(resim,acc);
      end if;
      x := Create(resre,resim);
    end if;
  end Simul_Div;

  procedure Simul_Div ( x : in out Multprec_Complex_Numbers.Complex_Number;
                        y : in Multprec_Complex_Numbers.Complex_Number ) is

    use Multprec_Floating_Numbers,Multprec_Floating_Numbers_io;
    use Multprec_Complex_Numbers;

    acc,avyre,avyim,resre,resim,xre,xim,yre,yim : Floating_Number;
    TWO : Floating_Number := Create(integer(2));

  begin
    put_line("Simulated division on multi-precision numbers");
    xre := REAL_PART(x); xim := IMAG_PART(x);
    yre := REAL_PART(y); yim := IMAG_PART(y);
    if Equal(yim,0.0) then
      put_line("imaginary part of divisor is zero");
      Div(xre,yre);
      Div(xim,yre);
      Clear(x);
      x := Create(xre,xim);
    elsif Equal(yre,0.0) then
      put_line("real part of divisor is zero");
      resre := xim/yim;
      resim := xre/yim; Min(resim);
      Clear(x); x := Create(resre,resim);
    else
      avyre := AbsVal(yre); avyim := AbsVal(yim);
      if avyre < avyim then
        put_line("imaginary part of divisor is larger than real part");
        acc := yre/yim;
        resre := xre*acc; Add(resre,xim);
        resim := xim*acc; Sub(resim,xre);
        Mul(acc,yre);
        Add(acc,yim);
        Div(resre,acc);
        Div(resim,acc);
        Clear(acc);
      elsif avyre > avyim then
        put_line("real part of divisor is larger than imaginary");
        acc := yim/yre;
        put("yim/yre : "); put(acc); new_line;
        resre := xim*acc; Add(resre,xre);
        put("resre : "); put(resre); new_line;
        resim := xre*acc; Sub(resim,xim); Min(resim);
        put("resim : "); put(resim); new_line;
        Mul(acc,yim); Add(acc,yre);
        put("acc : "); put(acc); new_line;
        Div(resre,acc);
        put("resre : "); put(resre); new_line;
        Div(resim,acc);
        put("resim : "); put(resim); new_line;
        Clear(acc);
      elsif Equal(yre,yim) then
        put_line("real part equals imaginary part");
        acc := TWO*yre;
        resre := xre + xim; Div(resre,acc);
        resim := xim - xre; Div(resim,acc);
        Clear(acc);
      else -- y.RE = -y.IM then
        put_line("real part is negative of imaginary part");
        acc := TWO*yre;
        resre := xre - xim; Div(resre,acc);
        resim := xim + xre; Div(resim,acc);
        Clear(acc);
      end if;
      Clear(avyre); Clear(avyim);
      Clear(x); x := Create(resre,resim);
      Clear(resre); Clear(resim);
    end if;
    Clear(xre); Clear(xim);
    Clear(yre); Clear(yim);
    Clear(TWO);
  end Simul_Div;

  procedure Interactive_Multiplication_and_Division is

    use Standard_Complex_Numbers,Standard_Complex_Numbers_io;
    use Multprec_Complex_Numbers,Multprec_Complex_Numbers_io;

    n1,n2,prod,quot : Multprec_Complex_Numbers.Complex_Number;
    stan1,stan2,staquot : Standard_Complex_Numbers.Complex_Number;
    ans : character;

  begin
    loop
      put("Give 1st number : "); get(n1);
      put("-> n1 : "); put(n1); new_line;
      put("Give 2nd number : "); get(n2);
      put("-> n2 : "); put(n2); new_line;
      prod := n1*n2;
      put("n1*n2 : "); put(prod); new_line;
      quot := prod/n2;
      put("(n1*n2)/n2  : "); put(quot); new_line;
      Clear(quot);
      quot := n1/n2;
      put("n1/n2 : "); put(quot); new_line;
      Copy(n1,quot);
      Div(quot,n2);
      put("n1/n2 : "); put(quot); new_line;
      Copy(n1,quot);
      put_line("Division with simulated procedure : ");
      Simul_Div(quot,n2);
      put("n1/n2 : "); put(quot); new_line;
      stan1 := Round(n1);
      stan2 := Round(n2);
      staquot := stan1;
      put_line("Shadowed division on standard numbers : ");
      Simul_Div(staquot,stan2);
      put("n1/n2 : "); put(staquot); new_line;
      Clear(n1); Clear(n2); Clear(prod); Clear(quot);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Interactive_Multiplication_and_Division;

  procedure Random_Multiplication_and_Division
               ( sz1,sz2 : in natural32; low,upp : in integer32;
                 acctol : in Multprec_Floating_Numbers.Floating_Number ) is

    use Multprec_Complex_Numbers,Multprec_Complex_Numbers_io;
    n1,n2,prod,quot,tmp : Complex_Number;

  begin
    n1 := Random(sz1,low,upp);
    n2 := Random(sz2,low,upp);
    prod := n1*n2;
    quot := prod/n2;
    if Equal(quot,n1) then
      put("n1*n2/n2 okay");
    else
      put("n1*n2/n2 Bug?");
      tmp := quot-n1;
      put("  Diff : "); put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Mul(quot,n2);
    if Equal(prod,quot) then
      put("  Mul okay");
    else
      put("  Mul Bug?");
      tmp := quot-prod;
      put("  Diff : "); put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Div(prod,n1);
    if Equal(prod,n2) then
      put("  Div okay"); new_line;
    else
      put("  Div Bug?");
      tmp := prod-n2;
      put("  Diff : "); put(tmp);
      Compare_with_Tolerance(n1,n2,tmp,acctol);
      Clear(tmp);
    end if;
    Clear(n1); Clear(n2);
    Clear(prod); Clear(quot);
  exception
    when CONSTRAINT_ERROR => put_line("input caused exception :");
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

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of standard and multi-precision "
                 & "complex numbers.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. Exit this program.                                     ");
      put_line("  1. Input/Output of standard complex numbers.              ");
      put_line("  2. Addition/subtraction on random standard numbers.       ");
      put_line("  3. Roots of unity of standard complex numbers.            ");
      put_line("  4.                of multi-precision complex numbers.     ");
      put_line("  5. Input/Output of multi-precision complex numbers.       ");
      put_line("  6. Addition/subtraction on random multi-precision numbers.");
      put_line("  7. Multiplication/division/remainder on random "
                                              & "multi-precision numbers.   ");
      put_line("  8. Multiplication/division on user-given numbers.         ");
      put("Type in your choice (0,1,2,3,4,5,6,7 or 8) : "); get(ans);
      exit when (ans = '0');
      new_line;
      case ans is
        when '1' => Test_Standard_io;
        when '2' => Standard_Additions_and_Subtractions_on_Randoms;
        when '3' => Test_Standard_Roots;
        when '4' => Test_Multprec_Roots;
        when '5' => Test_Multprec_io;
        when '6' => Multprec_Additions_and_Subtractions_on_Randoms;
        when '7' => Multiplications_and_Divisions_on_Randoms;
        when '8' => Interactive_Multiplication_and_Division;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Complex_Numbers;
