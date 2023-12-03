with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Standard_Complex_Series;
with Standard_Complex_Series_io;        use Standard_Complex_Series_io;
with Standard_Complex_Random_Series;
with Standard_Complex_Algebraic_Series;
with Standard_Complex_Series_Norms;
with Standard_Complex_Series_Functions;
with Standard_Complex_Series_Transforms;

package body Test_Standard_Complex_Series is

  procedure Standard_Construct is

    i : constant integer := 123;
    first : constant Standard_Complex_Series.Series
          := Standard_Complex_Series.Create(i);
    second : Standard_Complex_Series.Series(4);
    third : Standard_Complex_Series.Series(4);

  begin
    text_io.put_line("The first series, the constant 123 : ");
    Standard_Complex_Vectors_io.put_line(first.cff);
    text_io.put_line("The second series : ");
    second.cff(0) := Standard_Complex_Numbers.Create(integer32(1));
    second.cff(1) := Standard_Complex_Numbers.Create(integer32(-1));
    second.cff(2) := Standard_Complex_Numbers.Create(integer32(0));
    second.cff(3) := Standard_Complex_Numbers.Create(integer32(0));
    second.cff(4) := Standard_Complex_Numbers.Create(integer32(0));
    Standard_Complex_Vectors_io.put_line(second.cff);
    third := Standard_Complex_Series.Inverse(second);
    text_io.put_line("The inverse of the second series : ");
    Standard_Complex_Vectors_io.put_line(third.cff);
  end Standard_Construct;

  procedure Standard_Test_Creation ( degree : in integer32 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Series;

    s : constant Link_to_Series := Create(1,degree);
    t,x,y,z : Link_to_Series;

  begin
    put("One as series of degree "); put(degree,1); put_line(" :");
    put(s);
    t := Create(1,degree);
    t.cff(1) := Create(-1.0);
    put_line("The series 1 - t :"); put(t); 
    x := s/t;
    put_line("The series 1/(1-t) : "); put(x);
    y := x*t;
    put_line("Verifying multiplication with inverse : "); put(y);
    z := t*x;
    put_line("Verifying commutativity : "); put(z);
  end Standard_Test_Creation;

  procedure Standard_Test_Arithmetic ( degree : in integer32 ) is

    use Standard_Complex_Series;

    a,b,c : Series(degree);
    ans : character;

  begin
    a := Standard_Complex_Random_Series.Random_Series(degree);
    put_line("The first random series A :"); put(a);
    b := Standard_Complex_Random_Series.Random_Series(degree);
    put_line("The second random series B :"); put(b);
    c := a+b;
    put_line("The sum A + B :"); put(c);
    c := c-a;
    put_line("The sum A + B - A :"); put(c); 
    new_line;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      c := a*b;
      put_line("The product A*B :"); put(c);
      c := c/a;
      put_line("The product A*B/A :"); put(c);
    end if;
  end Standard_Test_Arithmetic;

  procedure Standard_Random_Test_sqrt ( degree : in integer32 ) is

    use Standard_Complex_Series;

    c : constant Series(degree)
      := Standard_Complex_Random_Series.Random_Series(degree);
    ans : character;
    x,y,z : Series(degree);
 
  begin
    put("Extra output during the computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("A random series c of degree "); put(degree,1); put_line(" :");
    put(c);
    if ans = 'y'
     then x := Standard_Complex_Algebraic_Series.sqrt(c,0,true);
     else x := Standard_Complex_Algebraic_Series.sqrt(c,0);
    end if;
    put_line("The square root x of the random series :"); put(x);
    y := x*x;
    put_line("The square y of the square root x : "); put(y);
    z := y-c;
    put_line("The equation x*x - c :"); put(z);
  end Standard_Random_Test_sqrt;

  procedure Standard_Random_Test_root ( degree : in integer32 ) is

    use Standard_Complex_Series;

    c : constant Series(degree)
      := Standard_Complex_Random_Series.Random_Series(degree);
    n,i : natural32 := 0;
    ans : character;
    x,y,z : Series(degree);
 
  begin
    put("Give the power n in x**n - c : "); get(n);
    put("Give the index i of the root : "); get(i);
    new_line;
    put("Extra output during the computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("A random series c of degree "); put(degree,1); put_line(" :");
    put(c);
    if ans = 'y'
     then x := Standard_Complex_Algebraic_Series.Root(c,n,i,true);
     else x := Standard_Complex_Algebraic_Series.Root(c,n,i);
    end if;
    put("The root x of index "); put(i,1);
    put_line(" of the random series :"); put(x);
    y := x**integer(n);
    put_line("The n-th power y of the root x : "); put(y);
    z := y-c;
    put_line("The equation x**n - c :"); put(z);
  end Standard_Random_Test_root;

  function Eval ( c : Standard_Complex_Vectors.Vector;
                  x : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the value of the polynomial with coefficients in c at x.

    use Standard_Complex_Numbers;

    res : Complex_Number := c(c'last);

  begin
    for i in reverse 0..c'last-1 loop
      res := res*x + c(i);
    end loop;
    return res;
  end Eval;

  procedure Standard_Random_Test_Poly_Root ( degree : in integer32 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Series;

    c : constant Series(degree)
      := Standard_Complex_Random_Series.Random_Series(degree);
    n : integer32 := 0;
    ans : character;

  begin
    put("Give the degree of the polynomial : "); get(n);
    declare
      p : Standard_Complex_Vectors.Vector(0..n)
        := Standard_Random_Vectors.Random_Vector(0,n);
      z0 : constant Complex_Number := Standard_Random_Numbers.Random1;
      pz0 : constant Complex_Number := Eval(p,z0);
      yz0 : Complex_Number;
      z,y : Series(degree);
      err : double_float;
    begin
      put_line("A random root : "); put(z0); new_line;
      p(0) := p(0) - pz0;
      yz0 := Eval(p,z0);
      put_line("Value at a random root : "); put(yz0);
      new_line;
      put("Extra output during the computation ? (y/n) ");
      Ask_Yes_or_No(ans);
      new_line;
      if ans = 'y'
       then z := Standard_Complex_Algebraic_Series.Poly_Root(p,z0,c,true);
       else z := Standard_Complex_Algebraic_Series.Poly_Root(p,z0,c);
      end if;
      put_line("The series expansion of the root : "); put(z);
      y := Standard_Complex_Algebraic_Series.Poly_Eval(p,z);
      put_line("The polynomial at the series expansion :"); put(y);
      put_line("The right hand side series :"); put(c);
      y := y - c;
      put_line("The error series :"); put(y);
      err := Standard_Complex_Series_Norms.Max_Norm(y);
      put("Max norm of the error :"); put(err,2); new_line;
    end;
  end Standard_Random_Test_Poly_Root; 

  procedure Standard_Test_Conjugate ( degree : in integer32 ) is

    use Standard_Complex_Series;
    use Standard_Complex_Series_Norms;

    s : constant Series(degree)
      := Standard_Complex_Random_Series.Random_Series(degree);
    c : constant Series(degree) := Conjugate(s);
    p,r,n,q,rq : Series(degree);

  begin
    put("A random series of degree "); put(degree,1); put_line(" :");
    put(s);
    put_line("Its conjugate : "); put(c);
    put_line("Conjugate(s)*s : "); put(c*s);
    put_line("s*Conjugate(s) : "); put(s*c);
    p := c*s;
    r := Standard_Complex_Algebraic_Series.sqrt(p,0);
    put_line("The square root r of Conjugate(s)*s :"); put(r);
    n := s/r;
    put_line("The normalized series s is s/r :"); put(n);
    q := Conjugate(n)*n;
    rq := Standard_Complex_Algebraic_Series.sqrt(q,0);
    put_line("The norm of the normalized series :"); put(rq);
  end Standard_Test_Conjugate;

  procedure Standard_Test_Norm ( degree : in integer32 ) is

    use Standard_Complex_Series;
    use Standard_Complex_Series_Norms;

    s : constant Series(degree)
      := Standard_Complex_Random_Series.Random_Series(degree);
    nrm : constant Series(degree) := Norm(s);
    ns : constant Series(degree) := Normalize(s);
    nrm2 : constant Series(degree) := Norm(ns);
  
  begin
    put("A random series of degree "); put(degree,1); put_line(" :");
    put(s);
    put_line("Its norm :"); put(nrm);
    put("The max-norm of the series : ");
    put(Max_Norm(s),3); new_line;
    put("The two-norm of the series : ");
    put(Two_Norm(s),3); new_line;
    put_line("The normalized series :"); put(ns);
    put_line("The norm of the normalized series :"); put(nrm2);
    put("The max-norm of the normalized series : ");
    put(Max_Norm(ns),3); new_line;
    put("The two-norm of the normalized series : ");
    put(Two_Norm(ns),3); new_line;
  end Standard_Test_Norm;

  procedure Standard_Test_Shift ( degree : in integer32 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Series;
    use Standard_Complex_Series_Functions;

    s : constant Series(degree)
      := Standard_Complex_Random_Series.Random_Series(degree);
    rc : double_float := 0.0;
    shifteds : Series(degree);
    cc,y,z : Complex_Number;

  begin
    put_line("Shifting the series parameter ...");
    put_line("on a random series s :"); put(s);
    put("Give a real constant for the shift : "); get(rc);
    shifteds := Shift(s,rc);
    y := Eval(s,-rc);
    z := Eval(shifteds,0.0);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
    new_line;
    put_line("Testing with a complex shift ...");
    put("Give a complex number for the shift : "); get(cc);
    shifteds := Shift(s,cc);
    y := Eval(s,-cc);
    z := Eval(shifteds,0.0);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
  end Standard_Test_Shift;

  procedure Standard_Test_Power ( degree : in integer32 ) is

    use Standard_Complex_Series;

    s : constant Series(degree)
      := Standard_Complex_Random_Series.Random_Series(degree);
    xs : constant Link_to_Series := new Series'(s);
    ls : Link_to_Series := new Series'(s);
    ps : Series(degree);
    pwr : integer32 := 0;

  begin
    put_line("Taking powers of a series ...");
    put_line("on a random series s :"); put(s);
    ps := s*s;
    put_line(" s**2 : "); put(ps);
    Mul(ls,ls);
    put_line("ls**2 : "); put(ls);
    put("Give a power : "); get(pwr);
    ps := s**integer(pwr);
    put(" s**"); put(pwr,1); put_line(" :"); put(ps);
    Power(ls,xs,integer(pwr));
    put("ls**"); put(pwr,1); put_line(" :"); put(ls);
  end Standard_Test_Power;

  procedure Standard_Test_Transform ( degree : in integer32 ) is

    use Standard_Complex_Series;
    use Standard_Complex_Series_Transforms;

    s : constant Series(degree)
      := Standard_Complex_Random_Series.Random_Series(degree);
    c : double_float := 0.0;
    sc : Series(degree);
    idx : integer32;
    mxc : double_float;

  begin
    put_line("a random series s :"); put(s);
    put("Give a double float : "); get(c);
    put("-> your constant c : "); put(c); new_line;
    sc := Standard_Complex_Series_Transforms.Scale(s,c);
    put_line("the series s scaled by c, s(c*t) :"); put(sc);
    Maximum_Coefficient_Modulus(sc,idx,mxc);
    put("the index : "); put(idx,1);
    put("  maximum modulus : "); put(mxc); new_line;
    Coefficient_Modulus_Transform(sc,idx,mxc);
    put_line("the transformed series :"); put(sc);
  end Standard_Test_Transform;

  procedure Main is

    ans : character;
    degree : integer32 := 0;

  begin
    new_line;
    put_line("MENU with testing operations :");
    put_line("  0. test the basic construct methods");
    put_line("  1. test the computation of 1/(1-t) for any degree");
    put_line("  2. test arithmetic");
    put_line("  3. square root of a random series");
    put_line("  4. p-th root of a random series");
    put_line("  5. series expansion of a root of a polynomial");
    put_line("  6. test complex conjugate of a series");
    put_line("  7. test the norm of a series");
    put_line("  8. test shift of series parameter");
    put_line("  9. test computation of powers");
    put_line("  A. test coefficient modulus transforms");
    put("Type 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, or A to select a test : ");
    Ask_Alternative(ans,"0123456789A");
    if ans /= '0' then
      new_line;
      put("Give the degree of the series : "); get(degree);
    end if;
    case ans is 
      when '0' => Standard_Construct;
      when '1' => Standard_Test_Creation(degree);
      when '2' => Standard_Test_Arithmetic(degree);
      when '3' => Standard_Random_Test_sqrt(degree);
      when '4' => Standard_Random_Test_root(degree);
      when '5' => Standard_Random_Test_Poly_Root(degree);
      when '6' => Standard_Test_Conjugate(degree);
      when '7' => Standard_Test_Norm(degree);
      when '8' => Standard_Test_Shift(degree);
      when '9' => Standard_Test_Power(degree);
      when 'A' => Standard_Test_Transform(degree);
      when others => null;
    end case;
  end Main;

end Test_Standard_Complex_Series;
