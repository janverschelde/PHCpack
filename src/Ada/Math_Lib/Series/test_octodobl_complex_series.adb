with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Octo_Double_Numbers;               use Octo_Double_Numbers;
with Octo_Double_Numbers_io;            use Octo_Double_Numbers_io;
with OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers_io;       use OctoDobl_Complex_Numbers_io;
with OctoDobl_Random_Numbers;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Vectors_io;
with OctoDobl_Random_Vectors;
with OctoDobl_Complex_Series;
with OctoDobl_Complex_Series_io;        use OctoDobl_Complex_Series_io;
with OctoDobl_Complex_Random_Series;
with OctoDobl_Complex_Algebraic_Series;
with OctoDobl_Complex_Series_Norms;
with OctoDobl_Complex_Series_Functions;
with OctoDobl_Complex_Series_Transforms;

package body Test_OctoDobl_Complex_Series is

  procedure OctoDobl_Construct is

    i : constant integer := 123;
    first : constant OctoDobl_Complex_Series.Series
          := OctoDobl_Complex_Series.Create(i);
    second : OctoDobl_Complex_Series.Series(4);
    third : OctoDobl_Complex_Series.Series(4);

  begin
    text_io.put_line("The first series, the constant 123 : ");
    OctoDobl_Complex_Vectors_io.put_line(first.cff);
    text_io.put_line("The second series : ");
    second.cff(0) := OctoDobl_Complex_Numbers.Create(integer32(1));
    second.cff(1) := OctoDobl_Complex_Numbers.Create(integer32(-1));
    second.cff(2) := OctoDobl_Complex_Numbers.Create(integer32(0));
    second.cff(3) := OctoDobl_Complex_Numbers.Create(integer32(0));
    second.cff(4) := OctoDobl_Complex_Numbers.Create(integer32(0));
    OctoDobl_Complex_Vectors_io.put_line(second.cff);
    third := OctoDobl_Complex_Series.Inverse(second);
    text_io.put_line("The inverse of the second series : ");
    OctoDobl_Complex_Vectors_io.put_line(third.cff);
  end OctoDobl_Construct;

  procedure OctoDobl_Test_Creation ( degree : in integer32 ) is

    use OctoDobl_Complex_Numbers;
    use OctoDobl_Complex_Series;

    s : constant Link_to_Series := Create(1,degree);
    t,x,y,z : Link_to_Series;
    minone : constant octo_double := create(-1.0);

  begin
    put("One as series of degree "); put(degree,1); put_line(" :");
    put(s);
    t := Create(1,degree);
    t.cff(1) := Create(minone);
    put_line("The series 1 - t :"); put(t); 
    x := s/t;
    put_line("The series 1/(1-t) : "); put(x);
    y := x*t;
    put_line("Verifying multiplication with inverse : "); put(y);
    z := t*x;
    put_line("Verifying commutativity : "); put(z);
  end OctoDobl_Test_Creation;

  procedure OctoDobl_Test_Arithmetic ( degree : in integer32 ) is

    use OctoDobl_Complex_Series;

    a,b,c : Series(degree);
    ans : character;

  begin
    a := OctoDobl_Complex_Random_Series.Random_Series(degree);
    put_line("The first random series A :"); put(a);
    b := OctoDobl_Complex_Random_Series.Random_Series(degree);
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
  end OctoDobl_Test_Arithmetic;

  procedure OctoDobl_Random_Test_sqrt ( degree : in integer32 ) is

    use OctoDobl_Complex_Series;

    c : constant Series(degree)
      := OctoDobl_Complex_Random_Series.Random_Series(degree);
    ans : character;
    x,y,z : Series(degree);
 
  begin
    put("Extra output during the computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("A random series c of degree "); put(degree,1); put_line(" :");
    put(c);
    if ans = 'y'
     then x := OctoDobl_Complex_Algebraic_Series.sqrt(c,0,true);
     else x := OctoDobl_Complex_Algebraic_Series.sqrt(c,0);
    end if;
    put_line("The square root x of the random series :"); put(x);
    y := x*x;
    put_line("The square y of the square root x : "); put(y);
    z := y-c;
    put_line("The equation x*x - c :"); put(z);
  end OctoDobl_Random_Test_sqrt;

  procedure OctoDobl_Random_Test_root ( degree : in integer32 ) is

    use OctoDobl_Complex_Series;

    c : constant Series(degree)
      := OctoDobl_Complex_Random_Series.Random_Series(degree);
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
     then x := OctoDobl_Complex_Algebraic_Series.Root(c,n,i,true);
     else x := OctoDobl_Complex_Algebraic_Series.Root(c,n,i);
    end if;
    put("The root x of index "); put(i,1);
    put_line(" of the random series :"); put(x);
    y := x**integer(n);
    put_line("The n-th power y of the root x : "); put(y);
    z := y-c;
    put_line("The equation x**n - c :"); put(z);
  end OctoDobl_Random_Test_root;

  function Eval ( c : OctoDobl_Complex_Vectors.Vector;
                  x : OctoDobl_Complex_Numbers.Complex_Number )
                return OctoDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the value of the polynomial with coefficients in c at x.

    use OctoDobl_Complex_Numbers;

    res : Complex_Number := c(c'last);

  begin
    for i in reverse 0..c'last-1 loop
      res := res*x + c(i);
    end loop;
    return res;
  end Eval;

  procedure OctoDobl_Random_Test_Poly_Root ( degree : in integer32 ) is

    use OctoDobl_Complex_Numbers;
    use OctoDobl_Complex_Series;

    c : constant Series(degree)
      := OctoDobl_Complex_Random_Series.Random_Series(degree);
    n : integer32 := 0;
    ans : character;

  begin
    put("Give the degree of the polynomial : "); get(n);
    declare
      p : OctoDobl_Complex_Vectors.Vector(0..n)
        := OctoDobl_Random_Vectors.Random_Vector(0,n);
      z0 : constant Complex_Number := OctoDobl_Random_Numbers.Random1;
      pz0 : constant Complex_Number := Eval(p,z0);
      yz0 : Complex_Number;
      z,y : Series(degree);
      err : octo_double;
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
       then z := OctoDobl_Complex_Algebraic_Series.Poly_Root(p,z0,c,true);
       else z := OctoDobl_Complex_Algebraic_Series.Poly_Root(p,z0,c);
      end if;
      put_line("The series expansion of the root : "); put(z);
      y := OctoDobl_Complex_Algebraic_Series.Poly_Eval(p,z);
      put_line("The polynomial at the series expansion :"); put(y);
      put_line("The right hand side series :"); put(c);
      y := y - c;
      put_line("The error series :"); put(y);
      err := OctoDobl_Complex_Series_Norms.Max_Norm(y);
      put("Max norm of the error : "); put(err,2); new_line;
    end;
  end OctoDobl_Random_Test_Poly_Root; 

  procedure OctoDobl_Test_Conjugate ( degree : in integer32 ) is

    use OctoDobl_Complex_Series;
    use OctoDobl_Complex_Series_Norms;

    s : constant Series(degree)
      := OctoDobl_Complex_Random_Series.Random_Series(degree);
    c : constant Series(degree) := Conjugate(s);
    p,r,n,q,rq : Series(degree);

  begin
    put("A random series of degree "); put(degree,1); put_line(" :");
    put(s);
    put_line("Its conjugate : "); put(c);
    put_line("Conjugate(s)*s : "); put(c*s);
    put_line("s*Conjugate(s) : "); put(s*c);
    p := c*s;
    r := OctoDobl_Complex_Algebraic_Series.sqrt(p,0);
    put_line("The square root r of Conjugate(s)*s :"); put(r);
    n := s/r;
    put_line("The normalized series s is s/r :"); put(n);
    q := Conjugate(n)*n;
    rq := OctoDobl_Complex_Algebraic_Series.sqrt(q,0);
    put_line("The norm of the normalized series :"); put(rq);
  end OctoDobl_Test_Conjugate;

  procedure OctoDobl_Test_Norm ( degree : in integer32 ) is

    use OctoDobl_Complex_Series;
    use OctoDobl_Complex_Series_Norms;

    s : constant Series(degree)
      := OctoDobl_Complex_Random_Series.Random_Series(degree);
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
  end OctoDobl_Test_Norm;

  procedure OctoDobl_Test_Shift ( degree : in integer32 ) is

    use OctoDobl_Complex_Numbers;
    use OctoDobl_Complex_Series;
    use OctoDobl_Complex_Series_Functions;

    s : constant Series(degree)
      := OctoDobl_Complex_Random_Series.Random_Series(degree);
    zero : constant octo_double := create(0.0);
    rc : octo_double := zero;
    shifteds : Series(degree);
    cc,y,z : Complex_Number;

  begin
    put_line("Shifting the series parameter ...");
    put_line("on a random series s :"); put(s);
    skip_line;
    put("Give a real constant for the shift : "); get(rc);
    shifteds := Shift(s,rc);
    y := Eval(s,-rc);
    z := Eval(shifteds,zero);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
    new_line;
    put_line("Testing with a complex shift ...");
    put("Give a complex number for the shift : "); get(cc);
    shifteds := Shift(s,cc);
    y := Eval(s,-cc);
    z := Eval(shifteds,zero);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
  end OctoDobl_Test_Shift;

  procedure OctoDobl_Test_Transform ( degree : in integer32 ) is

    use OctoDobl_Complex_Series;
    use OctoDobl_Complex_Series_Transforms;

    s : constant Series(degree)
      := OctoDobl_Complex_Random_Series.Random_Series(degree);
    c : octo_double := create(0.0);
    sc : Series(degree);
    idx : integer32;
    mxc : octo_double;

  begin
    put_line("a random series s :"); put(s);
    skip_line;
    put("Give an octo double : "); get(c);
    put("-> your constant c : "); put(c); new_line;
    sc := OctoDobl_Complex_Series_Transforms.Scale(s,c);
    put_line("the series s scaled by c, s(c*t) :"); put(sc);
    Maximum_Coefficient_Modulus(sc,idx,mxc);
    put("the index : "); put(idx,1);
    put("  maximum modulus : "); put(mxc); new_line;
    Coefficient_Modulus_Transform(sc,idx,mxc);
    put_line("the transformed series :"); put(sc);
  end OctoDobl_Test_Transform;

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
    put_line("  9. test coefficient modulus transforms");
    put("Type 0, 1, 2, 3, 4, 5, 6, 7, 8, or 9 to select a test : ");
    Ask_Alternative(ans,"0123456789");
    if ans /= '0' then
      new_line;
      put("Give the degree of the series : "); get(degree);
    end if;
    case ans is 
      when '0' => OctoDobl_Construct;
      when '1' => OctoDobl_Test_Creation(degree);
      when '2' => OctoDobl_Test_Arithmetic(degree);
      when '3' => OctoDobl_Random_Test_Sqrt(degree);
      when '4' => OctoDobl_Random_Test_root(degree);
      when '5' => OctoDobl_Random_Test_Poly_Root(degree);
      when '6' => OctoDobl_Test_Conjugate(degree);
      when '7' => OctoDobl_Test_Norm(degree);
      when '8' => OctoDobl_Test_Shift(degree);
      when '9' => OctoDobl_Test_Transform(degree);
      when others => null;
    end case;
  end Main;

end Test_OctoDobl_Complex_Series;
