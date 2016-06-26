with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with Standard_Complex_Vectors;
with Standard_Dense_Series;
with Standard_Dense_Series_io;           use Standard_Dense_Series_io;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_io;           use DoblDobl_Dense_Series_io;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_io;           use QuadDobl_Dense_Series_io;
with Standard_Random_Series;             use Standard_Random_Series;
with DoblDobl_Random_Series;             use DoblDobl_Random_Series;
with QuadDobl_Random_Series;             use QuadDobl_Random_Series;
with Standard_Algebraic_Series;
with DoblDobl_Algebraic_Series;
with QuadDobl_Algebraic_Series;
with Standard_Dense_Series_Norms;

procedure ts_series is

-- DESCRIPTION :
--   Tests the operations on truncated power series.

  procedure Standard_Test_Creation ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Verifies that 1/(1-t) = 1 + t + t^2 + ...
  --   for a truncated power series with coefficients
  --   in standard double precision.

    use Standard_Complex_Numbers;
    use Standard_Dense_Series;

    s : Series := Create(1,order);
    t : Series := s;
    x,y,z : Series;

  begin
    put("One as series of order "); put(order,1); put_line(" :");
    put(s);
    t.cff(1) := Create(-1.0);
    put_line("The series 1 - t :"); put(t); 
    x := s/t;
    put_line("The series 1/(1-t) : "); put(x);
    y := x*t;
    put_line("Verifying multiplication with inverse : "); put(y);
    z := t*x;
    put_line("Verifying commutativity : "); put(z);
  end Standard_Test_Creation;

  procedure DoblDobl_Test_Creation ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Verifies that 1/(1-t) = 1 + t + t^2 + ...
  --   for a truncated power series with coefficients
  --   in double double precision.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Dense_Series;

    s : Series := Create(1,order);
    t : Series := s;
    x,y,z : Series;
    minone : constant double_double := create(-1.0);

  begin
    put("One as series of order "); put(order,1); put_line(" :");
    put(s);
    t.cff(1) := Create(minone);
    put_line("The series 1 - t :"); put(t); 
    x := s/t;
    put_line("The series 1/(1-t) : "); put(x);
    y := x*t;
    put_line("Verifying multiplication with inverse : "); put(y);
    z := t*x;
    put_line("Verifying commutativity : "); put(z);
  end DoblDobl_Test_Creation;

  procedure QuadDobl_Test_Creation ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Verifies that 1/(1-t) = 1 + t + t^2 + ...
  --   for a truncated power series with coefficients
  --   in double double precision.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Dense_Series;

    s : Series := Create(1,order);
    t : Series := s;
    x,y,z : Series;
    minone : constant quad_double := create(-1.0);

  begin
    put("One as series of order "); put(order,1); put_line(" :");
    put(s);
    t.cff(1) := Create(minone);
    put_line("The series 1 - t :"); put(t); 
    x := s/t;
    put_line("The series 1/(1-t) : "); put(x);
    y := x*t;
    put_line("Verifying multiplication with inverse : "); put(y);
    z := t*x;
    put_line("Verifying commutativity : "); put(z);
  end QuadDobl_Test_Creation;

  procedure Standard_Random_Test_sqrt ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given order
  --   and tests the square root computation,
  --   in standard double precision.

    use Standard_Dense_Series;

    c : constant Series := Random_Series(order);
    ans : character;
    x,y,z : Series;
 
  begin
    put("Extra output during the computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("A random series c of order "); put(order,1); put_line(" :");
    put(c);
    if ans = 'y'
     then x := Standard_Algebraic_Series.sqrt(c,0,true);
     else x := Standard_Algebraic_Series.sqrt(c,0);
    end if;
    put_line("The square root x of the random series :"); put(x);
    y := x*x;
    put_line("The square y of the square root x : "); put(y);
    z := y-c;
    put_line("The equation x*x - c :"); put(z);
  end Standard_Random_Test_sqrt;

  procedure DoblDobl_Random_Test_sqrt ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given order
  --   and tests the square root computation,
  --   in double double precision.

    use DoblDobl_Dense_Series;

    c : constant Series := Random_Series(order);
    ans : character;
    x,y,z : Series;
 
  begin
    put("Extra output during the computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("A random series c of order "); put(order,1); put_line(" :");
    put(c);
    if ans = 'y'
     then x := DoblDobl_Algebraic_Series.sqrt(c,0,true);
     else x := DoblDobl_Algebraic_Series.sqrt(c,0);
    end if;
    put_line("The square root x of the random series :"); put(x);
    y := x*x;
    put_line("The square y of the square root x : "); put(y);
    z := y-c;
    put_line("The equation x*x - c :"); put(z);
  end DoblDobl_Random_Test_sqrt;

  procedure QuadDobl_Random_Test_sqrt ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given order
  --   and tests the square root computation,
  --   in quad double precision.

    use QuadDobl_Dense_Series;

    c : constant Series := Random_Series(order);
    ans : character;
    x,y,z : Series;
 
  begin
    put("Extra output during the computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("A random series c of order "); put(order,1); put_line(" :");
    put(c);
    if ans = 'y'
     then x := QuadDobl_Algebraic_Series.sqrt(c,0,true);
     else x := QuadDobl_Algebraic_Series.sqrt(c,0);
    end if;
    put_line("The square root x of the random series :"); put(x);
    y := x*x;
    put_line("The square y of the square root x : "); put(y);
    z := y-c;
    put_line("The equation x*x - c :"); put(z);
  end QuadDobl_Random_Test_sqrt;

  procedure Random_Test_root ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given order
  --   and tests the square root computation.

    use Standard_Dense_Series;

    c : constant Series := Random_Series(order);
    n,i : natural32 := 0;
    ans : character;
    x,y,z : Series;
 
  begin
    new_line;
    put("Give the power n in x**n - c : "); get(n);
    put("Give the index i of the root : "); get(i);
    new_line;
    put("Extra output during the computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("A random series c of order "); put(order,1); put_line(" :");
    put(c);
    if ans = 'y'
     then x := Standard_Algebraic_Series.Root(c,n,i,true);
     else x := Standard_Algebraic_Series.Root(c,n,i);
    end if;
    put("The root x of index "); put(i,1);
    put_line(" of the random series :"); put(x);
    y := x**n;
    put_line("The n-th power y of the root x : "); put(y);
    z := y-c;
    put_line("The equation x**n - c :"); put(z);
  end Random_Test_root;

  procedure Test_Conjugate ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given order
  --   and makes the product with its conjugate.

    use Standard_Dense_Series;

    s : constant Series := Random_Series(order);
    c : constant Series := Conjugate(s);
    p,r,n,q,rq : Series;

  begin
    new_line;
    put("A random series of order "); put(order,1); put_line(" :");
    put(s);
    put_line("Its conjugate : ");
    put(c);
    put_line("Conjugate(s)*s : ");
    put(c*s);
    put_line("s*Conjugate(s) : ");
    put(s*c);
    p := c*s;
    r := Standard_Algebraic_Series.sqrt(p,0);
    put_line("The square root r of Conjugate(s)*s :");
    put(r);
    n := s/r;
    put_line("The normalized series s is s/r :");
    put(n);
    q := Conjugate(n)*n;
    rq := Standard_Algebraic_Series.sqrt(q,0);
    put_line("The norm of the normalized series :");
    put(rq);
  end Test_Conjugate;

  procedure Test_Norm ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given order
  --   and computes its norm.

    use Standard_Dense_Series;
    use Standard_Dense_Series_Norms;

    s : constant Series := Random_Series(order);
    nrm : constant Series := Norm(s);
    ns : constant Series := Normalize(s);
    nrm2 : constant Series := Norm(ns);
  
  begin
    new_line;
    put("A random series of order "); put(order,1); put_line(" :");
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
  end Test_Norm;

  procedure Test_Division ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the division on random series of the given order.

    use Standard_Dense_Series;

    a,b,c : Series;

  begin
    new_line;
    put("Give "); put(order+1,1);
    put_line(" complex numbers for the first series : "); 
    a.order := order;
    for i in 0..a.order loop
      get(a.cff(i));
    end loop;
    put_line("The first series : "); put(a);
    new_line;
    put("Give "); put(order+1,1);
    put_line(" complex numbers for the second series : "); 
    b.order := order;
    for i in 0..b.order loop
      get(b.cff(i));
    end loop;
    put_line("The first series : "); put(b);
    c := a/b;
    new_line;
    put_line("The result of the division "); put(c);
  end Test_Division;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the order of the series.

    order : integer32 := 0;
    ans,prc : character;

  begin
    new_line;
    put_line("MENU with testing operations :");
    put_line("  0. test the computation of 1/(1-t)");
    put_line("  1. square root of a random series");
    put_line("  2. p-th root of a random series");
    put_line("  3. test complex conjugate of a series");
    put_line("  4. test the norm of a series");
    put_line("  5. test division operation");
    put("Type 0, 1, 2, 3, 4, or 5 to make your choice : ");
    Ask_Alternative(ans,"012345");
    new_line;
    put("Give the order of the series : "); get(order);
    new_line;
    put_line("MENU for the precision of the coefficients :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(prc,"012");
    new_line;
    case ans is
      when '0' =>
        case prc is
          when '0' => Standard_Test_Creation(order);
          when '1' => DoblDobl_Test_Creation(order);
          when '2' => QuadDobl_Test_Creation(order);
          when others => null;
        end case;
      when '1' =>
        case prc is
          when '0' => Standard_Random_Test_sqrt(order);
          when '1' => DoblDobl_Random_Test_sqrt(order);
          when '2' => QuadDobl_Random_Test_sqrt(order);
          when others => null;
        end case;
      when '2' => Random_Test_root(order);
      when '3' => Test_Conjugate(order);
      when '4' => Test_Norm(order);
      when '5' => Test_Division(order);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_series;
