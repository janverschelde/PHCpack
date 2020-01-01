with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
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
with Standard_Dense_Series_Norms;
with DoblDobl_Algebraic_Series;
with DoblDobl_Dense_Series_Norms;
with QuadDobl_Algebraic_Series;
with QuadDobl_Dense_Series_Norms;

procedure ts_series is

-- DESCRIPTION :
--   Tests the operations on truncated power series.

  procedure Standard_Test_Creation ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Verifies that 1/(1-t) = 1 + t + t^2 + ...
  --   for a truncated power series with coefficients
  --   in standard double precision.

    use Standard_Complex_Numbers;
    use Standard_Dense_Series;

    s : constant Series := Create(1,degree);
    t : Series := s;
    x,y,z : Series;

  begin
    put("One as series of degree "); put(degree,1); put_line(" :");
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

  procedure DoblDobl_Test_Creation ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Verifies that 1/(1-t) = 1 + t + t^2 + ...
  --   for a truncated power series with coefficients
  --   in double double precision.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Dense_Series;

    s : constant Series := Create(1,degree);
    t : Series := s;
    x,y,z : Series;
    minone : constant double_double := create(-1.0);

  begin
    put("One as series of degree "); put(degree,1); put_line(" :");
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

  procedure QuadDobl_Test_Creation ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Verifies that 1/(1-t) = 1 + t + t^2 + ...
  --   for a truncated power series with coefficients
  --   in double double precision.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Dense_Series;

    s : constant Series := Create(1,degree);
    t : Series := s;
    x,y,z : Series;
    minone : constant quad_double := create(-1.0);

  begin
    put("One as series of degree "); put(degree,1); put_line(" :");
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

  procedure Standard_Random_Test_sqrt ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and tests the square root computation,
  --   in standard double precision.

    use Standard_Dense_Series;

    c : constant Series := Random_Series(degree);
    ans : character;
    x,y,z : Series;
 
  begin
    put("Extra output during the computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("A random series c of degree "); put(degree,1); put_line(" :");
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

  procedure DoblDobl_Random_Test_sqrt ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and tests the square root computation,
  --   in double double precision.

    use DoblDobl_Dense_Series;

    c : constant Series := Random_Series(degree);
    ans : character;
    x,y,z : Series;
 
  begin
    put("Extra output during the computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("A random series c of degree "); put(degree,1); put_line(" :");
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

  procedure QuadDobl_Random_Test_sqrt ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and tests the square root computation,
  --   in quad double precision.

    use QuadDobl_Dense_Series;

    c : constant Series := Random_Series(degree);
    ans : character;
    x,y,z : Series;
 
  begin
    put("Extra output during the computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("A random series c of degree "); put(degree,1); put_line(" :");
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

  procedure Standard_Random_Test_root ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and tests the square root computation,
  --   in standard double precision.

    use Standard_Dense_Series;

    c : constant Series := Random_Series(degree);
    n,i : natural32 := 0;
    ans : character;
    x,y,z : Series;
 
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
     then x := Standard_Algebraic_Series.Root(c,n,i,true);
     else x := Standard_Algebraic_Series.Root(c,n,i);
    end if;
    put("The root x of index "); put(i,1);
    put_line(" of the random series :"); put(x);
    y := x**n;
    put_line("The n-th power y of the root x : "); put(y);
    z := y-c;
    put_line("The equation x**n - c :"); put(z);
  end Standard_Random_Test_root;

  procedure DoblDobl_Random_Test_root ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and tests the square root computation,
  --   in double double precision.

    use DoblDobl_Dense_Series;

    c : constant Series := Random_Series(degree);
    n,i : natural32 := 0;
    ans : character;
    x,y,z : Series;
 
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
     then x := DoblDobl_Algebraic_Series.Root(c,n,i,true);
     else x := DoblDobl_Algebraic_Series.Root(c,n,i);
    end if;
    put("The root x of index "); put(i,1);
    put_line(" of the random series :"); put(x);
    y := x**n;
    put_line("The n-th power y of the root x : "); put(y);
    z := y-c;
    put_line("The equation x**n - c :"); put(z);
  end DoblDobl_Random_Test_root;

  procedure QuadDobl_Random_Test_root ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and tests the square root computation,
  --   in quad double precision.

    use QuadDobl_Dense_Series;

    c : constant Series := Random_Series(degree);
    n,i : natural32 := 0;
    ans : character;
    x,y,z : Series;
 
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
     then x := QuadDobl_Algebraic_Series.Root(c,n,i,true);
     else x := QuadDobl_Algebraic_Series.Root(c,n,i);
    end if;
    put("The root x of index "); put(i,1);
    put_line(" of the random series :"); put(x);
    y := x**n;
    put_line("The n-th power y of the root x : "); put(y);
    z := y-c;
    put_line("The equation x**n - c :"); put(z);
  end QuadDobl_Random_Test_root;

  procedure Standard_Test_Conjugate ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and makes the product with its conjugate,
  --   in standard double precision.

    use Standard_Dense_Series;

    s : constant Series := Random_Series(degree);
    c : constant Series := Conjugate(s);
    p,r,n,q,rq : Series;

  begin
    put("A random series of degree "); put(degree,1); put_line(" :");
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
  end Standard_Test_Conjugate;

  procedure DoblDobl_Test_Conjugate ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and makes the product with its conjugate,
  --   in double double precision.

    use DoblDobl_Dense_Series;

    s : constant Series := Random_Series(degree);
    c : constant Series := Conjugate(s);
    p,r,n,q,rq : Series;

  begin
    put("A random series of degree "); put(degree,1); put_line(" :");
    put(s);
    put_line("Its conjugate : ");
    put(c);
    put_line("Conjugate(s)*s : ");
    put(c*s);
    put_line("s*Conjugate(s) : ");
    put(s*c);
    p := c*s;
    r := DoblDobl_Algebraic_Series.sqrt(p,0);
    put_line("The square root r of Conjugate(s)*s :");
    put(r);
    n := s/r;
    put_line("The normalized series s is s/r :");
    put(n);
    q := Conjugate(n)*n;
    rq := DoblDobl_Algebraic_Series.sqrt(q,0);
    put_line("The norm of the normalized series :");
    put(rq);
  end DoblDobl_Test_Conjugate;

  procedure QuadDobl_Test_Conjugate ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and makes the product with its conjugate,
  --   in quad double precision.

    use QuadDobl_Dense_Series;

    s : constant Series := Random_Series(degree);
    c : constant Series := Conjugate(s);
    p,r,n,q,rq : Series;

  begin
    put("A random series of degree "); put(degree,1); put_line(" :");
    put(s);
    put_line("Its conjugate : ");
    put(c);
    put_line("Conjugate(s)*s : ");
    put(c*s);
    put_line("s*Conjugate(s) : ");
    put(s*c);
    p := c*s;
    r := QuadDobl_Algebraic_Series.sqrt(p,0);
    put_line("The square root r of Conjugate(s)*s :");
    put(r);
    n := s/r;
    put_line("The normalized series s is s/r :");
    put(n);
    q := Conjugate(n)*n;
    rq := QuadDobl_Algebraic_Series.sqrt(q,0);
    put_line("The norm of the normalized series :");
    put(rq);
  end QuadDobl_Test_Conjugate;

  procedure Standard_Test_Norm ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and computes its norm, in standard double precision.

    use Standard_Dense_Series;
    use Standard_Dense_Series_Norms;

    s : constant Series := Random_Series(degree);
    nrm : constant Series := Norm(s);
    ns : constant Series := Normalize(s);
    nrm2 : constant Series := Norm(ns);
  
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

  procedure DoblDobl_Test_Norm ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and computes its norm, in double double precision.

    use DoblDobl_Dense_Series;
    use DoblDobl_Dense_Series_Norms;

    s : constant Series := Random_Series(degree);
    nrm : constant Series := Norm(s);
    ns : constant Series := Normalize(s);
    nrm2 : constant Series := Norm(ns);
  
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
  end DoblDobl_Test_Norm;

  procedure QuadDobl_Test_Norm ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given degree
  --   and computes its norm, in quad double precision.

    use QuadDobl_Dense_Series;
    use QuadDobl_Dense_Series_Norms;

    s : constant Series := Random_Series(degree);
    nrm : constant Series := Norm(s);
    ns : constant Series := Normalize(s);
    nrm2 : constant Series := Norm(ns);
  
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
  end QuadDobl_Test_Norm;

  procedure Standard_Test_Division ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the division on random series of the given degree,
  --   in standard double precision.

    use Standard_Dense_Series;

    a,b,c : Series;

  begin
    put("Give "); put(degree+1,1);
    put_line(" complex numbers for the first series : "); 
    a.deg := degree;
    for i in 0..a.deg loop
      get(a.cff(i));
    end loop;
    put_line("The first series : "); put(a);
    new_line;
    put("Give "); put(degree+1,1);
    put_line(" complex numbers for the second series : "); 
    b.deg := degree;
    for i in 0..b.deg loop
      get(b.cff(i));
    end loop;
    put_line("The first series : "); put(b);
    c := a/b;
    new_line;
    put_line("The result of the division "); put(c);
  end Standard_Test_Division;

  procedure DoblDobl_Test_Division ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the division on random series of the given degree,
  --   in double double precision.

    use DoblDobl_Dense_Series;

    a,b,c : Series;

  begin
    put("Give "); put(degree+1,1);
    put_line(" complex numbers for the first series : "); 
    a.deg := degree;
    for i in 0..a.deg loop
      get(a.cff(i));
    end loop;
    put_line("The first series : "); put(a);
    new_line;
    put("Give "); put(degree+1,1);
    put_line(" complex numbers for the second series : "); 
    b.deg := degree;
    for i in 0..b.deg loop
      get(b.cff(i));
    end loop;
    put_line("The first series : "); put(b);
    c := a/b;
    new_line;
    put_line("The result of the division "); put(c);
  end DoblDobl_Test_Division;

  procedure QuadDobl_Test_Division ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the division on random series of the given degree,
  --   in quad double precision.

    use QuadDobl_Dense_Series;

    a,b,c : Series;

  begin
    put("Give "); put(degree+1,1);
    put_line(" complex numbers for the first series : "); 
    a.deg := degree;
    for i in 0..a.deg loop
      get(a.cff(i));
    end loop;
    put_line("The first series : "); put(a);
    new_line;
    put("Give "); put(degree+1,1);
    put_line(" complex numbers for the second series : "); 
    b.deg := degree;
    for i in 0..b.deg loop
      get(b.cff(i));
    end loop;
    put_line("The first series : "); put(b);
    c := a/b;
    new_line;
    put_line("The result of the division "); put(c);
  end QuadDobl_Test_Division;

  procedure Standard_Test_Arithmetic ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Does a basic test on the arithmetic in standard double precision,
  --   on random series of the given degree.

    use Standard_Dense_Series;

    a,b,c : Series;
    ans : character;

  begin
    a := Random_Series(degree);
    put_line("The first random series A :"); put(a);
    b := Random_Series(degree);
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

  procedure DoblDobl_Test_Arithmetic ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Does a basic test on the arithmetic in double double precision,
  --   on random series of the given degree.

    use DoblDobl_Dense_Series;

    a,b,c : Series;
    ans : character;

  begin
    a := Random_Series(degree);
    put_line("The first random series A :"); put(a);
    b := Random_Series(degree);
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
  end DoblDobl_Test_Arithmetic;

  procedure QuadDobl_Test_Arithmetic ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Does a basic test on the arithmetic in quad double precision,
  --   on random series of the given degree.

    use QuadDobl_Dense_Series;

    a,b,c : Series;
    ans : character;

  begin
    a := Random_Series(degree);
    put_line("The first random series A :"); put(a);
    b := Random_Series(degree);
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
  end QuadDobl_Test_Arithmetic;

  procedure Standard_Test_Shift ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Does a basic test on shifting the series parameter
  --   on random series of the given degree.

    use Standard_Complex_Numbers;
    use Standard_Dense_Series;

    s : constant Series := Random_Series(degree);
    rc : double_float := 0.0;
    shifteds : Series;
    cc,y,z : Complex_Number;

  begin
    put_line("Shifting the series parameter ...");
    put_line("on a random series s :"); put(s);
    put("Give a real constant for the shift : "); get(rc);
    shifteds := Shift(s,rc);
    y := Eval(s,rc);
    z := Eval(shifteds,0.0);
    put("s(shift constant) : "); put(y); new_line;
    put("shifted series(0) : "); put(z); new_line;
    new_line;
    put_line("Testing with a complex shift ...");
    put("Give a complex number for the shift : "); get(cc);
    shifteds := Shift(s,cc);
    y := Eval(s,cc);
    z := Eval(shifteds,0.0);
    put("s(shift constant) : "); put(y); new_line;
    put("shifted series(0) : "); put(z); new_line;
  end Standard_Test_Shift;

  procedure DoblDobl_Test_Shift ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Does a basic test on shifting the series parameter
  --   on random series of the given degree.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Dense_Series;

    s : constant Series := Random_Series(degree);
    zero : constant double_double := create(0.0);
    rc : double_double := zero;
    shifteds : Series;
    cc,y,z : Complex_Number;

  begin
    put_line("Shifting the series parameter ...");
    put_line("on a random series s :"); put(s);
    put("Give a real constant for the shift : "); get(rc);
    shifteds := Shift(s,rc);
    y := Eval(s,rc);
    z := Eval(shifteds,zero);
    put("s(shift constant) : "); put(y); new_line;
    put("shifted series(0) : "); put(z); new_line;
    new_line;
    put_line("Testing with a complex shift ...");
    put("Give a complex number for the shift : "); get(cc);
    shifteds := Shift(s,cc);
    y := Eval(s,cc);
    z := Eval(shifteds,zero);
    put("s(shift constant) : "); put(y); new_line;
    put("shifted series(0) : "); put(z); new_line;
  end DoblDobl_Test_Shift;

  procedure QuadDobl_Test_Shift ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Does a basic test on shifting the series parameter
  --   on random series of the given degree.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Dense_Series;

    s : constant Series := Random_Series(degree);
    zero : constant quad_double := create(0.0);
    rc : quad_double := zero;
    shifteds : Series;
    cc,y,z : Complex_Number;

  begin
    put_line("Shifting the series parameter ...");
    put_line("on a random series s :"); put(s);
    put("Give a real constant for the shift : "); get(rc);
    shifteds := Shift(s,rc);
    y := Eval(s,rc);
    z := Eval(shifteds,zero);
    put("s(shift constant) : "); put(y); new_line;
    put("shifted series(0) : "); put(z); new_line;
    new_line;
    put_line("Testing with a complex shift ...");
    put("Give a complex number for the shift : "); get(cc);
    shifteds := Shift(s,cc);
    y := Eval(s,cc);
    z := Eval(shifteds,zero);
    put("s(shift constant) : "); put(y); new_line;
    put("shifted series(0) : "); put(z); new_line;
  end QuadDobl_Test_Shift;
 
  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree of the series.

    degree : integer32 := 0;
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
    put_line("  6. test arithmetic");
    put_line("  7. test shift of series parameter");
    put("Type 0, 1, 2, 3, 4, 5, 6, or 7 to make your choice : ");
    Ask_Alternative(ans,"01234567");
    new_line;
    put("Give the degree of the series : "); get(degree);
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
          when '0' => Standard_Test_Creation(degree);
          when '1' => DoblDobl_Test_Creation(degree);
          when '2' => QuadDobl_Test_Creation(degree);
          when others => null;
        end case;
      when '1' =>
        case prc is
          when '0' => Standard_Random_Test_sqrt(degree);
          when '1' => DoblDobl_Random_Test_sqrt(degree);
          when '2' => QuadDobl_Random_Test_sqrt(degree);
          when others => null;
        end case;
      when '2' =>
        case prc is
          when '0' => Standard_Random_Test_root(degree);
          when '1' => DoblDobl_Random_Test_root(degree);
          when '2' => QuadDobl_Random_Test_root(degree);
          when others => null;
        end case;
      when '3' =>
        case prc is
          when '0' => Standard_Test_Conjugate(degree);
          when '1' => DoblDobl_Test_Conjugate(degree);
          when '2' => QuadDobl_Test_Conjugate(degree);
          when others => null;
        end case;
      when '4' =>
        case prc is 
          when '0' => Standard_Test_Norm(degree);
          when '1' => DoblDobl_Test_Norm(degree);
          when '2' => QuadDobl_Test_Norm(degree);
          when others => null;
        end case;
      when '5' =>
        case prc is
          when '0' => Standard_Test_Division(degree);
          when '1' => DoblDobl_Test_Division(degree);
          when '2' => QuadDobl_Test_Division(degree);
          when others => null;
        end case;
      when '6' =>
        case prc is
          when '0' => Standard_Test_Arithmetic(degree);
          when '1' => DoblDobl_Test_Arithmetic(degree);
          when '2' => QuadDobl_Test_Arithmetic(degree);
          when others => null;
        end case;
      when '7' => 
        case prc is
          when '0' => Standard_Test_Shift(degree);
          when '1' => DoblDobl_Test_Shift(degree);
          when '2' => QuadDobl_Test_Shift(degree);
          when others => null;
        end case;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_series;
