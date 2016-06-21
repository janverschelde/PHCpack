with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with Standard_Complex_Vectors;
with Standard_Dense_Series;              use Standard_Dense_Series;
with Standard_Dense_Series_io;           use Standard_Dense_Series_io;
with Standard_Random_Series;             use Standard_Random_Series;
with Standard_Algebraic_Series;

procedure ts_series is

-- DESCRIPTION :
--   Tests the operations on truncated power series.

  procedure Test_Creation ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Verifies that 1/(1-t) = 1 + t + t^2 + ...

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
  end Test_Creation;

  procedure Random_Test_sqrt ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given order
  --   and tests the square root computation.

    c : constant Series := Random_Series(order);
    ans : character;
    x,y,z : Series;
 
  begin
    new_line;
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
  end Random_Test_sqrt;

  procedure Random_Test_root ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given order
  --   and tests the square root computation.

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

    s : constant Series := Random_Series(order);
    c : constant Series := Conjugate(s);
    p,r,n : Series;

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
    put_line("The normed series s is s/r :");
    put(n);
  end Test_Conjugate;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the order of the series.

    order : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the order of the series : "); get(order);
    new_line;
    put_line("MENU with testing operations :");
    put_line("  0. test the computation of 1/(1-t)");
    put_line("  1. square root of a random series");
    put_line("  2. p-th root of a random series");
    put_line("  3. test complex conjugate of a series");
    put("Type 0, 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' => Test_Creation(order);
      when '1' => Random_Test_sqrt(order);
      when '2' => Random_Test_root(order);
      when '3' => Test_Conjugate(order);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_series;
