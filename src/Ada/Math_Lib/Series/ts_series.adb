with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Dense_Series;              use Standard_Dense_Series;
with Standard_Dense_Series_io;           use Standard_Dense_Series_io;

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

  function Random_Series ( order : integer32 ) return Series is

  -- DESCRIPTION :
  --   Returns a series of the given order, with random coefficient,
  --   on the unit circle on the complex plane.

    cff : Standard_Complex_Vectors.Vector(0..order)
        := Random_Vector(0,order);

  begin
    return Create(cff);
  end Random_Series;

  function sqrt ( c : Series ) return Series is

  -- ALGORITHM : apply Newton's method to x^2 - c = 0,
  --   compute dx = (x^2 - c)/(2*x), followed by x = x - dx.
  --   Stops when x.order > c.order.

    x,wrk,dx : Series;
    fac : constant Complex_Number := Create(0.5);

  begin
    x.order := 1;
    x.cff(0) := Standard_Complex_Numbers_Polar.Root(c.cff(0),2,1);
    x.cff(1) := Create(0.0);
    loop
      for k in 1..3+x.order loop
        wrk := x*x - c;
        dx := fac*wrk/x;
        x := x - dx;
      end loop;
      exit when x.order >= c.order;
      x := Create(x,2*x.order);
    end loop;
    return x;
  end sqrt;

  procedure Random_Test_sqrt ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given order
  --   and tests the square root computation.

    c : constant Series := Random_Series(order);
    x,y,z : Series;
 
  begin
    put("A random series c of order "); put(order,1); put_line(" :");
    put(c);
    x := Sqrt(c);
    put_line("The square root x of the random series :"); put(x);
    y := x*x;
    put_line("The square y of the square root x : "); put(y);
    z := y-c;
    put_line("The equation x*x - c :"); put(z);
  end Random_Test_sqrt;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the order of the series.

    order : integer32 := 0;

  begin
    new_line;
    put("Give the order of the series : "); get(order);
    -- Test_Creation(order);
    Random_Test_sqrt(order);
  end Main;

begin
  Main;
end ts_series;
