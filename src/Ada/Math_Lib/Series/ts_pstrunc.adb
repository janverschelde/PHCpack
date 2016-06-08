with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Truncated_Series;         use Standard_Truncated_Series;

procedure ts_pstrunc is

-- DESCRIPTION :
--   Test on truncated power series.

  procedure One_Test_Inverse ( order : in integer32 ) is

  -- DESCRIPTION :
  --   We have that 1/(1-t) = 1 + t + t^2 + t^3 + ...
  --   This procedure tests this identity.

    c : Standard_Complex_Vectors.Vector(0..order);
    x,r : Standard_Complex_Vectors.Vector(0..order);

  begin
    c(0) := Create(1.0);
    c(1) := Create(-1.0);
    c(2..order) := Standard_Complex_Vectors.Vector'(2..order => Create(0.0));
    put("Coefficients of a series of order ");
    put(order,1); put_line(" :");
    put_line(c);
    x := Inverse(c);
    put_line("Coefficients of its inverse : ");
    put_line(x);
    r := c*x;
    put_line("Coefficients of product with its inverse : ");
    put_line(r);
    r := x*c;
    put_line("Check on the commutativity of the product : ");
    put_line(r);
  end One_Test_Inverse;

  procedure Random_Test_Inverse ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a truncated series with random coefficients
  --   of the given order and computes its inverse.
  --   The test succeeds if the product of the first random series
  --   with its computed inverse results in one.

    c : Standard_Complex_Vectors.Vector(0..order) := Random_Vector(0,order);
    x,r : Standard_Complex_Vectors.Vector(0..order);

  begin
    put("Coefficients of a series of order ");
    put(order,1); put_line(" :");
    put_line(c);
    x := Inverse(c);
    put_line("Coefficients of its inverse : ");
    put_line(x);
    r := c*x;
    put_line("Coefficients of product with its inverse : ");
    put_line(r);
    r := x*c;
    put_line("Check on the commutativity of the product : ");
    put_line(r);
  end Random_Test_Inverse;

  procedure Random_Test_Sqrt ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a truncated series with random coefficients
  --   of the given order and computes its square root.

    c : Standard_Complex_Vectors.Vector(0..order) := Random_Vector(0,order);
    x,y,dyc : Standard_Complex_Vectors.Vector(0..order);

  begin
    put("Coefficients of a series of order ");
    put(order,1); put_line(" :");
    put_line(c);
    x := Sqrt(c);
    put_line("Coefficients of its square root : ");
    put_line(x);
    y := x*x;
    put_line("The square of the square root : ");
    put_line(y);
    dyc := Standard_Complex_Vectors."-"(y,c);
    put_line("The coefficients of the difference : ");
    put_line(dyc);
  end Random_Test_Sqrt;

  procedure Convergence_Test ( order : in integer32 ) is

  -- DESCRIPTION :
  --   Consider x^2 - c = 0 for a complex number c,
  --   and a series solution x of the given order.
  --   Evaluation of the series x at t = 0, 0.1, 0.01, etc
  --   should give accurate solutions for c with error as O(t^order).

    c : Standard_Complex_Vectors.Vector(0..order) := Random_Vector(0,order);
    x,y,dyc : Standard_Complex_Vectors.Vector(0..order);
    zero : constant Complex_Number := Create(0.0);
    t : constant double_float := 0.1;
    xt,ct,yt : Complex_Number;
    dytct : double_float;

  begin
    c(2..order) := Standard_Complex_Vectors.Vector'(2..order => zero);
    put("Coefficients of a series of order ");
    put(order,1); put_line(" :");
    put_line(c);
    x := Sqrt(c);
    put_line("Coefficients of its square root : ");
    put_line(x);
    y := x*x;
    put_line("The square of the square root : ");
    put_line(y);
    dyc := Standard_Complex_Vectors."-"(y,c);
    put_line("The coefficients of the difference : ");
    put_line(dyc);
    put("Evaluation at t = "); put(t); put_line(" :");
    ct := Eval(c,t);
    xt := Eval(x,t);
    yt := xt*xt;
    put("      x(t) : "); put(xt); new_line;
    put("      c(t) : "); put(ct); new_line;
    put(" x(t)*x(t) : "); put(yt); new_line;
    dytct := AbsVal(yt - ct);
    put("difference : "); put(dytct,2,3,3); new_line;
  end Convergence_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the order of the truncated series
  --   and runs some test on series with random coefficients.

    order : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the truncation order : "); get(order);
    new_line;
    put_line("MENU to test truncated power series :");
    put_line("  0. test inverse of 1/(1-t)");
    put_line("  1. test inverse of a random series");
    put_line("  2. compute square root of a series");
    put_line("  3. perform a convergence test");
    put("Type 0, 1, 2, or 3 to choose : ");
    Ask_Alternative(ans,"0123");
    new_line;
    case ans is
      when '0' => One_Test_Inverse(order);
      when '1' => Random_Test_Inverse(order);
      when '2' => Random_Test_Sqrt(order);
      when '3' => Convergence_Test(order);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_pstrunc;
