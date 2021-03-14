with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;

procedure ts_pslaur is

-- DESCRIPTION :
--   Performs some basic tests on truncated Laurent power series,
--   defined by a leading exponent (which may be negative)
--   and a complex coefficient vector.

  procedure Write ( e : in integer32;
                    c : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the power series with leading exponent e
  --   and coefficients in c.

  begin
    for i in c'range loop
      if i > c'first
       then put(" + (");
       else put("   (");
      end if;
      put(c(i)); put(")*t^"); put(e+i,1); new_line;
    end loop;
  end Write;

  procedure Multiply ( d,xe,ye : in integer32;
                       xc,yc : in Standard_Complex_Vectors.Vector;
                       ze : out integer32;
                       zc : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Multiplies two Laurent series.

  -- REQUIRED :
  --   All coefficient vectors range between 0 and d,
  --   where d is the same constant for all series.

  -- ON ENTRY :
  --   d       only coefficients in the range 0 to d are considered;
  --   xe      leading exponent of the first series x;
  --   ye      leading exponent of the second series y;
  --   xc      coefficient vector of the first series x;
  --   yc      coefficient vector of the second series y.

  -- ON RETURN :
  --   ze      leading exponent of the product of x with y;
  --   zc      coefficient vector of the product of x with y.

  begin
    ze := xe + ye;
    for i in 0..d loop
       zc(i) := xc(0)*yc(i);
       for j in 1..i loop
         zc(i) := zc(i) + xc(j)*yc(i-j);
       end loop;
    end loop;
  end Multiply;

  procedure Inverse ( d,xe : in integer32;
                      xc : in Standard_Complex_Vectors.Vector;
                      ye : out integer32;
                      yc : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the inverse of a Laurent series.
  --   Any nonzero Laurent series has a nonzero leading coefficient.

  -- REQUIRED : xc'first = 0 = yc'first and xc'last = yc'last.

  -- ON ENTRY :
  --   d       truncation degree of the series;
  --   xe      leading exponent of a nonzero series x;
  --   xc      coefficient vector of the series x.

  -- ON RETURN :
  --   ye      leading exponent of the inverse of x;
  --   yc      coefficient vector of the inverse of x.

  begin
    ye := -xe;
    yc(0) := 1.0/xc(0);
    for i in 1..d loop
      yc(i) := -xc(1)*yc(i-1);
      for j in 2..i loop
        yc(i) := yc(i) - xc(j)*yc(i-j);
      end loop;
      yc(i) := yc(i)/xc(0);
    end loop;
  end Inverse;

  procedure Divide ( d,xe,ye : in integer32;
                     xc,yc : in Standard_Complex_Vectors.Vector;
                     ze : out integer32;
                     zc : out Standard_Complex_Vectors.Vector;
                     iyc : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Divides two Laurent series.

  -- REQUIRED :
  --   All coefficient vectors range between 0 and d,
  --   where d is the same constant for all series.
  --   The second series y should be nonzero.

  -- ON ENTRY :
  --   d       only coefficients in the range 0 to d are considered;
  --   xe      leading exponent of the first series x;
  --   ye      leading exponent of the second series y;
  --   xc      coefficient vector of the first series x;
  --   yc      coefficient vector of the second series y;

  -- ON RETURN :
  --   ze      leading exponent of the product of x with y;
  --   zc      coefficient vector of the product of x with y;
  --   iyc     coefficient vector of the inverse of y, as work space.

    iye : integer32;

  begin
    Inverse(d,ye,yc,iye,iyc);
    Multiply(d,xe,iye,xc,iyc,ze,zc);
  end Divide;

  procedure Normalize ( d : in integer32; e : in out integer32;
                        c : in out Standard_Complex_Vectors.Vector;
                        tol : in double_float := 1.0E-15 ) is

  -- DESCRIPTION :
  --   Normalizes the representation of the Laurent series
  --   so that the leading coefficient is larger than tol.

  -- ON ENTRY :
  --   d       only coefficients in the range 0 to d are considered;
  --   e       leading exponent of the series;
  --   c       coefficient vector of the series;
  --   tol     tolerance to decide the leading degree
  --           when the first coefficient becomes too small.

  -- ON RETURN :
  --   e       augmented exponent for each shift,
  --           if all coefficients are zero, then e is zero as well;
  --   c       shifted coefficient vector so the leading coefficient
  --           is nonzero, or zero for a zero series.

    allzero : boolean := true;

  begin
    for i in 0..d loop
      if (AbsVal(c(0)) > tol)
       then allzero := false; exit;
      end if;
      e := e + 1;
      for k in 1..d-i loop -- shift the coefficients
        c(k-1) := c(k);
      end loop;
    end loop;
    if allzero
     then e := 0;
    end if;
  end Normalize;

  procedure Add ( d,xe,ye : in integer32;
                  xc,yc : in Standard_Complex_Vectors.Vector;
                  ze : out integer32;
                  zc : out Standard_Complex_Vectors.Vector;
                  tol : in double_float := 1.0E-15 ) is

  -- DESCRIPTION :
  --   Adds two Laurent series, with a tolerance to determine
  --   the leading degree when the leading coefficient is too small.

  -- REQUIRED :
  --   All coefficient vectors range between 0 and d,
  --   where d is the same constant for all series.

  -- ON ENTRY :
  --   d       only coefficients in the range 0 to d are considered;
  --   xe      leading exponent of the first series x;
  --   ye      leading exponent of the second series y;
  --   xc      coefficient vector of the first series x;
  --   yc      coefficient vector of the second series y;
  --   tol     tolerance to decide the leading degree
  --           when the first coefficient becomes too small.

  -- ON RETURN :
  --   ze      leading exponent of the sum of x and y;
  --   zc      coefficient vector of the sum of x and y.

    gap : integer32; -- gap between the leading coefficients

  begin
    if xe < ye then
      ze := xe;
      gap := abs(ye - xe);
      for i in 0..gap-1 loop
        zc(i) := xc(i);
      end loop;
      for i in gap..d loop
        zc(i) := xc(i) + yc(i-gap);
      end loop;
    elsif xe > ye then
      ze := ye;
      gap := abs(xe - ye);
      for i in 0..gap-1 loop
        zc(i) := yc(i);
      end loop;
      for i in gap..d loop
        zc(i) := yc(i) + xc(i-gap);
      end loop;
    else -- xe = ye
      ze := xe;
      for i in 0..d loop
        zc(i) := xc(i) + yc(i);
      end loop;
      Normalize(d,ze,zc,tol);
    end if;
  end Add;

  procedure Subtract ( d,xe,ye : in integer32;
                       xc,yc : in Standard_Complex_Vectors.Vector;
                       ze : out integer32;
                       zc : out Standard_Complex_Vectors.Vector;
                       tol : in double_float := 1.0E-15 ) is

  -- DESCRIPTION :
  --   Subtracts two Laurent series, with a tolerance to determine
  --   the leading degree when the leading coefficient is too small.

  -- REQUIRED :
  --   All coefficient vectors range between 0 and d,
  --   where d is the same constant for all series.

  -- ON ENTRY :
  --   d       only coefficients in the range 0 to d are considered;
  --   xe      leading exponent of the first series x;
  --   ye      leading exponent of the second series y;
  --   xc      coefficient vector of the first series x;
  --   yc      coefficient vector of the second series y;
  --   tol     tolerance to decide the leading degree
  --           when the first coefficient becomes too small.

  -- ON RETURN :
  --   ze      leading exponent of the difference x - y;
  --   zc      coefficient vector of the difference x - y.

    gap : integer32; -- gap between the leading coefficients

  begin
    if xe < ye then
      ze := xe;
      gap := abs(ye - xe);
      for i in 0..gap-1 loop
        zc(i) := xc(i);
      end loop;
      for i in gap..d loop
        zc(i) := xc(i) - yc(i-gap);
      end loop;
    elsif xe > ye then
      ze := ye;
      gap := abs(xe - ye);
      for i in 0..gap-1 loop
        zc(i) := -yc(i);
      end loop;
      for i in gap..d loop
        zc(i) := xc(i-gap) - yc(i);
      end loop;
    else -- xe = ye
      ze := xe;
      for i in 0..d loop
        zc(i) := xc(i) - yc(i);
      end loop;
      Normalize(d,ze,zc,tol);
    end if;
  end Subtract;

  procedure Test_Multiply_Inverse_Divide ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Generates coefficient vectors for truncated Laurent series
  --   to degree deg, to test multiplication, inverse, and division.

    ale,ble,prodle,invble,quotle : integer32 := 0;
    a : constant Standard_Complex_Vectors.Vector(0..deg)
      := Standard_Random_Vectors.Random_Vector(0,deg);
    b : constant Standard_Complex_Vectors.Vector(0..deg)
      := Standard_Random_Vectors.Random_Vector(0,deg);
    prod : Standard_Complex_Vectors.Vector(0..deg);
    invb : Standard_Complex_Vectors.Vector(0..deg);
    quot : Standard_Complex_Vectors.Vector(0..deg);

  begin
    new_line;
    put("Give the leading exponent of a : "); get(ale);
    put("Give the leading exponent of b : "); get(ble);
    new_line;
    put_line("A random series a :"); Write(ale,a);
    new_line;
    put_line("A random series b :"); Write(ble,b);
    Multiply(deg,ale,ble,a,b,prodle,prod);
    new_line;
    put_line("The product of a with b :"); Write(prodle,prod);
    Divide(deg,prodle,ble,prod,b,quotle,quot,invb);
    new_line;
    put_line("The quotient of a*b with b :"); Write(quotle,quot);
    invble := -ble;
    Multiply(deg,ble,invble,b,invb,prodle,prod);
    put_line("The product of b with 1/b :"); Write(prodle,prod);
  end Test_Multiply_Inverse_Divide;

  procedure Test_Add_and_Subtract ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Generates coefficient vectors for truncated Laurent series
  --   to degree deg, to test addition and subtraction.

    ale,ble,sumle,difle : integer32 := 0;
    a : constant Standard_Complex_Vectors.Vector(0..deg)
      := Standard_Random_Vectors.Random_Vector(0,deg);
    b : constant Standard_Complex_Vectors.Vector(0..deg)
      := Standard_Random_Vectors.Random_Vector(0,deg);
    sum : Standard_Complex_Vectors.Vector(0..deg);
    dif : Standard_Complex_Vectors.Vector(0..deg);

  begin
    new_line;
    put("Give the leading exponent of a : "); get(ale);
    put("Give the leading exponent of b : "); get(ble);
    new_line;
    put_line("A random series a :"); Write(ale,a);
    new_line;
    put_line("A random series b :"); Write(ble,b);
    Add(deg,ale,ble,a,b,sumle,sum);
    new_line;
    put_line("The sum of a and b :"); Write(sumle,sum);
    Subtract(deg,sumle,ble,sum,b,difle,dif);
    new_line;
    put_line("The result of (a + b) - b :"); Write(difle,dif);
    Subtract(deg,ale,ale,a,a,difle,dif);
    new_line;
    put_line("The result of a - a :"); Write(difle,dif);
  end Test_Add_and_Subtract;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the degree of truncation,
  --   and tests some basic arithmetical operations.

     d : integer32 := 0;

  begin
    new_line;
    put("Give the truncation degree : "); get(d);
    Test_Multiply_Inverse_Divide(d);
    Test_Add_and_Subtract(d);
  end Main;

begin
  Main;
end ts_pslaur;
