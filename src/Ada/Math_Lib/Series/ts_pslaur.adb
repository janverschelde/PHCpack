with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the degree of truncation, leading exponents,
  --   and tests some basic arithmetical operations.

     d,ale,ble,prodle,invble,quotle : integer32 := 0;

  begin
    new_line;
    put("Give the truncation degree : "); get(d);
    put("Give the leading exponent of a : "); get(ale);
    put("Give the leading exponent of b : "); get(ble);
    declare
      a : constant Standard_Complex_Vectors.Vector(0..d)
        := Standard_Random_Vectors.Random_Vector(0,d);
      b : constant Standard_Complex_Vectors.Vector(0..d)
        := Standard_Random_Vectors.Random_Vector(0,d);
      prod : Standard_Complex_Vectors.Vector(0..d);
      invb : Standard_Complex_Vectors.Vector(0..d);
      quot : Standard_Complex_Vectors.Vector(0..d);
    begin
      new_line;
      put_line("A random series a :"); Write(ale,a);
      new_line;
      put_line("A random series b :"); Write(ble,b);
      Multiply(d,ale,ble,a,b,prodle,prod);
      new_line;
      put_line("The product of a with b :"); Write(prodle,prod);
      Divide(d,prodle,ble,prod,b,quotle,quot,invb);
      new_line;
      put_line("The quotient of a*b with b :"); Write(quotle,quot);
      invble := -ble;
      Multiply(d,ble,invble,b,invb,prodle,prod);
      put_line("The product of b with 1/b :"); Write(prodle,prod);
    end;
  end Main;

begin
  Main;
end ts_pslaur;
