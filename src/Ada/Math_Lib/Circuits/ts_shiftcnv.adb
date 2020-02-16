with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Series;
with Standard_Complex_Series_io;         use Standard_Complex_Series_io;
with Standard_Complex_Series_Functions;
with Standard_Complex_Random_Series;
with Binomial_Coefficients;

procedure ts_shiftcnv is

-- DESCRIPTION :
--   Tests the development of the procedures to shift convolution circuits.

  procedure Shift ( c,wrk : in out Standard_Complex_Vectors.Vector;
                    t : in double_float ) is

  -- DESCRIPTION :
  --   The coefficients in c of the series x(s) are shifted to 
  --   correspond to the coefficients of the series x(s-t).

  -- ON ENTRY :
  --   c        coefficients of a power series;
  --   wrk      work space vector of same range as c;
  --   t        value of the shift.

  -- ON RETURN :
  --   c        shifted coefficients;
  --   wrk      work space equals c.

    use Standard_Complex_Numbers;
    use Binomial_Coefficients;

    bcf : double_float;
    sgn : integer32;

  begin
    for i in 0..c'last loop
      wrk(i) := Create(0.0);
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
        bcf := double_float(sgn*binomial(i,j));
        bcf := bcf*(t**(natural(i-j)));
        wrk(j) := wrk(j) + c(i)*bcf;
        sgn := -sgn;
      end loop;
    end loop;
    c := wrk;
  end Shift;

  procedure Shift ( c,wrk : in out Standard_Complex_Vectors.Vector;
                    t : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   The coefficients in c of the series x(s) are shifted to 
  --   correspond to the coefficients of the series x(s-t).

  -- ON ENTRY :
  --   c        coefficients of a power series;
  --   wrk      work space vector of same range as c;
  --   t        value of the shift.

  -- ON RETURN :
  --   c        shifted coefficients;
  --   wrk      work space equals c.

    use Standard_Complex_Numbers;
    use Binomial_Coefficients;

    bcf : double_float;
    rcf : Complex_Number;
    sgn : integer32;

  begin
    for i in 0..c'last loop
      wrk(i) := Create(0.0);
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
        bcf := double_float(sgn*binomial(i,j));
        rcf := bcf*(t**(natural(i-j)));
        wrk(j) := wrk(j) + c(i)*rcf;
        sgn := -sgn;
      end loop;
    end loop;
    c := wrk;
  end Shift;

  procedure Standard_Test ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the shifting of the series parameter on a random series
  --   of the given degree, in double precision.

    use Standard_Complex_Numbers;
    use Standard_Complex_Series;
    use Standard_Complex_Series_Functions;

    s : constant Series(degree)
      := Standard_Complex_Random_Series.Random_Series(degree);
    rc : double_float := 0.0;
    shifteds : Series(degree);
    cc,y,z : Complex_Number;
    cff,wrk : Standard_Complex_Vectors.Vector(0..degree);

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
    put_line("The shifted coefficients :"); put(shifteds);
    cff := s.cff;
    Shift(cff,wrk,rc);
    put_line("The shifted coefficient vector :"); put_line(cff);
    new_line;
    put_line("Testing with a complex shift ...");
    put("Give a complex number for the shift : "); get(cc);
    shifteds := Shift(s,cc);
    y := Eval(s,-cc);
    z := Eval(shifteds,0.0);
    put("s(-shift constant) : "); put(y); new_line;
    put(" shifted series(0) : "); put(z); new_line;
    put_line("The shifted coefficients :"); put(shifteds);
    cff := s.cff;
    Shift(cff,wrk,cc);
    put_line("The shifted coefficient vector :"); put_line(cff);
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a degree and then launches the test.

    deg : integer32 := 0;

  begin
    new_line;
    put_line("Testing the shifting of coefficients of circuits ...");
    new_line;
    put("Give the degree : "); get(deg);
    Standard_Test(deg);
  end Main;


begin
  Main;
end ts_shiftcnv;
