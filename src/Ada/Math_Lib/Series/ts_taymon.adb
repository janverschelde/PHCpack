with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Binomial_Coefficients;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
-- with Standard_Random_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;

procedure ts_taymon is

-- DESCRIPTION :
--   Tests on Taylor series developments of monomials
--   with positive real exponents.

  function Double_Taylor_Coefficients 
             ( deg : integer32; alpha, point : double_float )
             return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector of coefficients of the Taylor series
  --   in the shifted basis.

  -- ON ENTRY :
  --   deg     truncation degree of the Taylor series;
  --   alpha   positive real power of the monomial;
  --   point   base point of the development.

    res : Standard_Floating_Vectors.Vector(0..deg) := (0..deg => 0.0);
    cff : double_float := alpha;

  begin
    res(0) := point**alpha;
    res(1) := alpha*(point**(alpha-1.0));
    for k in 2..deg loop
      cff := cff*(alpha-1.0)/double_float(k);
      res(k) := cff*(point**(alpha-double_float(k)));
    end loop;
    return res;
  end Double_Taylor_Coefficients;

  function Double_Taylor_Value
             ( cff : Standard_Floating_Vectors.Vector;
               arg, point : double_float ) return double_float is

  -- DESCRIPTION :
  --   Returns the value of the Taylor expansion about the point,
  --   given the coefficients of the series in the shifted basis.

  -- ON ENTRY :
  --   deg     truncation degree of the Taylor series;
  --   arg     argument where to evaluate the series;
  --   point   base point of the development.

    res : double_float := cff(0);
    inc : double_float := 1.0;
    argminpoint : constant double_float := arg - point;

  begin
    for k in 1..cff'last loop
      inc := inc*argminpoint;
      res := res + cff(k)*inc;
    end loop;
    return res;
  end Double_Taylor_Value;

  function Double_Taylor_Value
             ( cff : Standard_Floating_Vectors.Vector;
               arg : double_float ) return double_float is

  -- DESCRIPTION :
  --   Returns the value of the Taylor expansion about the point,
  --   given the coefficients of the series in the monomial basis.

  -- ON ENTRY :
  --   deg     truncation degree of the Taylor series;
  --   arg     argument where to evaluate the series.

    res : double_float := cff(0);
    inc : double_float := 1.0;

  begin
    for k in 1..cff'last loop
      inc := inc*arg;
      res := res + cff(k)*inc;
    end loop;
    return res;
  end Double_Taylor_Value;

  function Double_Taylor_Expansion
             ( cff : Standard_Floating_Vectors.Vector;
               point : double_float )
             return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the Taylor expansion about the point,
  --   expressed in the standard monomial basis.

  -- ON ENTRY :
  --   deg     truncation degree of the Taylor series;
  --   point   base point of the development.

    res : Standard_Floating_Vectors.Vector(cff'range) := (cff'range => 0.0);
    bincff : double_float;
    expcff : natural;

  begin
    for k in cff'range loop
      for ell in 0..k loop
        expcff := natural(k - ell);
        bincff := Binomial_Coefficients.binomial(k,ell)*(point**expcff);
        if expcff mod 2 = 0
         then res(ell) := res(ell) + cff(k)*bincff;
         else res(ell) := res(ell) - cff(k)*bincff;
        end if;
      end loop;
    end loop;
    return res;
  end Double_Taylor_Expansion;

  procedure Double_Test
              ( deg : in integer32; alpha, point : in double_float ) is

  -- DESCRIPTION :
  --   Runs a test on the double Taylor coefficients.
  --
  -- ON ENTRY :
  --   deg      truncation degree of the Taylor series;
  --   alpha    positive real power of the monomial;
  --   point    base point of the development.

    cff : constant Standard_Floating_Vectors.Vector(0..deg)
        := Double_Taylor_Coefficients(deg,alpha,point);
    tcf : constant Standard_Floating_Vectors.Vector(0..deg)
        := Double_Taylor_Expansion(cff,point);
    apt : double_float := point - 1.0e-2; -- Standard_Random_Numbers.Random;
    mval,tval1,tval2,err : double_float;

  begin
    new_line;
    put_line("The coefficients in the shifted basis :");
    put_line(cff);
    put_line("The coefficients in the monomial basis :");
    put_line(tcf);
    if apt < 0.0
     then apt := -apt;
    end if;
    put("The argument to evaluate the series :"); put(apt); new_line;
    tval1 := Double_Taylor_Value(cff,apt,point);
    put("The value in the shifted basis  :"); put(tval1); new_line;
    tval2 := Double_Taylor_Value(tcf,apt);
    put("The value in the monomial basis :"); put(tval2); new_line;
    mval := apt**alpha;
    put("The value of the monomial       :"); put(mval); new_line;
    err := abs(mval - tval1);
    put("-> the error on the values :"); put(err,3); new_line;
  end Double_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the input parameters.

    deg : integer32 := 0;
    alpha,point : double_float := 0.0;

  begin
    new_line;
    put("Give the truncation degree : "); get(deg);
    put("Give the positive real power : "); get(alpha);
    put("Give the positive real point : "); get(point);
    new_line;
    put("-> the truncation degree : "); put(deg,1); new_line;
    put("-> power of the monomial :"); put(alpha); new_line;
    put("-> point of development  :"); put(point); new_line;
    Double_Test(deg,alpha,point);
  end Main;

begin
  Main;
end ts_taymon;
