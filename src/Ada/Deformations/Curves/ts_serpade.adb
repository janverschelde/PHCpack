with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Rational_Approximations;

procedure ts_serpade is

-- DESCRIPTION :
--   Interactive development of the rational approximation of a function,
--   given the power series at the origin.

  function log1plusx ( dim : integer32 )
                     return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the first dim+1 coefficients of the series of log(1+x)
  --   as a vector of range 0..dim.

    res : Standard_Complex_Vectors.Vector(0..dim);

  begin
    res(0) := Create(0.0);
    for k in 1..dim loop
      res(k) := Create(1.0/double_float(k));
      if k mod 2 = 0
       then Min(res(k));
      end if;
    end loop;
    return res;
  end log1plusx;

  function invfactorial ( n : integer32 ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns 1/n! where n! is the factorial,
  --   stored as a complex number.

    res : Complex_Number;
    fac : integer32 := 1;
    invfac : double_float;

  begin
    for k in 2..n loop
      fac := fac*k;
    end loop;
    invfac := 1.0/double_float(fac);
    res := Create(invfac);
    return res;
  end invfactorial;

  function expseries ( dim : integer32 )
                     return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of exp(x) at x = 0.

    res : Standard_Complex_Vectors.Vector(0..dim);
    plus : boolean := true;

  begin
    for k in res'range loop
      res(k) := invfactorial(k);
    end loop;
    return res;
  end expseries;

  function sinseries ( dim : integer32 )
                     return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of sin(x) at x = 0.

    res : Standard_Complex_Vectors.Vector(0..dim);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 0 then
        res(k) := Create(0.0);
      else
        res(k) := invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end sinseries;

  function cosseries ( dim : integer32 )
                     return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of cos(x) at x = 0.

    res : Standard_Complex_Vectors.Vector(0..dim);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 1 then
        res(k) := Create(0.0);
      else
        res(k) := invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end cosseries;

  procedure Standard_Logarithm_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the natural logarithm of 1 + x.

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := log1plusx(dim);
    num : Standard_Complex_Vectors.Vector(0..numdeg);
    den : Standard_Complex_Vectors.Vector(0..dendeg);
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.LN(1.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den);
    put_line("The coefficients of the numerator :"); put_line(num);
    put_line("The coefficients of the denominator :"); put_line(den);
    eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
    put("The value at 1.1      :"); put(eva); new_line;
    put("The value of log(1.1) :"); put(chkpnt); new_line;
  end Standard_Logarithm_Test;

  procedure Standard_Sin_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the series of sin(x) at x = 0.

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := sinseries(dim);
    num : Standard_Complex_Vectors.Vector(0..numdeg);
    den : Standard_Complex_Vectors.Vector(0..dendeg);
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.SIN(0.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den);
    put_line("The coefficients of the numerator :"); put_line(num);
    put_line("The coefficients of the denominator :"); put_line(den);
    eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
    put("The value at 0.1      :"); put(eva); new_line;
    put("The value of sin(0.1) :"); put(chkpnt); new_line;
  end Standard_Sin_Test;

  procedure Standard_Exp_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the series of exp(x) at x = 0.

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := expseries(dim);
    num : Standard_Complex_Vectors.Vector(0..numdeg);
    den : Standard_Complex_Vectors.Vector(0..dendeg);
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.EXP(0.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den);
    put_line("The coefficients of the numerator :"); put_line(num);
    put_line("The coefficients of the denominator :"); put_line(den);
    eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
    put("The value at 0.1      :"); put(eva); new_line;
    put("The value of sin(0.1) :"); put(chkpnt); new_line;
  end Standard_Exp_Test;

  procedure Standard_Cos_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the series of cos(x) at x = 0.

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := cosseries(dim);
    num : Standard_Complex_Vectors.Vector(0..numdeg);
    den : Standard_Complex_Vectors.Vector(0..dendeg);
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.COS(0.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den);
    put_line("The coefficients of the numerator :"); put_line(num);
    put_line("The coefficients of the denominator :"); put_line(den);
    eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
    put("The value at 0.1      :"); put(eva); new_line;
    put("The value of cos(0.1) :"); put(chkpnt); new_line;
  end Standard_Cos_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension
  --   and the coefficient vector of the series.

    degnum : integer32 := 0; -- degree numerator
    degden : integer32 := 0; -- degree denominator
    dim : integer32 := 0;    -- degnum + degden
    ans : character;

  begin
    new_line;
    put_line("MENU to test rational approximations for a function ...");
    put_line("  1. series for log(1+x)");
    put_line("  2. series for exp(x)");
    put_line("  3. series for sin(x)");
    put_line("  4. series for cos(x)");
    put("Type 1, 2, 3, or 4 to choose the test : ");
    Ask_Alternative(ans,"1234");
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    dim := degnum + degden;
    put("The dimension : "); put(dim,1); new_line;
    case ans is 
      when '1' => Standard_Logarithm_Test(degnum,degden);
      when '2' => Standard_Exp_Test(degnum,degden);
      when '3' => Standard_Sin_Test(degnum,degden);
      when '4' => Standard_Cos_Test(degnum,degden);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serpade;
