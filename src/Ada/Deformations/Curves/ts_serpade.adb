with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Mathematical_Functions;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;

procedure ts_serpade is

-- DESCRIPTION :
--   Interactive development of the rational approximation of a function,
--   given the power series at the origin.

  function standard_log_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the first dim+1 coefficients of the series of log(1+x)
  --   as a vector of range 0..dim as a vector of complex numbers in
  --   standard double precision.

    use Standard_Complex_Numbers;

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
  end standard_log_series;

  function dobldobl_log_series
             ( dim : integer32 ) return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the first dim+1 coefficients of the series of log(1+x)
  --   as a vector of range 0..dim as a vector of complex numbers in
  --   double double precision.

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Vectors.Vector(0..dim);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    val : double_double;

  begin
    res(0) := Create(zero);
    for k in 1..dim loop
      val := create(k);
      val := one/val;
      res(k) := Create(val);
      if k mod 2 = 0
       then Min(res(k));
      end if;
    end loop;
    return res;
  end dobldobl_log_series;

  function standard_invfactorial
             ( n : integer32 )
             return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns 1/n! where n! is the factorial,
  --   stored as a complex number in standard double precision.

    use Standard_Complex_Numbers;

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
  end standard_invfactorial;

  function dobldobl_invfactorial
             ( n : integer32 )
             return DoblDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns 1/n! where n! is the factorial,
  --   stored as a complex number in double double precision.

    use DoblDobl_Complex_Numbers;

    res : Complex_Number;
    fac : integer32 := 1;
    dd_fac : double_double;
    invfac : double_double;

  begin
    for k in 2..n loop
      fac := fac*k;
    end loop;
    dd_fac := create(fac);
    invfac := 1.0/dd_fac;
    res := Create(invfac);
    return res;
  end dobldobl_invfactorial;

  function standard_exp_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of exp(x) at x = 0, as a vector
  --   of complex numbers in standard double precision.

    res : Standard_Complex_Vectors.Vector(0..dim);
    plus : boolean := true;

  begin
    for k in res'range loop
      res(k) := standard_invfactorial(k);
    end loop;
    return res;
  end standard_exp_series;

  function dobldobl_exp_series
             ( dim : integer32 ) return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of exp(x) at x = 0, as a vector
  --   of complex numbers in double double precision.

    res : DoblDobl_Complex_Vectors.Vector(0..dim);
    plus : boolean := true;

  begin
    for k in res'range loop
      res(k) := dobldobl_invfactorial(k);
    end loop;
    return res;
  end dobldobl_exp_series;

  function standard_sin_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of sin(x) at x = 0, as a vector
  --   of complex numbers in standard double precision.

    use Standard_Complex_Numbers;

    res : Standard_Complex_Vectors.Vector(0..dim);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 0 then
        res(k) := Create(0.0);
      else
        res(k) := standard_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end standard_sin_series;

  function dobldobl_sin_series
             ( dim : integer32 ) return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of sin(x) at x = 0, as a vector
  --   of complex numbers in double double precision.

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Vectors.Vector(0..dim);
    zero : constant double_double := create(0.0);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 0 then
        res(k) := Create(zero);
      else
        res(k) := dobldobl_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end dobldobl_sin_series;

  function standard_cos_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of cos(x) at x = 0, as a vector
  --   of complex numbers in standard double precision.

    use Standard_Complex_Numbers;

    res : Standard_Complex_Vectors.Vector(0..dim);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 1 then
        res(k) := Create(0.0);
      else
        res(k) := standard_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end standard_cos_series;

  function dobldobl_cos_series
             ( dim : integer32 ) return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of cos(x) at x = 0, as a vector
  --   of complex numbers in double double precision.

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Vectors.Vector(0..dim);
    zero : constant double_double := create(0.0);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 1 then
        res(k) := Create(zero);
      else
        res(k) := dobldobl_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end dobldobl_cos_series;

  procedure Standard_log_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the natural logarithm of 1 + x.

    use Standard_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := standard_log_series(dim);
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
  end Standard_log_Test;

  procedure DoblDobl_log_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in double double arithmetic
  --   on the natural logarithm of 1 + x.

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
        := DoblDobl_log_series(dim);
    num : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    den : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    nbr : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(nbr);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.LN(1.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den);
    put_line("The coefficients of the numerator :"); put_line(num);
    put_line("The coefficients of the denominator :"); put_line(den);
    eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
    put("The value at 1.1      :"); put(eva); new_line;
    put("The value of log(1.1) :"); put(chkpnt); new_line;
  end DoblDobl_log_Test;

  procedure Standard_sin_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the series of sin(x) at x = 0.

    use Standard_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := standard_sin_series(dim);
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

  procedure DoblDobl_sin_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in double double arithmetic
  --   on the series of sin(x) at x = 0.

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
        := dobldobl_sin_series(dim);
    num : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    den : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    arg : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant double_double
           := DoblDobl_Mathematical_Functions.SIN(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den);
    put_line("The coefficients of the numerator :"); put_line(num);
    put_line("The coefficients of the denominator :"); put_line(den);
    eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
    put("The value at 0.1      :"); put(eva); new_line;
    put("The value of sin(0.1) :"); put(chkpnt); new_line;
  end DoblDobl_Sin_Test;

  procedure Standard_exp_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the series of exp(x) at x = 0.

    use Standard_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := standard_exp_series(dim);
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
  end Standard_exp_Test;

  procedure DoblDobl_exp_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in double double arithmetic
  --   on the series of exp(x) at x = 0.

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
        := dobldobl_exp_series(dim);
    num : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    den : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    arg : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.EXP(0.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den);
    put_line("The coefficients of the numerator :"); put_line(num);
    put_line("The coefficients of the denominator :"); put_line(den);
    eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
    put("The value at 0.1      :"); put(eva); new_line;
    put("The value of sin(0.1) :"); put(chkpnt); new_line;
  end DoblDobl_exp_Test;

  procedure Standard_cos_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the series of cos(x) at x = 0.

    use Standard_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := standard_cos_series(dim);
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
  end Standard_cos_Test;

  procedure DoblDobl_cos_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in double double arithmetic
  --   on the series of cos(x) at x = 0.

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
        := dobldobl_cos_series(dim);
    num : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    den : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    arg : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant double_double
           := DoblDobl_Mathematical_Functions.COS(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den);
    put_line("The coefficients of the numerator :"); put_line(num);
    put_line("The coefficients of the denominator :"); put_line(den);
    eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
    put("The value at 0.1      :"); put(eva); new_line;
    put("The value of cos(0.1) :"); put(chkpnt); new_line;
  end DoblDobl_cos_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension
  --   and the coefficient vector of the series.

    degnum : integer32 := 0; -- degree numerator
    degden : integer32 := 0; -- degree denominator
    dim : integer32 := 0;    -- degnum + degden
    ans,prc : character;

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
    new_line;
    put_line("MENU for the precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put("Type 0 or 1 to select the precision : ");
    Ask_Alternative(prc,"01");
    case prc is
      when '0' =>
        case ans is 
          when '1' => Standard_log_Test(degnum,degden);
          when '2' => Standard_exp_Test(degnum,degden);
          when '3' => Standard_sin_Test(degnum,degden);
          when '4' => Standard_cos_Test(degnum,degden);
          when others => null;
        end case;
      when '1' =>
        case ans is 
          when '1' => DoblDobl_log_Test(degnum,degden);
          when '2' => DoblDobl_exp_Test(degnum,degden);
          when '3' => DoblDobl_sin_Test(degnum,degden);
          when '4' => DoblDobl_cos_Test(degnum,degden);
          when others => null;
        end case;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serpade;
