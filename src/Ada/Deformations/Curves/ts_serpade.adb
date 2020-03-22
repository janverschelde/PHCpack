with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
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
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Mathematical_Functions;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_Vectors_io;
with Standard_CSeries_Vector_Functions;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors_io;
with DoblDobl_CSeries_Vector_Functions;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors_io;
with QuadDobl_CSeries_Vector_Functions;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Homotopy_Series_Readers;
with Standard_Pade_Approximants;
with Standard_Pade_Approximants_io;
with DoblDobl_Pade_Approximants;
with DoblDobl_Pade_Approximants_io;
with QuadDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants_io;
with Homotopy_Pade_Approximants;

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

  function quaddobl_log_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the first dim+1 coefficients of the series of log(1+x)
  --   as a vector of range 0..dim as a vector of complex numbers in
  --   quad double precision.

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Vectors.Vector(0..dim);
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    val : quad_double;

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
  end quaddobl_log_series;

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

  function quaddobl_invfactorial
             ( n : integer32 )
             return QuadDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns 1/n! where n! is the factorial,
  --   stored as a complex number in quad double precision.

    use QuadDobl_Complex_Numbers;

    res : Complex_Number;
    fac : integer32 := 1;
    dd_fac : quad_double;
    invfac : quad_double;

  begin
    for k in 2..n loop
      fac := fac*k;
    end loop;
    dd_fac := create(fac);
    invfac := 1.0/dd_fac;
    res := Create(invfac);
    return res;
  end quaddobl_invfactorial;

  function standard_exp_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of exp(x) at x = 0, as a vector
  --   of complex numbers in standard double precision.

    res : Standard_Complex_Vectors.Vector(0..dim);

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

  begin
    for k in res'range loop
      res(k) := dobldobl_invfactorial(k);
    end loop;
    return res;
  end dobldobl_exp_series;

  function quaddobl_exp_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of exp(x) at x = 0, as a vector
  --   of complex numbers in quad double precision.

    res : QuadDobl_Complex_Vectors.Vector(0..dim);

  begin
    for k in res'range loop
      res(k) := quaddobl_invfactorial(k);
    end loop;
    return res;
  end quaddobl_exp_series;

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

  function quaddobl_sin_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of sin(x) at x = 0, as a vector
  --   of complex numbers in quad double precision.

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Vectors.Vector(0..dim);
    zero : constant quad_double := create(0.0);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 0 then
        res(k) := Create(zero);
      else
        res(k) := quaddobl_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end quaddobl_sin_series;

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

  function quaddobl_cos_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of cos(x) at x = 0, as a vector
  --   of complex numbers in quad double precision.

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Vectors.Vector(0..dim);
    zero : constant quad_double := create(0.0);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 1 then
        res(k) := Create(zero);
      else
        res(k) := quaddobl_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end quaddobl_cos_series;

  procedure Standard_line_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction of the Pade approximant in standard
  --   double precision for the line x = t.

    use Standard_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : Standard_Complex_Vectors.Vector(0..dim);
    num : Standard_Complex_Vectors.Vector(0..numdeg);
    den : Standard_Complex_Vectors.Vector(0..dendeg);
    info : integer32;

  begin
    cff := (0..dim => Create(0.0));
    cff(1) := Create(1.0);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
    end if;
  end Standard_line_Test;

  procedure DoblDobl_line_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction of the Pade approximant in double
  --   double precision for the line x = t.

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : DoblDobl_Complex_Vectors.Vector(0..dim);
    num : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    den : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    info : integer32;

  begin
    cff := (0..dim => Create(zero));
    cff(1) := Create(one);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
    end if;
  end DoblDobl_line_Test;

  procedure QuadDobl_line_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction of the Pade approximant in quad
  --   double precision for the line x = t.

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : QuadDobl_Complex_Vectors.Vector(0..dim);
    num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    info : integer32;

  begin
    cff := (0..dim => Create(zero));
    cff(1) := Create(one);
    QuadDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
    end if;
  end QuadDobl_line_Test;

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
    info : integer32;
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.LN(1.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 1.1      :"); put(eva); new_line;
      put("The value of log(1.1) :"); put(chkpnt); new_line;
    end if;
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
    info : integer32;
    nbr : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(nbr);
    eva : Complex_Number;
    arg : constant double_double := create(1.1);
    chkpnt : constant double_double := Double_Double_Numbers.log(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 1.1      : "); put(eva); new_line;
      put("The value of log(1.1) : "); put(chkpnt); new_line;
    end if;
  end DoblDobl_log_Test;

  procedure QuadDobl_log_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in quad double arithmetic
  --   on the natural logarithm of 1 + x.

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
        := QuadDobl_log_series(dim);
    num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    nbr : constant quad_double := create(0.1);
    pnt : constant Complex_Number := Create(nbr);
    eva : Complex_Number;
    arg : constant quad_double := create(1.1);
    chkpnt : constant quad_double := Quad_Double_Numbers.log(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    QuadDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := QuadDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 1.1      : "); put(eva); new_line;
      put("The value of log(1.1) : "); put(chkpnt); new_line;
    end if;
  end QuadDobl_log_Test;

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
    info : integer32;
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.SIN(0.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      :"); put(eva); new_line;
      put("The value of sin(0.1) :"); put(chkpnt); new_line;
    end if;
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
    info : integer32;
    arg : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant double_double
           := DoblDobl_Mathematical_Functions.SIN(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of sin(0.1) : "); put(chkpnt); new_line;
    end if;
  end DoblDobl_sin_Test;

  procedure QuadDobl_sin_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in quad double arithmetic
  --   on the series of sin(x) at x = 0.

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
        := quaddobl_sin_series(dim);
    num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    arg : constant quad_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant quad_double
           := QuadDobl_Mathematical_Functions.SIN(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    QuadDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := QuadDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of sin(0.1) : "); put(chkpnt); new_line;
    end if;
  end QuadDobl_Sin_Test;

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
    info : integer32;
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.EXP(0.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      :"); put(eva); new_line;
      put("The value of exp(0.1) :"); put(chkpnt); new_line;
    end if;
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
    info : integer32;
    arg : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant double_double := Double_Double_Numbers.exp(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of exp(0.1) : "); put(chkpnt); new_line;
    end if;
  end DoblDobl_exp_Test;

  procedure QuadDobl_exp_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in quad double arithmetic
  --   on the series of exp(x) at x = 0.

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
        := quaddobl_exp_series(dim);
    num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    arg : constant quad_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant quad_double := Quad_Double_Numbers.exp(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    QuadDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := QuadDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of exp(0.1) : "); put(chkpnt); new_line;
    end if;
  end QuadDobl_exp_Test;

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
    info : integer32;
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.COS(0.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      :"); put(eva); new_line;
      put("The value of cos(0.1) :"); put(chkpnt); new_line;
    end if;
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
    info : integer32;
    arg : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant double_double
           := DoblDobl_Mathematical_Functions.COS(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of cos(0.1) : "); put(chkpnt); new_line;
    end if;
  end DoblDobl_cos_Test;

  procedure QuadDobl_cos_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in quad double arithmetic
  --   on the series of cos(x) at x = 0.

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
        := quaddobl_cos_series(dim);
    num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    arg : constant quad_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant quad_double
           := QuadDobl_Mathematical_Functions.COS(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    QuadDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := QuadDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of cos(0.1) : "); put(chkpnt); new_line;
    end if;
  end QuadDobl_cos_Test;

  procedure Standard_Pade_Approximation
              ( nbequ,nbsteps : in integer32;
                srv : in Standard_Complex_Series_Vectors.Vector;
                pv : in Standard_Pade_Approximants.Pade_Vector ) is

  -- DESCRIPTION :
  --   The Pade approximant pv and the series srv are evaluated in 
  --   as many points as the value of nbsteps.

    arg : double_float;
    pnt : Standard_Complex_Numbers.Complex_Number;
    value : Standard_Complex_Vectors.Vector(1..nbequ);
    approx : Standard_Complex_Vectors.Vector(pv'range);
    valsrv : Standard_Complex_Vectors.Vector(srv'range);

  begin
    for k in 1..nbsteps loop
      arg := double_float(k)*0.04;
      pnt := Standard_Complex_Numbers.Create(arg);
      valsrv := Standard_CSeries_Vector_Functions.Eval(srv,pnt);
      put_line("The value of the series approximation :");
      put_line(valsrv);
      approx := Standard_Pade_Approximants.Eval(pv,pnt);
      put_line("The value of the rational approximation :");
      put_line(approx);
      value := Standard_Homotopy.Eval(valsrv,pnt);
      put("Series approximation evaluated at "); put(arg,3); put_line(" :");
      put_line(value);
      value := Standard_Homotopy.Eval(approx,pnt);
      put("Pade approximation evaluated at "); put(arg,3); put_line(" :");
      put_line(value);
    end loop;
  end Standard_Pade_Approximation;

  procedure DoblDobl_Pade_Approximation
              ( nbequ,nbsteps : in integer32;
                srv : in DoblDobl_Complex_Series_Vectors.Vector;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector ) is

  -- DESCRIPTION :
  --   The Pade approximant pv and the series srv are evaluated in
  --   as many points as the value of nbsteps.

    arg : double_double;
    pnt : DoblDobl_Complex_Numbers.Complex_Number;
    value : DoblDobl_Complex_Vectors.Vector(1..nbequ);
    approx : DoblDobl_Complex_Vectors.Vector(pv'range);
    valsrv : DoblDobl_Complex_Vectors.Vector(srv'range);

  begin
    for k in 1..nbsteps loop
      arg := Double_Double_Numbers.create(double_float(k)*0.04);
      pnt := DoblDobl_Complex_Numbers.Create(arg);
      valsrv := DoblDobl_CSeries_Vector_Functions.Eval(srv,pnt);
      put_line("The value of the series approximation :");
      put_line(valsrv);
      approx := DoblDobl_Pade_Approximants.Eval(pv,pnt);
      put_line("The value of the rational approximation :");
      put_line(approx);
      value := DoblDobl_Homotopy.Eval(valsrv,pnt);
      put("Series approximation evaluated at "); put(arg,3); put_line(" :");
      put_line(value);
      value := DoblDobl_Homotopy.Eval(approx,pnt);
      put("Pade approximant evaluated at "); put(arg,3); put_line(" :");
      put_line(value);
    end loop;
  end DoblDobl_Pade_Approximation;

  procedure QuadDobl_Pade_Approximation
              ( nbequ,nbsteps : in integer32;
                srv : in QuadDobl_Complex_Series_Vectors.Vector;
                pv : in QuadDobl_Pade_Approximants.Pade_Vector ) is

  -- DESCRIPTION :
  --   The Pade approximant pv and the series srv are evaluated in 
  --   as many points as the value of nbsteps.

    dd_arg : double_double;
    qd_arg : quad_double;
    pnt : QuadDobl_Complex_Numbers.Complex_Number;
    value : QuadDobl_Complex_Vectors.Vector(1..nbequ);
    approx : QuadDobl_Complex_Vectors.Vector(pv'range);
    valsrv : QuadDobl_Complex_Vectors.Vector(srv'range);

  begin
    for k in 1..nbsteps loop
      dd_arg := Double_Double_Numbers.create(double_float(k)*0.04);
      qd_arg := Quad_Double_Numbers.create(dd_arg);
      pnt := QuadDobl_Complex_Numbers.Create(qd_arg);
      valsrv := QuadDobl_CSeries_Vector_Functions.Eval(srv,pnt);
      put_line("The value of the series approximation :");
      put_line(valsrv);
      approx := QuadDobl_Pade_Approximants.Eval(pv,pnt);
      put_line("The value of the rational approximation :");
      put_line(approx);
      value := QuadDobl_Homotopy.Eval(valsrv,pnt);
      put("Series approximation evaluated at "); put(qd_arg,3); put_line(" :");
      put_line(value);
      value := QuadDobl_Homotopy.Eval(approx,pnt);
      put("Pade approximation evaluated at "); put(qd_arg,3); put_line(" :");
      put_line(value);
    end loop;
  end QuadDobl_Pade_Approximation;

  procedure Standard_Test_Homotopy is

  -- DESCRIPTION :
  --   Stores the test homotopy, starting at x^2 - 1, and with
  --   target system 3*x^2 - 3/2 in the Standard_Homotopy data.
  --   The test homotopy (1-t)*(x^2 - 1) + t*(3*x^2 - 3/2) = 0
  --   expands into x^2 - 1 - t*x^2 + t + t*3*x^2 - 3/2*t = 0
  --   which leads to (1+2*t)*x^2 = 1 + 1/2*t and thus defines
  --   the function x(t) = ((1 + 1/2*t)/(1 + 2*t))^(1/2).

    use Standard_Complex_Polynomials;

    start,target : Standard_Complex_Poly_Systems.Poly_Sys(1..1);
    startpol,targetpol : Poly;
    trm : Term;
    tpow : constant natural32 := 1;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Complex_Numbers.Create(1.0);

  begin
    trm.cf := Standard_Complex_Numbers.Create(1.0);
    trm.dg := new Standard_Natural_Vectors.Vector'(1..1 => 0);
    trm.dg(1) := 2;
    startpol := Create(trm);
    trm.dg(1) := 0;
    Sub(startpol,trm);
    start(1) := startpol;
    trm.dg(1) := 2;
    trm.cf := Standard_Complex_Numbers.Create(3.0);
    targetpol := Create(trm);
    trm.dg(1) := 0;
    trm.cf := Standard_Complex_Numbers.Create(1.5);
    Sub(targetpol,trm);
    target(1) := targetpol;
    Standard_Homotopy.Create(target,start,tpow,gamma);
  end Standard_Test_Homotopy;

  procedure DoblDobl_Test_Homotopy is

  -- DESCRIPTION :
  --   Stores the test homotopy, starting at x^2 - 1, and with
  --   target system 3*x^2 - 3/2 in the DoblDobl_Homotopy data.
  --   The test homotopy (1-t)*(x^2 - 1) + t*(3*x^2 - 3/2) = 0
  --   expands into x^2 - 1 - t*x^2 + t + t*3*x^2 - 3/2*t = 0
  --   which leads to (1+2*t)*x^2 = 1 + 1/2*t and thus defines
  --   the function x(t) = ((1 + 1/2*t)/(1 + 2*t))^(1/2).

    use DoblDobl_Complex_Polynomials;

    start,target : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..1);
    startpol,targetpol : Poly;
    trm : Term;
    tpow : constant natural32 := 1;
    nbr : double_double := create(1.0);
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Complex_Numbers.Create(nbr);

  begin
    trm.cf := DoblDobl_Complex_Numbers.Create(nbr);
    trm.dg := new Standard_Natural_Vectors.Vector'(1..1 => 0);
    trm.dg(1) := 2;
    startpol := Create(trm);
    trm.dg(1) := 0;
    Sub(startpol,trm);
    start(1) := startpol;
    trm.dg(1) := 2;
    nbr := create(3.0);
    trm.cf := DoblDobl_Complex_Numbers.Create(nbr);
    targetpol := Create(trm);
    trm.dg(1) := 0;
    nbr := create(1.5);
    trm.cf := DoblDobl_Complex_Numbers.Create(nbr);
    Sub(targetpol,trm);
    target(1) := targetpol;
    DoblDobl_Homotopy.Create(target,start,tpow,gamma);
  end DoblDobl_Test_Homotopy;

  procedure QuadDobl_Test_Homotopy is

  -- DESCRIPTION :
  --   Stores the test homotopy, starting at x^2 - 1, and with
  --   target system 3*x^2 - 3/2 in the QuadDobl_Homotopy data.
  --   The test homotopy (1-t)*(x^2 - 1) + t*(3*x^2 - 3/2) = 0
  --   expands into x^2 - 1 - t*x^2 + t + t*3*x^2 - 3/2*t = 0
  --   which leads to (1+2*t)*x^2 = 1 + 1/2*t and thus defines
  --   the function x(t) = ((1 + 1/2*t)/(1 + 2*t))^(1/2).

    use QuadDobl_Complex_Polynomials;

    start,target : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..1);
    startpol,targetpol : Poly;
    trm : Term;
    tpow : constant natural32 := 1;
    nbr : quad_double := create(1.0);
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Complex_Numbers.Create(nbr);

  begin
    trm.cf := QuadDobl_Complex_Numbers.Create(nbr);
    trm.dg := new Standard_Natural_Vectors.Vector'(1..1 => 0);
    trm.dg(1) := 2;
    startpol := Create(trm);
    trm.dg(1) := 0;
    Sub(startpol,trm);
    start(1) := startpol;
    trm.dg(1) := 2;
    nbr := create(3.0);
    trm.cf := QuadDobl_Complex_Numbers.Create(nbr);
    targetpol := Create(trm);
    trm.dg(1) := 0;
    nbr := create(1.5);
    trm.cf := QuadDobl_Complex_Numbers.Create(nbr);
    Sub(targetpol,trm);
    target(1) := targetpol;
    QuadDobl_Homotopy.Create(target,start,tpow,gamma);
  end QuadDobl_Test_Homotopy;

  procedure Standard_Test_Start_Solutions
              ( sols : out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Returns in sols the two start solutions +1 and -1
  --   for the test homotopy, in standard double precision.

    sol : Standard_Complex_Solutions.Solution(1);
    two,one : Standard_Complex_Solutions.Link_to_Solution;

  begin
    sol.t := Standard_Complex_Numbers.Create(0.0);
    sol.m := 1;
    sol.v(1) := Standard_Complex_Numbers.Create(-1.0);
    sol.err := 0.0; sol.rco := 1.0; sol.res := 0.0; 
    two := new Standard_Complex_Solutions.Solution'(sol);
    Standard_Complex_Solutions.Construct(two,sols);
    sol.v(1) := Standard_Complex_Numbers.Create(1.0);
    one := new Standard_Complex_Solutions.Solution'(sol);
    Standard_Complex_Solutions.Construct(one,sols);
  end Standard_Test_Start_Solutions;

  procedure DoblDobl_Test_Start_Solutions
              ( sols : out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Returns in sols the two start solutions +1 and -1
  --   for the test homotopy, in double double precision.

    sol : DoblDobl_Complex_Solutions.Solution(1);
    nbr : double_double;
    two,one : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    nbr := create(0.0);
    sol.t := DoblDobl_Complex_Numbers.Create(nbr);
    sol.m := 1;
    nbr := create(-1.0);
    sol.v(1) := DoblDobl_Complex_Numbers.Create(nbr);
    sol.err := create(0.0);
    sol.rco := create(1.0);
    sol.res := create(0.0); 
    two := new DoblDobl_Complex_Solutions.Solution'(sol);
    DoblDobl_Complex_Solutions.Construct(two,sols);
    nbr := create(1.0);
    sol.v(1) := DoblDobl_Complex_Numbers.Create(nbr);
    one := new DoblDobl_Complex_Solutions.Solution'(sol);
    DoblDobl_Complex_Solutions.Construct(one,sols);
  end DoblDobl_Test_Start_Solutions;

  procedure QuadDobl_Test_Start_Solutions
              ( sols : out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Returns in sols the two start solutions +1 and -1
  --   for the test homotopy, in quad double precision.

    sol : QuadDobl_Complex_Solutions.Solution(1);
    nbr : quad_double;
    two,one : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    nbr := create(0.0);
    sol.t := QuadDobl_Complex_Numbers.Create(nbr);
    sol.m := 1;
    nbr := create(-1.0);
    sol.v(1) := QuadDobl_Complex_Numbers.Create(nbr);
    sol.err := create(0.0);
    sol.rco := create(1.0);
    sol.res := create(0.0); 
    two := new QuadDobl_Complex_Solutions.Solution'(sol);
    QuadDobl_Complex_Solutions.Construct(two,sols);
    nbr := create(1.0);
    sol.v(1) := QuadDobl_Complex_Numbers.Create(nbr);
    one := new QuadDobl_Complex_Solutions.Solution'(sol);
    QuadDobl_Complex_Solutions.Construct(one,sols);
  end QuadDobl_Test_Start_Solutions;

  procedure Standard_Test_Case
              ( sols : out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Defines the homotopy for an example test case,
  --   in standard double precision.

  begin
    Standard_Test_Homotopy;
    Standard_Test_Start_Solutions(sols);
  end Standard_Test_Case;

  procedure DoblDobl_Test_Case
              ( sols : out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Defines the homotopy for an example test case,
  --   in double double precision.

  begin
    DoblDobl_Test_Homotopy;
    DoblDobl_Test_Start_Solutions(sols);
  end DoblDobl_Test_Case;

  procedure QuadDobl_Test_Case
              ( sols : out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Defines the homotopy for an example test case,
  --   in quad double precision.

  begin
    QuadDobl_Test_Homotopy;
    QuadDobl_Test_Start_Solutions(sols);
  end QuadDobl_Test_Case;

  procedure Standard_Pade_Homotopy
              ( numdeg,dendeg,nbeq,nbsteps : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Computes Pade approximations for the solution paths
  --   defined by an artificial parameter homotopy,
  --   in standard double precision.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   nbeq     number of equations;
  --   nbsteps  number of steps;
  --   sols     start solutions for an artificial-parameter homotopy.

    lnk : constant Standard_Complex_Solutions.Link_to_Solution
        := Standard_Complex_Solutions.Head_Of(sols);
    sol : constant Standard_Complex_Solutions.Solution := lnk.all;
    nbt : constant natural32 := natural32(numdeg+dendeg+1);
    nit : constant natural32 := 4*nbt;
    srv : Standard_Complex_Series_Vectors.Vector(sol.v'range);
    eva : Standard_Complex_Series_Vectors.Vector(1..nbeq);
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range);

  begin
    Homotopy_Pade_Approximants.Standard_Pade_Approximant
      (sol.v,nbeq+1,nbeq,numdeg,dendeg,nit,srv,eva,pv);
    put_line("The solution series :");
    Standard_Complex_Series_Vectors_io.put(srv);
    put_line("The evaluated solution series :");
    Standard_Complex_Series_Vectors_io.put(eva);
    Standard_Pade_Approximation(nbeq,nbsteps,srv,pv);
    put_line("The Pade approximant :");
    for i in pv'range loop
      put_line(Standard_Pade_Approximants_io.Write(pv(i)));
    end loop;
    Standard_Pade_Approximants.Clear(pv);
  end Standard_Pade_Homotopy;

  procedure DoblDobl_Pade_Homotopy
              ( numdeg,dendeg,nbeq,nbsteps : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Computes Pade approximations for the solution paths
  --   defined by an artificial parameter homotopy,
  --   in double double precision.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   nbeq     number of equations;
  --   nbsteps  number of steps;
  --   sols     start solutions for an artificial-parameter homotopy.

    lnk : constant DoblDobl_Complex_Solutions.Link_to_Solution
        := DoblDobl_Complex_Solutions.Head_Of(sols);
    sol : constant DoblDobl_Complex_Solutions.Solution := lnk.all;
    nbt : constant natural32 := natural32(numdeg+dendeg+1);
    nit : constant natural32 := 4*nbt;
    srv : DoblDobl_Complex_Series_Vectors.Vector(sol.v'range);
    eva : DoblDobl_Complex_Series_Vectors.Vector(1..nbeq);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range);

  begin
    Homotopy_Pade_Approximants.DoblDobl_Pade_Approximant
      (sol.v,nbeq+1,nbeq,numdeg,dendeg,nit,srv,eva,pv);
    put_line("The solution series :");
    DoblDobl_Complex_Series_Vectors_io.put(srv);
    put_line("The evaluated solution series :");
    DoblDobl_Complex_Series_Vectors_io.put(eva);
    DoblDobl_Pade_Approximation(nbeq,nbsteps,srv,pv);
    put_line("The Pade approximant :");
    for i in pv'range loop
      put_line(DoblDobl_Pade_Approximants_io.Write(pv(i)));
    end loop;
    DoblDobl_Pade_Approximants.Clear(pv);
  end DoblDobl_Pade_Homotopy;

  procedure QuadDobl_Pade_Homotopy
              ( numdeg,dendeg,nbeq,nbsteps : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Computes Pade approximations for the solution paths
  --   defined by an artificial parameter homotopy,
  --   in quad double precision.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   nbeq     number of equations;
  --   nbsteps  number of steps;
  --   sols     start solutions for an artificial-parameter homotopy.

    lnk : constant QuadDobl_Complex_Solutions.Link_to_Solution
        := QuadDobl_Complex_Solutions.Head_Of(sols);
    sol : constant QuadDobl_Complex_Solutions.Solution := lnk.all;
    nbt : constant natural32 := natural32(numdeg+dendeg+1);
    nit : constant natural32 := 4*nbt;
    srv : QuadDobl_Complex_Series_Vectors.Vector(sol.v'range);
    eva : QuadDobl_Complex_Series_Vectors.Vector(1..nbeq);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range);

  begin
    Homotopy_Pade_Approximants.QuadDobl_Pade_Approximant
      (sol.v,nbeq+1,nbeq,numdeg,dendeg,nit,srv,eva,pv);
    put_line("The solution series :");
    QuadDobl_Complex_Series_Vectors_io.put(srv);
    put_line("The evaluated solution series :");
    QuadDobl_Complex_Series_Vectors_io.put(eva);
    QuadDobl_Pade_Approximation(nbeq,nbsteps,srv,pv);
    put_line("The Pade approximant :");
    for i in pv'range loop
      put_line(QuadDobl_Pade_Approximants_io.Write(pv(i)));
    end loop;
    QuadDobl_Pade_Approximants.Clear(pv);
  end QuadDobl_Pade_Homotopy;

  procedure Standard_Homotopy_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts for a target, start system, and start solutions.
  --   Applies Newton's method for a series development of the first
  --   start solution, in standard double precision.

    nbeq : integer32;
    nbsteps : integer32 := 5;
    sols : Standard_Complex_Solutions.Solution_List;
    ans : character;

  begin
    new_line;
    put("Run an example test case ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y' then
      Standard_Test_Case(sols); nbeq := 1; nbsteps := 50;
    else
      put("Random gamma ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Homotopy_Series_Readers.Standard_Reader(nbeq,sols,tpow=>1);
      else
        declare
          gamma : constant Standard_Complex_Numbers.Complex_Number
                := Standard_Complex_Numbers.Create(1.0);
        begin
          Homotopy_Series_Readers.Standard_Reader(nbeq,sols,1,gamma);
        end;
      end if;
      new_line;
      put("Give the number of steps : "); get(nbsteps);
      new_line;
    end if;
    Standard_Pade_Homotopy(numdeg,dendeg,nbeq,nbsteps,sols);
  end Standard_Homotopy_Test;

  procedure DoblDobl_Homotopy_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts for a target, start system, and start solutions.
  --   Applies Newton's method for a series development of the first
  --   start solution, in double double precision.

    nbeq : integer32;
    nbsteps : integer32 := 5;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    ans : character;

  begin
    new_line;
    put("Run an example test case ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y' then
      DoblDobl_Test_Case(sols); nbeq := 1; nbsteps := 50;
    else
      put("Random gamma ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols,tpow=>1);
      else
        declare
          one : constant double_double := create(1.0);
          gamma : constant DoblDobl_Complex_Numbers.Complex_Number
                := DoblDobl_Complex_Numbers.Create(one);
        begin
          Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols,1,gamma);
        end;
      end if;
      new_line;
      put("Give the number of steps : "); get(nbsteps);
      new_line;
    end if;
    DoblDobl_Pade_Homotopy(numdeg,dendeg,nbeq,nbsteps,sols);
  end DoblDobl_Homotopy_Test;

  procedure QuadDobl_Homotopy_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts for a target, start system, and start solutions.
  --   Applies Newton's method for a series development of the first
  --   start solution, in quad double precision.

    nbeq : integer32;
    nbsteps : integer32 := 5;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    ans : character;

  begin
    new_line;
    put("Run an example test case ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y' then
      QuadDobl_Test_Case(sols); nbeq := 1; nbsteps := 50;
    else
      put("Random gamma ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols,tpow=>1);
      else
        declare
          one : constant quad_double := create(1.0);
          gamma : constant QuadDobl_Complex_Numbers.Complex_Number
                := QuadDobl_Complex_Numbers.Create(one);
        begin
          Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols,1,gamma);
        end;
      end if;
      new_line;
      put("Give the number of steps : "); get(nbsteps);
      new_line;
    end if;
    QuadDobl_Pade_Homotopy(numdeg,dendeg,nbeq,nbsteps,sols);
  end QuadDobl_Homotopy_Test;

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
    put_line("  0. simple test on line x = t");
    put_line("  1. series for log(1+x)");
    put_line("  2. series for exp(x)");
    put_line("  3. series for sin(x)");
    put_line("  4. series for cos(x)");
    put_line("  5. artificial parameter homotopy");
    put("Type 0, 1, 2, 3, 4, or 5 to choose the test : ");
    Ask_Alternative(ans,"012345");
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    dim := degnum + degden;
    put("The dimension : "); put(dim,1); new_line;
    prc := Prompt_for_Precision;
    case prc is
      when '0' =>
        case ans is 
          when '0' => Standard_line_Test(degnum,degden);
          when '1' => Standard_log_Test(degnum,degden);
          when '2' => Standard_exp_Test(degnum,degden);
          when '3' => Standard_sin_Test(degnum,degden);
          when '4' => Standard_cos_Test(degnum,degden);
          when '5' => Standard_Homotopy_Test(degnum,degden);
          when others => null;
        end case;
      when '1' =>
        case ans is 
          when '0' => DoblDobl_line_Test(degnum,degden);
          when '1' => DoblDobl_log_Test(degnum,degden);
          when '2' => DoblDobl_exp_Test(degnum,degden);
          when '3' => DoblDobl_sin_Test(degnum,degden);
          when '4' => DoblDobl_cos_Test(degnum,degden);
          when '5' => DoblDobl_Homotopy_Test(degnum,degden);
          when others => null;
        end case;
      when '2' =>
        case ans is 
          when '0' => QuadDobl_line_Test(degnum,degden);
          when '1' => QuadDobl_log_Test(degnum,degden);
          when '2' => QuadDobl_exp_Test(degnum,degden);
          when '3' => QuadDobl_sin_Test(degnum,degden);
          when '4' => QuadDobl_cos_Test(degnum,degden);
          when '5' => QuadDobl_Homotopy_Test(degnum,degden);
          when others => null;
        end case;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serpade;
