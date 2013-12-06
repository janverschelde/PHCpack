with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Multprec_Natural_Coefficients;
with Multprec_Natural64_Coefficients;

package body Multprec_Random_Numbers is

  function Random ( size : natural32 )
                  return Multprec_Natural_Numbers.Natural_Number is

    use Multprec_Natural_Coefficients,Multprec_Natural_Numbers;

    res : Natural_Number;
    the_base : constant natural32 := Multprec_Natural_Coefficients.Base;
    arn : Array_of_Naturals(0..size);

  begin
    for i in arn'range loop
      arn(i) := natural32(Standard_Random_Numbers.Random
                            (0,integer32(the_base)-1));
    end loop;
    res := Create(arn);
    return res;
  end Random;

  function Random ( size : natural32 )
                  return Multprec_Natural64_Numbers.Natural_Number is

    use Multprec_Natural64_Coefficients,Multprec_Natural64_Numbers;

    res : Natural_Number;
    the_base : constant natural64 := Multprec_Natural64_Coefficients.Base;
    arn : Array_of_Naturals(0..size);
    r : integer64;

  begin
    for i in arn'range loop
      r := Standard_Random_Numbers.Random(0,integer64(the_base-1));
      arn(i) := natural64(r);
    end loop;
    res := Create(arn);
    return res;
  end Random;

  function Random ( size : natural32 )
                  return Multprec_Integer_Numbers.Integer_Number is

    use Multprec_Natural_Numbers,Multprec_Integer_Numbers;

    res : Integer_Number;
    sig : constant double_float := Standard_Random_Numbers.Random;
    rnd : Natural_Number := Random(size);

  begin
    res := Create(rnd);
    if sig < 0.0
     then Min(res);
    end if;
    Clear(rnd);
    return res;
  end Random;

  function Random ( size : natural32 )
                  return Multprec_Integer64_Numbers.Integer_Number is

    use Multprec_Natural64_Numbers,Multprec_Integer64_Numbers;

    res : Integer_Number;
    sig : constant double_float := Standard_Random_Numbers.Random;
    rnd : Natural_Number := Random(size);

  begin
    res := Create(rnd);
    if sig < 0.0
     then Min(res);
    end if;
    Clear(rnd);
    return res;
  end Random;

  function Random ( size : natural32 )
                  return Multprec_Floating_Numbers.Floating_Number is

    use Multprec_Natural_Numbers,Multprec_Integer_Numbers;
    use Multprec_Floating_Numbers;

    res : Floating_Number;
    naturfrac : Natural_Number := Random(size);
    frac,expo : Integer_Number;
    szfrac : integer32;

  begin
    frac := Create(naturfrac);
    Clear(naturfrac);
    szfrac := integer32(Decimal_Places(frac));
    expo := Create(-szfrac);
    res := Create(frac,expo);
    Clear(frac); Clear(expo);
    return res;
  end Random;

  function Random ( size : natural32 )
                  return Multprec_Floating64_Numbers.Floating_Number is

    use Multprec_Natural64_Numbers,Multprec_Integer64_Numbers;
    use Multprec_Floating64_Numbers;

    res : Floating_Number;
    naturfrac : Natural_Number := Random(size);
    frac,expo : Integer_Number;

  begin
    frac := Create(naturfrac);
    Clear(naturfrac);
    expo := Create64(-integer64(Decimal_Places(frac)));
    res := Create(frac,expo);
    Clear(frac); Clear(expo);
    return res;
  end Random;

  function Random ( size : natural32; lower,upper : integer32 )
                  return Multprec_Floating_Numbers.Floating_Number is

    use Multprec_Floating_Numbers;

    res : Floating_Number := Random(size);
    exp : constant integer32 := Standard_Random_Numbers.Random(lower,upper);

  begin
    if exp > 0 then
      for i in 1..exp loop
        Mul(res,10.0);
      end loop;
    elsif exp < 0 then
      for i in 1..(-exp) loop
        Div(res,10.0);
      end loop;
    end if;
    return res;
  end Random;

  function Random ( size : natural32; lower,upper : integer64 )
                  return Multprec_Floating64_Numbers.Floating_Number is

    use Multprec_Floating64_Numbers;

    res : Floating_Number := Random(size);
    exp : constant integer64 := Standard_Random_Numbers.Random(lower,upper);

  begin
    if exp > 0 then
      for i in 1..exp loop
        Mul(res,10.0);
      end loop;
    elsif exp < 0 then
      for i in 1..(-exp) loop
        Div(res,10.0);
      end loop;
    end if;
    return res;
  end Random;

  function Random ( size : natural32 ) return Complex_Number is

    res : Complex_Number;
    re : Multprec_Floating_Numbers.Floating_Number := Random(size);
    im : Multprec_Floating_Numbers.Floating_Number := Random(size);

  begin
    res := Create(re,im);
    Multprec_Floating_Numbers.Clear(re);
    Multprec_Floating_Numbers.Clear(im);
    return res;
  end Random;

  function Random ( size : natural32; lower,upper : integer32 )
                  return Multprec_Complex_Numbers.Complex_Number is

    use Multprec_Floating_Numbers;

    res : Multprec_Complex_Numbers.Complex_Number;
    resre,resim : Multprec_Floating_Numbers.Floating_Number;

  begin
    resre := Random(size,lower,upper);
    resim := Random(size,lower,upper);
    res := Create(resre,resim);
    Clear(resre); Clear(resim);
    return res;
  end Random;

end Multprec_Random_Numbers;
