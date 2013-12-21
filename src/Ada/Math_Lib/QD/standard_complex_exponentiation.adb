with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Complex_Numbers_Polar;
with Multprec_Floating_Constants;
with Double_Double_Constants;
with Quad_Double_Constants;
with Multprec_DoblDobl_Convertors;       use Multprec_DoblDobl_Convertors;
with Multprec_QuadDobl_Convertors;       use Multprec_QuadDobl_Convertors;

-- for testing purposes:
--with text_io,integer_io;           use text_io,integer_io;
--with Standard_Floating_Numbers_io; use Standard_Floating_Numbers_io;
--with Multprec_Integer_Numbers_io;  use Multprec_Integer_Numbers_io;
--with Multprec_Floating_Numbers_io; use Multprec_Floating_Numbers_io;
--with Double_Double_Numbers_io;     use Double_Double_Numbers_io;
--with Quad_Double_Numbers_io;       use Quad_Double_Numbers_io;

package body Standard_Complex_Exponentiation is

  dftwopi : constant double_float := 2.0*Standard_Mathematical_Functions.Pi;
  ddtwopi : constant double_double := Double_Double_Constants.twopi;
  qdtwopi : constant quad_double := Quad_Double_Constants.twopi;

  procedure PosDivModTwoPi
               ( x : in double_float;
                 q : out integer32; r : out double_float ) is

  -- DESCRIPTION :
  --   Computes quotient and remainder modulo 2*Pi for x >= 2*Pi.

  -- NOTE : calculations are executed with double double arithmetic
  --   because otherwise the logarithm of the error on the remainder 
  --   could be about the same as the number of decimal places of x,
  --   as subtraction of numbers of equal magnitude is instable.

    dd_x : constant double_double := create(x);
    fq : constant double_double := dd_x/ddtwopi;
    dd_q,dd_r : double_double;
    iq : integer32;

  begin
    iq := integer32(hi_part(fq));
    if double_float(iq) > hi_part(fq)
     then iq := iq - 1;
    end if;
    dd_q := Double_Double_Numbers.create(iq);
    dd_r := (fq - dd_q)*ddtwopi;
    r := hi_part(dd_r);
    if r < 0.0
     then r := r + dftwopi; iq := iq - 1;
    end if;
    q := iq;
  end PosDivModTwoPi;

  procedure PosDivModTwoPi
               ( x : in double_double;
                 q : out Integer_Number; r : out double_double ) is

  -- DESCRIPTION :
  --   Computes quotient and remainder modulo 2*Pi for x >= 2*Pi.

  -- NOTE : calculations are executed with quad double arithmetic
  --   because otherwise the logarithm of the error on the remainder 
  --   could be about the same as the number of decimal places of x,
  --   as subtraction of numbers of equal magnitude is instable.

    qd_x : constant quad_double := create(x);
    fq : constant quad_double := qd_x/qdtwopi;
    fn : Floating_Number := to_floating_number(fq);
    iq : Integer_Number;
    qd_q,qd_r : quad_double;

  begin
    iq := Truncate_to_Nearest_Integer(fn); Clear(fn);
    qd_q := to_quad_double(iq);
    qd_r := (fq - qd_q)*qdtwopi;
    r := hi_part(qd_r);
    if r < 0.0
     then r := r + ddtwopi; Sub(iq,1);
    end if;
    q := iq;
  end PosDivModTwoPi;

  procedure PosDivModTwoPi
               ( x : in quad_double;
                 q : out Integer_Number; r : out quad_double ) is

  -- DESCRIPTION :
  --   Computes quotient and remainder modulo 2*Pi for x >= 2*Pi.

  -- NOTE : calculations are executed with extended arithmetic
  --   because otherwise the logarithm of the error on the remainder 
  --   could be about the same as the number of decimal places of x,
  --   as subtraction of numbers of equal magnitude is instable.

    fq,mp_twopi,mp_q,mp_r : Floating_Number;
    size : natural32;
    iq : Integer_Number;

  begin
    fq := to_floating_number(x); size := Size_Fraction(fq);
    mp_twopi := to_floating_number(qdtwopi);
    Set_Size(fq,2*size); Set_Size(mp_twopi,2*size);
    Div(fq,mp_twopi);
    iq := Truncate_to_Nearest_Integer(fq);
    mp_q := Create(iq);
    mp_r := fq - mp_q;
    Mul(mp_r,mp_twopi);
    r := to_quad_double(mp_r);
    if r < 0.0
     then r := r + qdtwopi; Sub(iq,1);
    end if;
    q := iq;
    Clear(fq); Clear(mp_twopi);
    Clear(mp_q); Clear(mp_r);
  end PosDivModTwoPi;

  procedure PosDivModTwoPi
              ( x,twopi : in Floating_Number; d : in natural32;
                q : out Integer_Number; r : out Floating_Number ) is

  -- DESCRIPTION :
  --   Returns quotient q and remainder r after division of x by 2*Pi,
  --   calculating with 3*d decimal places, for x > 0.

    fq,iq : Floating_Number;

  begin
    Copy(x,fq); Set_Size(fq,3*d);
    Div(fq,twopi);
    q := Truncate_to_Nearest_Integer(fq);
    iq := Create(q);
    r := fq - iq; Set_Size(r,3*d);
    Mul(r,twopi);
    if r < 0.0
     then Add(r,twopi); Sub(q,1);
    end if;
    Set_Size(r,d);
    Clear(fq); Clear(iq);
  end PosDivModTwoPi;

  procedure DivModTwoPi ( x : in double_float;
                          q : out integer32; r : out double_float ) is
  begin
    if x > -dftwopi and x < dftwopi then
      q := 0; r := x;
    elsif x < 0.0 then
      PosDivModTwoPi(-x,q,r);
      q := -q; r := -r;
    else
      PosDivModTwoPi(x,q,r);
    end if;
  end DivModTwoPi;

  procedure DivModTwoPi ( x : in double_double;
                          q : out Integer_Number; r : out double_double ) is
  begin
    if x > -ddtwopi and x < ddtwopi then
      q := Create(integer(0)); r := x;
    elsif x < 0.0 then
      PosDivModTwoPi(-x,q,r);
      Min(q); r := -r;
    else
      PosDivModTwoPi(x,q,r);
    end if;
  end DivModTwoPi;

  procedure DivModTwoPi ( x : in quad_double;
                          q : out Integer_Number; r : out quad_double ) is
  begin
    if x > -qdtwopi and x < qdtwopi then
      q := Create(integer(0)); r := x;
    elsif x < 0.0 then
      PosDivModTwoPi(-x,q,r);
      Min(q); r := -r;
    else
      PosDivModTwoPi(x,q,r);
    end if;
  end DivModTwoPi;

  procedure DivModTwoPi ( x : in Floating_Number; d : in natural32;
                          q : out Integer_Number; r : out Floating_Number ) is

    twopi : Floating_Number := Multprec_Floating_Constants.TwoPi(2*d);
    mintwopi : Floating_Number := -twopi;
    minx : Floating_Number;

  begin
    if x > mintwopi and x < twopi then
      q := Create(integer(0)); Copy(x,r);
    elsif x < 0.0 then
      minx := -x;
      PosDivModTwoPi(minx,twopi,d,q,r);
      Min(q); Min(r); Clear(minx);
    else
      PosDivModTwoPi(x,twopi,d,q,r);
    end if;
    Clear(twopi); Clear(mintwopi);
  end DivModTwoPi;

  function Polar_Exponentiation_ModTwoPi
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    r : constant double_float := Standard_Complex_Numbers_Polar.Radius(x);
    a : constant double_float := Standard_Complex_Numbers_Polar.Angle(x);
    s : constant double_float := r**integer(e);
    f : constant double_float := double_float(e)*a;
    q : integer32;
    b,re,im : double_float;
    res : Complex_Number;

  begin
    DivModTwoPi(f,q,b);
    re := s*COS(b);
    im := s*SIN(b);
    res := Create(re,im);
    return res;
  end Polar_Exponentiation_ModTwoPi;

  function Polar_Exponentiation_ModTwoPi_of_Unit
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    a : constant double_float := Standard_Complex_Numbers_Polar.Angle(x);
    f : constant double_float := double_float(e)*a;
    q : integer32;
    b,re,im : double_float;
    res : Complex_Number;

  begin
    DivModTwoPi(f,q,b);
    re := COS(b);
    im := SIN(b);
    res := Create(re,im);
    return res;
  end Polar_Exponentiation_ModTwoPi_of_Unit;

  function DoblDobl_Polar_Exponentiation_ModTwoPi_of_Unit
             ( x : Complex_Number; e : Integer_Number )
             return Complex_Number is

  -- DESCRIPTION :
  --   Auxiliary function to Polar_Exponentiation_ModTwoPi_of_Unit
  --   in case the number of decimal places of e is less than 32,
  --   using double double numbers in the mod 2*pi calculation.

    a : constant double_float := Standard_Complex_Numbers_Polar.Angle(x);
    dd_e : constant double_double := to_double_double(e);
    f : constant double_double := dd_e*a;
    q : Integer_Number;
    dd_b : double_double;
    b,re,im : double_float;
    res : Complex_Number;

  begin
    DivModTwoPi(f,q,dd_b);
    b := Double_Double_Numbers.hi_part(dd_b);
    re := COS(b);
    im := SIN(b);
    res := Create(re,im);
    Clear(q);
    return res;
  end DoblDobl_Polar_Exponentiation_ModTwoPi_of_Unit;

  function QuadDobl_Polar_Exponentiation_ModTwoPi_of_Unit
             ( x : Complex_Number; e : Integer_Number )
             return Complex_Number is

  -- DESCRIPTION :
  --   Auxiliary function to Polar_Exponentiation_ModTwoPi_of_Unit
  --   in case the number of decimal places of e is less than 64,
  --   using quad double numbers in the mod 2*pi calculation.
  --   Although the function also applies for larger values of e,
  --   the accuracy of the result is no longer certain.

    a : constant double_float := Standard_Complex_Numbers_Polar.Angle(x);
    qd_e : constant quad_double := to_quad_double(e);
    f : constant quad_double := qd_e*a;
    q : Integer_Number;
    qd_b : quad_double;
    b,re,im : double_float;
    res : Complex_Number;

  begin
    DivModTwoPi(f,q,qd_b);
    b := Quad_Double_Numbers.hihi_part(qd_b);
    re := COS(b);
    im := SIN(b);
    res := Create(re,im);
    Clear(q);
    return res;
  end QuadDobl_Polar_Exponentiation_ModTwoPi_of_Unit;

  function Multprec_Polar_Exponentiation_ModTwoPi_of_Unit
             ( x : Complex_Number; e : Integer_Number; d : natural32 )
             return Complex_Number is

  -- DESCRIPTION :
  --   Auxiliary function to Polar_Exponentiation_ModTwoPi_of_Unit
  --   in case the number of decimal places of e equals d,
  --   using multiprecision numbers in the mod 2*pi calculation.

    a : constant double_float := Standard_Complex_Numbers_Polar.Angle(x);
    mp_e : Floating_Number := Create(e);
    f : Floating_Number := mp_e*a;
    q : Integer_Number;
    r : Floating_Number;
    b,re,im : double_float;
    res : Complex_Number;

  begin
    DivModTwoPi(f,d,q,r);
    b := Round(r);
    re := COS(b);
    im := SIN(b);
    res := Create(re,im);
    Clear(mp_e); Clear(f); Clear(q);
    return res;
  end Multprec_Polar_Exponentiation_ModTwoPi_of_Unit;

  function Polar_Exponentiation_ModTwoPi_of_Unit
             ( x : Complex_Number; e : Integer_Number )
             return Complex_Number is

    dp : constant natural32 := Decimal_Places(e);
    dp2 : constant natural32 := 2*dp; -- very important for accuracy !

  begin
    if dp2 < 32 then
      return DoblDobl_Polar_Exponentiation_ModTwoPi_of_Unit(x,e);
    elsif dp2 < 64 then
      return QuadDobl_Polar_Exponentiation_ModTwoPi_of_Unit(x,e);
    else
      return Multprec_Polar_Exponentiation_ModTwoPi_of_Unit(x,e,dp);
    end if;
  end Polar_Exponentiation_ModTwoPi_of_Unit;

end Standard_Complex_Exponentiation;
