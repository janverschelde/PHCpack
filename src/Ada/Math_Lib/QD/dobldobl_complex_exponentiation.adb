with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Double_Double_Constants;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Mathematical_Functions;    use DoblDobl_Mathematical_Functions;
with DoblDobl_Complex_Numbers_Polar;
with Quad_Double_Constants;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Constants;
with Multprec_DoblDobl_Convertors;       use Multprec_DoblDobl_Convertors;
with Multprec_QuadDobl_Convertors;       use Multprec_QuadDobl_Convertors;
with Standard_Complex_Exponentiation;

-- for testing purposes:
--with text_io,integer_io;           use text_io,integer_io;
--with Standard_Floating_Numbers_io; use Standard_Floating_Numbers_io;
--with Multprec_Integer_Numbers_io;  use Multprec_Integer_Numbers_io;
--with Multprec_Floating_Numbers_io; use Multprec_Floating_Numbers_io;
--with Double_Double_Numbers_io;     use Double_Double_Numbers_io;
--with Quad_Double_Numbers_io;       use Quad_Double_Numbers_io;

package body DoblDobl_Complex_Exponentiation is

  ddtwopi : constant double_double := Double_Double_Constants.twopi;
  qdtwopi : constant quad_double := Quad_Double_Constants.twopi;

  function Polar_Exponentiation_ModTwoPi
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    r : constant double_double := DoblDobl_Complex_Numbers_Polar.Radius(x);
    a : constant double_double := DoblDobl_Complex_Numbers_Polar.Angle(x);
    s : constant double_double := r**integer(e);
    f : constant double_double := Double_Double_Numbers.create(e)*a;
    q : Integer_Number;
    b,re,im : double_double;
    res : Complex_Number;

  begin
    Standard_Complex_Exponentiation.DivModTwoPi(f,q,b);
    re := s*COS(b);
    im := s*SIN(b);
    res := Create(re,im);
    Clear(q);
    return res;
  end Polar_Exponentiation_ModTwoPi;

  function Polar_Exponentiation_ModTwoPi_of_Unit
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    a : constant double_double := DoblDobl_Complex_Numbers_Polar.Angle(x);
    f : constant double_double := Double_Double_Numbers.create(e)*a;
    q : Integer_Number;
    b,re,im : double_double;
    res : Complex_Number;

  begin
    Standard_Complex_Exponentiation.DivModTwoPi(f,q,b);
    re := COS(b);
    im := SIN(b);
    res := Create(re,im);
    Clear(q);
    return res;
  end Polar_Exponentiation_ModTwoPi_of_Unit;

  function DoblDobl_Polar_Exponentiation_ModTwoPi_of_Unit
             ( x : Complex_Number; e : Integer_Number )
             return Complex_Number is

  -- DESCRIPTION :
  --   Auxiliary function to Polar_Exponentiation_ModTwoPi_of_Unit
  --   in case the number of decimal places of e is less than 32,
  --   using double double numbers in the mod 2*pi calculation.

    a : constant double_double := DoblDobl_Complex_Numbers_Polar.Angle(x);
    dd_e : constant double_double := to_double_double(e);
    f : constant double_double := dd_e*a;
    q : Integer_Number;
    b,re,im : double_double;
    res : Complex_Number;

  begin
    Standard_Complex_Exponentiation.DivModTwoPi(f,q,b);
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

    a : constant double_double := DoblDobl_Complex_Numbers_Polar.Angle(x);
    qd_e : constant quad_double := to_quad_double(e);
    f : constant quad_double := qd_e*a;
    q : Integer_Number;
    qd_b : quad_double;
    b,re,im : double_double;
    res : Complex_Number;

  begin
    Standard_Complex_Exponentiation.DivModTwoPi(f,q,qd_b);
    b := Quad_Double_Numbers.to_double_double(qd_b);
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

    a : constant double_double := DoblDobl_Complex_Numbers_Polar.Angle(x);
    mp_e : Floating_Number := Create(e);
    mp_a : Floating_Number := to_floating_number(a);
    f : Floating_Number := mp_e*mp_a;
    q : Integer_Number;
    r : Floating_Number;
    b,re,im : double_double;
    res : Complex_Number;

  begin
    Standard_Complex_Exponentiation.DivModTwoPi(f,d,q,r);
    b := to_double_double(r);
    re := COS(b);
    im := SIN(b);
    res := Create(re,im);
    Clear(mp_e); Clear(mp_a); Clear(f); Clear(q);
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

end DoblDobl_Complex_Exponentiation;
