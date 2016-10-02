with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Standard_Complex_Numbers_Polar is

  function Radius ( c : Complex_Number ) return double_float is
  begin
    return SQRT(REAL_PART(c)**2 + IMAG_PART(c)**2);
  end Radius;

  function Angle ( c : Complex_Number ) return double_float is
  begin
    return Standard_Mathematical_Functions.Angle(IMAG_PART(c),REAL_PART(c));
  end Angle;

  function Root ( c : Complex_Number; n,i : natural32 ) return Complex_Number is

    arg,radius_c,angle_c : double_float;
    tmp : Complex_Number;

  begin
    arg := (2.0 * PI * double_float(i)) / double_float(n);
    radius_c := RADIUS(c)**(1.0/double_float(n));
    angle_c := ANGLE(c)/double_float(n);
    tmp := Create(radius_c)*Create(COS(angle_c),SIN(angle_c));
    return Create(COS(arg),SIN(arg))*tmp;
  end Root;

  function Polar_Exponentiation
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    r : constant double_float := Radius(x);
    a : constant double_float := Angle(x);
    s : constant double_float := r**integer(e);
    b : constant double_float := a*double_float(e);
    re : constant double_float := s*COS(b);
    im : constant double_float := s*SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation;

  function Polar_Exponentiation
             ( x : Complex_Number; e : double_float )
             return Complex_Number is

    r : constant double_float := Radius(x);
    a : constant double_float := Angle(x);
    s : constant double_float := r**e;
    b : constant double_float := a*e;
    re : constant double_float := s*COS(b);
    im : constant double_float := s*SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation;

  function Polar_Exponentiation_of_Unit
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    a : constant double_float := Angle(x);
    b : constant double_float := a*double_float(e);
    re : constant double_float := COS(b);
    im : constant double_float := SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation_of_Unit;

end Standard_Complex_Numbers_Polar;
