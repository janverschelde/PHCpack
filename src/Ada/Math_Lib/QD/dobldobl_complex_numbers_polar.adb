with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Constants;            use Double_Double_Constants;
with DoblDobl_Mathematical_Functions;    use DoblDobl_Mathematical_Functions;

package body DoblDobl_Complex_Numbers_Polar is

  function Radius ( c : Complex_Number ) return double_double is
  begin
    return SQRT(sqr(REAL_PART(c)) + sqr(IMAG_PART(c)));
  end Radius;

  function Angle ( c : Complex_Number ) return double_double is
  begin
    return DoblDobl_Mathematical_Functions.Angle(IMAG_PART(c),REAL_PART(c));
  end Angle;

  function Root ( c : Complex_Number; n,i : natural32 )
                return Complex_Number is

    arg,radius_c,angle_c : double_double;
    tmp : Complex_Number;
    one : constant double_double := create(1.0);
    df_i : constant double_float := double_float(i);
    df_n : constant double_float := double_float(n);
    dd_i : constant double_double := create(df_i);
    dd_n : constant double_double := create(df_n);

  begin
    arg := (2.0 * PI * dd_i) / dd_n;
    radius_c := RADIUS(c)**(one/dd_n);
    angle_c := ANGLE(c)/dd_n;
    tmp := Create(radius_c)*Create(COS(angle_c),SIN(angle_c));
    return Create(COS(arg),SIN(arg))*tmp;
  end Root;

  function Polar_Exponentiation
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    r : constant double_double := Radius(x);
    a : constant double_double := Angle(x);
    s : constant double_double := r**integer(e);
    dd_e : constant double_double := create(e);
    b : constant double_double := a*dd_e;
    re : constant double_double := s*COS(b);
    im : constant double_double := s*SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation;

  function Polar_Exponentiation
             ( x : Complex_Number; e : double_double )
             return Complex_Number is

    r : constant double_double := Radius(x);
    a : constant double_double := Angle(x);
    s : constant double_double := r**e;
    b : constant double_double := a*e;
    re : constant double_double := s*COS(b);
    im : constant double_double := s*SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation;

end DoblDobl_Complex_Numbers_Polar;
