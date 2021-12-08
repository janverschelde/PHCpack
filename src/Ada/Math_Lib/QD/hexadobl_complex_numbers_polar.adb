with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Hexa_Double_Constants;              use Hexa_Double_Constants;
with HexaDobl_Mathematical_Functions;    use HexaDobl_Mathematical_Functions;

package body HexaDobl_Complex_Numbers_Polar is

  function Radius ( c : Complex_Number ) return hexa_double is
  begin
    return SQRT(sqr(REAL_PART(c)) + sqr(IMAG_PART(c)));
  end Radius;

  function Angle ( c : Complex_Number ) return hexa_double is
  begin
    return HexaDobl_Mathematical_Functions.Angle(IMAG_PART(c),REAL_PART(c));
  end Angle;

  function Root ( c : Complex_Number; n,i : natural32 )
                return Complex_Number is

    arg,radius_c,angle_c : hexa_double;
    tmp : Complex_Number;
    one : constant hexa_double := create(1.0);
    df_i : constant double_float := double_float(i);
    df_n : constant double_float := double_float(n);
    td_i : constant hexa_double := create(df_i);
    td_n : constant hexa_double := create(df_n);

  begin
    arg := (2.0 * PI * td_i) / td_n;
    radius_c := RADIUS(c)**(one/td_n);
    angle_c := ANGLE(c)/td_n;
    tmp := Create(radius_c)*Create(COS(angle_c),SIN(angle_c));
    return Create(COS(arg),SIN(arg))*tmp;
  end Root;

  function Polar_Exponentiation
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    r : constant hexa_double := Radius(x);
    a : constant hexa_double := Angle(x);
    s : constant hexa_double := r**integer(e);
    td_e : constant hexa_double := create(e);
    b : constant hexa_double := a*td_e;
    re : constant hexa_double := s*COS(b);
    im : constant hexa_double := s*SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation;

  function Polar_Exponentiation
             ( x : Complex_Number; e : hexa_double )
             return Complex_Number is

    r : constant hexa_double := Radius(x);
    a : constant hexa_double := Angle(x);
    s : constant hexa_double := r**e;
    b : constant hexa_double := a*e;
    re : constant hexa_double := s*COS(b);
    im : constant hexa_double := s*SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation;

end HexaDobl_Complex_Numbers_Polar;
