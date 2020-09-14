with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Deca_Double_Constants;              use Deca_Double_Constants;
with DecaDobl_Mathematical_Functions;    use DecaDobl_Mathematical_Functions;

package body DecaDobl_Complex_Numbers_Polar is

  function Radius ( c : Complex_Number ) return deca_double is
  begin
    return SQRT(sqr(REAL_PART(c)) + sqr(IMAG_PART(c)));
  end Radius;

  function Angle ( c : Complex_Number ) return deca_double is
  begin
    return DecaDobl_Mathematical_Functions.Angle(IMAG_PART(c),REAL_PART(c));
  end Angle;

  function Root ( c : Complex_Number; n,i : natural32 )
                return Complex_Number is

    arg,radius_c,angle_c : deca_double;
    tmp : Complex_Number;
    one : constant deca_double := create(1.0);
    df_i : constant double_float := double_float(i);
    df_n : constant double_float := double_float(n);
    td_i : constant deca_double := create(df_i);
    td_n : constant deca_double := create(df_n);

  begin
    arg := (2.0 * PI * td_i) / td_n;
    radius_c := RADIUS(c)**(one/td_n);
    angle_c := ANGLE(c)/td_n;
    tmp := Create(radius_c)*Create(COS(angle_c),SIN(angle_c));
    return Create(COS(arg),SIN(arg))*tmp;
  end Root;

  function Polar_Exponentiation
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    r : constant deca_double := Radius(x);
    a : constant deca_double := Angle(x);
    s : constant deca_double := r**integer(e);
    td_e : constant deca_double := create(e);
    b : constant deca_double := a*td_e;
    re : constant deca_double := s*COS(b);
    im : constant deca_double := s*SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation;

  function Polar_Exponentiation
             ( x : Complex_Number; e : deca_double )
             return Complex_Number is

    r : constant deca_double := Radius(x);
    a : constant deca_double := Angle(x);
    s : constant deca_double := r**e;
    b : constant deca_double := a*e;
    re : constant deca_double := s*COS(b);
    im : constant deca_double := s*SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation;

end DecaDobl_Complex_Numbers_Polar;
