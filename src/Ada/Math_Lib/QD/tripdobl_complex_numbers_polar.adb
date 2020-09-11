with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Triple_Double_Constants;            use Triple_Double_Constants;
with TripDobl_Mathematical_Functions;    use TripDobl_Mathematical_Functions;

package body TripDobl_Complex_Numbers_Polar is

  function Radius ( c : Complex_Number ) return triple_double is
  begin
    return SQRT(sqr(REAL_PART(c)) + sqr(IMAG_PART(c)));
  end Radius;

  function Angle ( c : Complex_Number ) return triple_double is
  begin
    return TripDobl_Mathematical_Functions.Angle(IMAG_PART(c),REAL_PART(c));
  end Angle;

  function Root ( c : Complex_Number; n,i : natural32 )
                return Complex_Number is

    arg,radius_c,angle_c : triple_double;
    tmp : Complex_Number;
    one : constant triple_double := create(1.0);
    df_i : constant double_float := double_float(i);
    df_n : constant double_float := double_float(n);
    td_i : constant triple_double := create(df_i);
    td_n : constant triple_double := create(df_n);

  begin
    arg := (2.0 * PI * td_i) / td_n;
    radius_c := RADIUS(c)**(one/td_n);
    angle_c := ANGLE(c)/td_n;
    tmp := Create(radius_c)*Create(COS(angle_c),SIN(angle_c));
    return Create(COS(arg),SIN(arg))*tmp;
  end Root;

  function Polar_Exponentiation
             ( x : Complex_Number; e : integer32 ) return Complex_Number is

    r : constant triple_double := Radius(x);
    a : constant triple_double := Angle(x);
    s : constant triple_double := r**integer(e);
    td_e : constant triple_double := create(e);
    b : constant triple_double := a*td_e;
    re : constant triple_double := s*COS(b);
    im : constant triple_double := s*SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation;

  function Polar_Exponentiation
             ( x : Complex_Number; e : triple_double )
             return Complex_Number is

    r : constant triple_double := Radius(x);
    a : constant triple_double := Angle(x);
    s : constant triple_double := r**e;
    b : constant triple_double := a*e;
    re : constant triple_double := s*COS(b);
    im : constant triple_double := s*SIN(b);
    res : constant Complex_Number := Create(re,im);

  begin
    return res;
  end Polar_Exponentiation;

end TripDobl_Complex_Numbers_Polar;
