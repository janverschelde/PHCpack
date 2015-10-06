with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Constants;              use Quad_Double_Constants;
with QuadDobl_Mathematical_Functions;    use QuadDobl_Mathematical_Functions;

package body QuadDobl_Complex_Numbers_Polar is

  function Radius ( c : Complex_Number ) return quad_double is
  begin
    return SQRT(sqr(REAL_PART(c)) + sqr(IMAG_PART(c)));
  end Radius;

  function Angle ( c : Complex_Number ) return quad_double is
  begin
    return QuadDobl_Mathematical_Functions.Angle(IMAG_PART(c),REAL_PART(c));
  end Angle;

  function Root ( c : Complex_Number; n,i : natural32 )
                return Complex_Number is

    arg,radius_c,angle_c : quad_double;
    tmp : Complex_Number;
    one : constant quad_double := create(1.0);
    df_i : constant double_float := double_float(i);
    df_n : constant double_float := double_float(n);
    dd_i : constant quad_double := create(df_i);
    dd_n : constant quad_double := create(df_n);

  begin
    arg := (2.0 * PI * dd_i) / dd_n;
    radius_c := RADIUS(c)**(one/dd_n);
    angle_c := ANGLE(c)/dd_n;
    tmp := Create(radius_c)*Create(COS(angle_c),SIN(angle_c));
    return Create(COS(arg),SIN(arg))*tmp;
  end Root;

end QuadDobl_Complex_Numbers_Polar;
