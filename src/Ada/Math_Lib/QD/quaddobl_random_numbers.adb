with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Standard_Mathematical_Functions;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Constants;
with QuadDobl_Mathematical_Functions;

package body QuadDobl_Random_Numbers is

  function Random return quad_double is

    res : quad_double;
    rf : constant double_float := Standard_Random_Numbers.Random;

  begin
    res := Create(rf);
    return res;
  end Random;

  function Random_Magnitude ( m : natural32 ) return quad_double is

    res : quad_double := Random;
    r : constant double_float := Standard_Random_Numbers.Random_Magnitude(m);

  begin
    res := r*res;
    return res;
  end Random_Magnitude;

  function Random return Complex_Number is

    res : Complex_Number;
    rf1 : constant double_float := Standard_Random_Numbers.Random;
    rf2 : constant double_float := Standard_Random_Numbers.Random;
    rlp : constant quad_double := create(rf1);
    imp : constant quad_double := create(rf2);

  begin
    res := create(rlp,imp);
    return res;
  end Random;

  function Random1 return Complex_Number is

    res : Complex_Number;
    arg : quad_double := QuadDobl_Random_Numbers.Random;
    cs,sn : quad_double;

  begin
    arg := arg*Quad_Double_Constants.pi;
    cs := QuadDobl_Mathematical_Functions.cos(arg);
    sn := QuadDobl_Mathematical_Functions.sin(arg);
    res := create(cs,sn);
    return res;
  end Random1;

  function Random_Magnitude ( m : natural32 ) return Complex_Number is

    res : Complex_Number := Random1;
    r : constant double_float := Standard_Random_Numbers.Random_Magnitude(m);
    r_dd : constant double_double := create(r);
    r_qd : constant quad_double := create(r_dd);

  begin
    res := r_qd*res;
    return res;
  end Random_Magnitude;

end QuadDobl_Random_Numbers;
