with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Standard_Mathematical_Functions;

package body DoblDobl_Random_Numbers is

  function Random return double_double is

    res : double_double;
    rf : constant double_float := Standard_Random_Numbers.Random;

  begin
    res := Create(rf);
    return res;
  end Random;

  function Random_Magnitude ( m : natural32 ) return double_double is

    res : double_double := Random;
    r : constant double_float := Standard_Random_Numbers.Random_Magnitude(m);

  begin
    res := r*res;
    return res;
  end Random_Magnitude;

  function Random return Complex_Number is

    res : Complex_Number;
    rf1 : constant double_float := Standard_Random_Numbers.Random;
    rf2 : constant double_float := Standard_Random_Numbers.Random;
    rlp : constant double_double := create(rf1);
    imp : constant double_double := create(rf2);

  begin
    res := create(rlp,imp);
    return res;
  end Random;

  function Random1 return Complex_Number is

    res : Complex_Number;
    arg : double_float := Standard_Random_Numbers.Random;
    cs,sn : double_float;
    rlp,imp : double_double;

  begin
    arg := arg*Standard_Mathematical_Functions.PI;
    cs := Standard_Mathematical_Functions.cos(arg);
    sn := Standard_Mathematical_Functions.sin(arg);
    rlp := create(cs);
    imp := create(sn);
    res := create(rlp,imp);
    return res;
  end Random1;

  function Random_Magnitude ( m : natural32 ) return Complex_Number is

    res : Complex_Number := Random1;
    r : constant double_float := Standard_Random_Numbers.Random_Magnitude(m);
    r_dd : constant double_double := create(r);

  begin
    res := r_dd*res;
    return res;
  end Random_Magnitude;

end DoblDobl_Random_Numbers;
