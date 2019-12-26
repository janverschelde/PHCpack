with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Double_Double_Constants;
with DoblDobl_Mathematical_Functions;

package body DoblDobl_Random_Numbers is

  function Random return double_double is

    res : double_double;
    rf : constant double_float := Standard_Random_Numbers.Random;

  begin
    res := Create(rf);
    return res;
  end Random;

  procedure Random_Double_Double
              ( seed : in out integer32; f : out double_double ) is

    rdf : double_float;

  begin
    Standard_Random_Numbers.Random_Double_Float(seed,rdf);
    f := Create(rdf);
  end Random_Double_Double;

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

  procedure Random_Complex_Number
              ( seed : in out integer32; c : out Complex_Number ) is

    rpc,ipc : double_double;

  begin
    Random_Double_Double(seed,rpc);
    Random_Double_Double(seed,ipc);
    c := Create(rpc,ipc);
  end Random_Complex_Number;

  function Random1 return Complex_Number is

    res : Complex_Number;
    arg : double_double := DoblDobl_Random_Numbers.Random;
    cs,sn : double_double;

  begin
    arg := arg*Double_Double_Constants.pi;
    cs := DoblDobl_Mathematical_Functions.cos(arg);
    sn := DoblDobl_Mathematical_Functions.sin(arg);
    res := create(cs,sn);
    return res;
  end Random1;

  procedure Random1_Complex_Number
              ( seed : in out integer32; c : out Complex_Number ) is

    arg,cs,sn : double_double;

  begin
    Random_Double_Double(seed,arg);
    arg := arg*Double_Double_Constants.pi;
    cs := DoblDobl_Mathematical_Functions.cos(arg);
    sn := DoblDobl_Mathematical_Functions.sin(arg);
    c := create(cs,sn);
  end Random1_Complex_Number;

  function Random_Magnitude ( m : natural32 ) return Complex_Number is

    res : Complex_Number := Random1;
    r : constant double_float := Standard_Random_Numbers.Random_Magnitude(m);
    r_dd : constant double_double := create(r);

  begin
    res := r_dd*res;
    return res;
  end Random_Magnitude;

end DoblDobl_Random_Numbers;
