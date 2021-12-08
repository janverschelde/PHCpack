with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Hexa_Double_Constants;
with HexaDobl_Mathematical_Functions;

package body HexaDobl_Random_Numbers is

  function Random return hexa_double is

    res : hexa_double;
    first : constant double_float := Standard_Random_Numbers.Random; 
    second : double_float := Standard_Random_Numbers.Random; 
    eps : constant double_float := 2.0**(-52);
    multiplier : double_float := eps;

  begin
    res := create(first);
    res := res + eps*second;
    for k in 3..16 loop
      multiplier := eps*multiplier;
      second := Standard_Random_Numbers.Random;
      res := res + multiplier*second;
    end loop;
    return res;
  end Random;

  procedure Random_Hexa_Double
              ( seed : in out integer32; f : out hexa_double ) is

    res : hexa_double;
    first,second : double_float;
    eps : constant double_float := 2.0**(-52);
    multiplier : double_float := eps;

  begin
    Standard_Random_Numbers.Random_Double_Float(seed,first);
    Standard_Random_Numbers.Random_Double_Float(seed,second);
    res := create(first);
    res := res + eps*second;
    for k in 3..16 loop
      multiplier := eps*multiplier;
      Standard_Random_Numbers.Random_Double_Float(seed,second);
      res := res + multiplier*second;
    end loop;
    f := res;
  end Random_Hexa_Double;

  function Random return Complex_Number is

    res : Complex_Number;
    realpart : constant hexa_double := Random; 
    imagpart : constant hexa_double := Random; 

  begin
    res := Create(realpart,imagpart);
    return res;
  end Random;

  procedure Random_Complex_Number
              ( seed : in out integer32; c : out Complex_Number ) is

    realpart,imagpart : hexa_double;

  begin
    Random_Hexa_Double(seed,realpart);
    Random_Hexa_Double(seed,imagpart);
    c := Create(realpart,imagpart);
  end Random_Complex_Number;

  function Random1 return Complex_Number is

    res : Complex_Number;
    arg : Hexa_double := HexaDobl_Random_Numbers.Random;
    cs,sn : hexa_double;

  begin
    arg := arg*Hexa_Double_Constants.pi;
    cs := HexaDobl_Mathematical_Functions.cos(arg);
    sn := HexaDobl_Mathematical_Functions.sin(arg);
    res := create(cs,sn);
    return res;
  end Random1;

  procedure Random1_Complex_Number
              ( seed : in out integer32; c : out Complex_Number ) is


    arg,cs,sn : hexa_double;

  begin
    Random_Hexa_Double(seed,arg);
    arg := arg*Hexa_Double_Constants.pi;
    cs := HexaDobl_Mathematical_Functions.cos(arg);
    sn := HexaDobl_Mathematical_Functions.sin(arg);
    c := create(cs,sn);
  end Random1_Complex_Number;

end HexaDobl_Random_Numbers; 
