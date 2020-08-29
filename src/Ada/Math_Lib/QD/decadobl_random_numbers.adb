with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;

package body DecaDobl_Random_Numbers is

  function Random return deca_double is

    res : deca_double;
    first : constant double_float := Standard_Random_Numbers.Random; 
    second : double_float := Standard_Random_Numbers.Random; 
    eps : constant double_float := 2.0**(-52);
    multiplier : double_float := eps;

  begin
    res := create(first);
    res := res + eps*second;
    for k in 3..10 loop
      multiplier := eps*multiplier;
      second := Standard_Random_Numbers.Random;
      res := res + multiplier*second;
    end loop;
    return res;
  end Random;

  procedure Random_Deca_Double
              ( seed : in out integer32; f : out Deca_double ) is

    res : deca_double;
    first,second : double_float;
    eps : constant double_float := 2.0**(-52);
    multiplier : double_float := eps;

  begin
    Standard_Random_Numbers.Random_Double_Float(seed,first);
    Standard_Random_Numbers.Random_Double_Float(seed,second);
    res := create(first);
    res := res + eps*second;
    for k in 3..10 loop
      multiplier := eps*multiplier;
      Standard_Random_Numbers.Random_Double_Float(seed,second);
      res := res + multiplier*second;
    end loop;
    f := res;
  end Random_Deca_Double;

  function Random return Complex_Number is

    res : Complex_Number;
    realpart : constant deca_double := Random; 
    imagpart : constant deca_double := Random; 

  begin
    res := Create(realpart,imagpart);
    return res;
  end Random;

  procedure Random_Complex_Number
              ( seed : in out integer32; c : out Complex_Number ) is

    realpart,imagpart : deca_double;

  begin
    Random_Deca_Double(seed,realpart);
    Random_Deca_Double(seed,imagpart);
    c := Create(realpart,imagpart);
  end Random_Complex_Number;

end DecaDobl_Random_Numbers; 
