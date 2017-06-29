package body standard_floating_numbers is

  function Create ( i : integer ) return single_float is
  begin
    return single_float(i);
  end Create;

  function Create ( i : integer ) return double_float is
  begin
    return double_float(i);
  end Create;

  function Create ( i : integer32 ) return single_float is
  begin
    return single_float(i);
  end Create;

  function Create ( i : integer32 ) return double_float is
  begin
    return double_float(i);
  end Create;

  function Create ( n : natural32 ) return single_float is
  begin
    return single_float(n);
  end Create;

  function Create ( n : natural32 ) return double_float is
  begin
    return double_float(n);
  end Create;

  function Equal ( a,b : single_float ) return boolean is
  begin
    return (a = b);
  end equal;

  function Equal ( a,b : double_float ) return boolean is
  begin
    return (a = b);
  end Equal;

  function AbsVal ( a : single_float ) return single_float is
  begin
    return ABS(a);
  end AbsVal;

  function AbsVal ( a : double_float ) return double_float is
  begin
    return ABS(a);
  end AbsVal;

  procedure Copy ( a : in single_float; b : in out single_float ) is
  begin
    b := a;
  end Copy;

  procedure Copy ( a : in double_float; b : in out double_float ) is
  begin
    b := a;
  end Copy;

  procedure Add ( a : in out single_float; b : in single_float ) is
  begin
    a := a + b;
  end Add;

  procedure Add ( a : in out double_float; b : in double_float ) is
  begin
    a := a + b;
  end Add;

  procedure Sub ( a : in out single_float; b : in single_float ) is
  begin
    a := a - b;
  end Sub;

  procedure Sub ( a : in out double_float; b : in double_float ) is
  begin
    a := a - b;
  end Sub;

  procedure Min  ( a : in out single_float ) is
  begin
    a := -a;
  end Min;

  procedure Min ( a : in out double_float ) is
  begin
    a := -a;
  end Min;

  procedure Mul ( a : in out single_float; b : in single_float ) is
  begin
    a := a*b;
  end Mul;

  procedure Mul ( a : in out double_float; b : in double_float ) is
  begin
    a := a*b;
  end Mul;

  procedure Div ( a : in out single_float; b : in single_float ) is
  begin
    a := a/b;
  end Div;

  procedure Div ( a : in out double_float; b : in double_float ) is
  begin
    a := a/b;
  end Div;

  function Is_Valid ( x : single_float ) return boolean is
  begin
    if x > 0.0 then
      return true;
    elsif x < 0.0 then
      return true;
    elsif x = 0.0 then
      return true;
    else
      return false;
    end if;
  end Is_Valid;

  function Is_Valid ( x : double_float ) return boolean is
  begin
    if x > 0.0 then
      return true;
    elsif x < 0.0 then
      return true;
    elsif x = 0.0 then
      return true;
    else
      return false;
    end if;
  end Is_Valid;

  procedure Clear ( a : in out single_float ) is
  begin
    a := 0.0;
  end Clear;

  procedure Clear ( a : in out double_float ) is
  begin
    a := 0.0;
  end Clear;

end Standard_Floating_Numbers;
