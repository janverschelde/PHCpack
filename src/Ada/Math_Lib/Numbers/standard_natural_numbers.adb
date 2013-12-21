package body Standard_Natural_Numbers is

  function Create ( i : integer ) return natural32 is
  begin
    return natural32(i);
  end Create;

  function Create ( i : integer ) return natural64 is
  begin
    return natural64(i);
  end Create;

  function Equal ( a,b : natural32 ) return boolean is
  begin
    return (a = b);
  end Equal;

  function Equal ( a,b : natural64 ) return boolean is
  begin
    return (a = b);
  end Equal;

  procedure Copy ( a : in natural32; b : in out natural32 ) is
  begin
    b := a;
  end Copy;

  procedure Copy ( a : in natural64; b : in out natural64 ) is
  begin
    b := a;
  end Copy;

  procedure Add ( a : in out natural32; b : in natural32 ) is
  begin
    a := a+b;
  end Add;

  procedure Add ( a : in out natural64; b : in natural64 ) is
  begin
    a := a+b;
  end Add;

  procedure Sub ( a : in out natural32; b : in natural32 ) is
  begin
    a := a-b;
  end Sub;

  procedure Sub ( a : in out natural64; b : in natural64 ) is
  begin
    a := a-b;
  end Sub;

  procedure Min ( a : in out natural32 ) is
  begin
    a := -a;
  end Min;

  procedure Min ( a : in out natural64 ) is
  begin
    a := -a;
  end Min;

  procedure Mul ( a : in out natural32; b : in natural32 ) is
  begin
    a := a*b;
  end Mul;

  procedure Mul ( a : in out natural64; b : in natural64 ) is
  begin
    a := a*b;
  end Mul;

  procedure Clear ( a : in out natural32 ) is
  begin
    a := 0;
  end Clear;

  procedure Clear ( a : in out natural64 ) is
  begin
    a := 0;
  end Clear;

end Standard_Natural_Numbers;
