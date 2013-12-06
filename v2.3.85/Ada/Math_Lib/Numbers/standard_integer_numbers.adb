package body Standard_Integer_Numbers is

  function Create ( i : integer ) return integer32 is
  begin
    return integer32(i);
  end Create;

  function Create ( i : integer ) return integer64 is
  begin
    return integer64(i);
  end Create;

  function Equal ( a,b : integer32 ) return boolean is
  begin
    return (a = b);
  end Equal;

  function Equal ( a,b : integer64 ) return boolean is
  begin
    return (a = b);
  end Equal;

  procedure Copy ( a : in integer32; b : in out integer32 ) is
  begin
    b := a;
  end Copy;

  procedure Copy ( a : in integer64; b : in out integer64 ) is
  begin
    b := a;
  end Copy;

  procedure Add ( a : in out integer32; b : in integer32 ) is
  begin
    a := a+b;
  end Add;

  procedure Add ( a : in out integer64; b : in integer64 ) is
  begin
    a := a+b;
  end Add;

  procedure Sub ( a : in out integer32; b : in integer32 ) is
  begin
    a := a-b;
  end Sub;

  procedure Sub ( a : in out integer64; b : in integer64 ) is
  begin
    a := a-b;
  end Sub;

  procedure Min ( a : in out integer32 ) is
  begin
    a := -a;
  end Min;

  procedure Min ( a : in out integer64 ) is
  begin
    a := -a;
  end Min;

  procedure Mul ( a : in out integer32; b : in integer32 ) is
  begin
    a := a*b;
  end Mul;

  procedure Mul ( a : in out integer64; b : in integer64 ) is
  begin
    a := a*b;
  end Mul;

  function Rmd ( a,b : integer32 ) return integer32 is
  begin
    return a mod b;
  end Rmd;

  function Rmd ( a,b : integer64 ) return integer64 is
  begin
    return a mod b;
  end Rmd;

  procedure Rmd ( a : in out integer32; b : in integer32 ) is
  begin
    a := a mod b;
  end Rmd;

  procedure Rmd ( a : in out integer64; b : in integer64 ) is
  begin
    a := a mod b;
  end Rmd;

  procedure Div ( a : in out integer32; b : in integer32 ) is
  begin
    a := a/b;
  end Div;

  procedure Div ( a : in out integer64; b : in integer64 ) is
  begin
    a := a/b;
  end Div;

  procedure Div ( a,b : in integer32; q : out integer32; r : out integer32 ) is
  begin
    q := a/b;
    r := a mod b;
  end Div;

  procedure Div ( a,b : in integer64; q : out integer64; r : out integer64 ) is
  begin
    q := a/b;
    r := a mod b;
  end Div;

  procedure Div ( a : in out integer32;
                  b : in integer32; r : out integer32 ) is
  begin
    r := a mod b;
    a := a/b;
  end Div;

  procedure Div ( a : in out integer64;
                  b : in integer64; r : out integer64 ) is
  begin
    r := a mod b;
    a := a/b;
  end Div;

  procedure Clear ( a : in out integer32 ) is
  begin
    a := 0;
  end Clear;

  procedure Clear ( a : in out integer64 ) is
  begin
    a := 0;
  end Clear;

end Standard_Integer_Numbers;
