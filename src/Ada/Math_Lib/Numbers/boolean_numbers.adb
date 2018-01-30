package body Boolean_Numbers is

  function Create ( i : integer ) return boolean is
  begin
    if i = 0
     then return false;
     else return true;
    end if;
  end Create;

  function Equal ( a,b : boolean ) return boolean is
  begin
    return (a = b);
  end Equal;

  procedure Copy ( a : in boolean; b : in out boolean ) is
  begin
    b := a;
  end Copy;

  function "+" ( a,b : boolean ) return boolean is
  begin
    if not b then  -- b is zero
      return a;
    else           -- b is one
      if a
       then return false; -- 1 + 1 mod 2 is 0
       else return true;  -- 0 + 1 is 1
      end if;
    end if;
  end "+";

  function "+" ( a : boolean )   return boolean is
  begin
    return a;
  end "+";

  function "-" ( a,b : boolean ) return boolean is
  begin
    if not b then
      return a;
    else
      if a
       then return false;
       else return true;
      end if;
    end if;
  end "-";

  function "-" ( a : boolean )   return boolean is
  begin
    return a;
  end "-";

  function "*" ( a,b : boolean ) return boolean is
  begin
    if b
     then return a;
     else return false;
    end if;
  end "*";

  procedure Add ( a : in out boolean; b : in boolean ) is
  begin
    if b then
      if a
       then a := false; -- 1+1 mod 2 is 0
       else a := true;  -- 1+0 mod 2 is 1
      end if;
    end if;
  end Add;

  procedure Sub ( a : in out boolean; b : in boolean ) is
  begin
    if b then
      if a
       then a := false; -- 1-1 = 0
       else a := true;  -- 0-1 mod 2 = 1
      end if;
    end if;
  end Sub;

  procedure Min ( a : in out boolean ) is
  begin
    null; -- as -1 modulo 2 equals 1 and 0 modulo 2 is 0
  end Min;

  procedure Mul ( a : in out boolean; b : in boolean ) is
  begin
    if not b
     then a := false;
    end if;
  end Mul;

  function Rmd ( a,b : boolean ) return boolean is
  begin
    return a; -- undefined if b is false
  end Rmd;

  procedure Rmd ( a : in out boolean; b : in boolean ) is
  begin
    null;
  end Rmd;

  procedure Div ( a : in out boolean; b : in boolean ) is
  begin
    null;
  end Div;

  procedure Div ( a,b : in boolean; q : out boolean; r : out boolean ) is
  begin
    q := a;
    r := a;
  end Div;

  procedure Div ( a : in out boolean;
                  b : in boolean; r : out boolean ) is
  begin
    r := a;
  end Div;

  procedure Clear ( a : in out boolean ) is
  begin
    a := false;
  end Clear;

end Boolean_Numbers;
