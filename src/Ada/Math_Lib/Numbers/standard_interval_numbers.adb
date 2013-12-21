package body Standard_Interval_Numbers is

-- AUXILIARIES :

  function min ( a,b : double_float ) return double_float is

  -- DESCRIPTION :
  --   Returns the minimum of a and b.

  begin
    if a < b
     then return a;
     else return b;
    end if;
  end min;

  function max ( a,b : double_float ) return double_float is

  -- DESCRIPTION :
  --   Returns the maximum of a and b.

  begin
    if a > b
     then return a;
     else return b;
    end if;
  end max;

  function min ( a,b,c,d : double_float ) return double_float is

  -- DESCRIPTION :
  --   Returns the minimum of a, b, c, and d.

  begin
    return min(min(min(a,b),c),d);
  end min;

  function max ( a,b,c,d : double_float ) return double_float is

  -- DESCRIPTION :
  --   Returns the maximum of a, b, c, and d.

  begin
    return max(max(max(a,b),c),d);
  end max;

-- CREATORS :

  function Create ( i : integer ) return Interval is
  begin
    return Create(double_float(i),double_float(i));
  end Create;

  function Create ( n : natural32 ) return Interval is
  begin
    return Create(double_float(n),double_float(n));
  end Create;

  function Create ( i : integer32 ) return Interval is
  begin
    return Create(double_float(i),double_float(i));
  end Create;

  function Create ( x : double_float ) return Interval is
  begin
    return Create(x,x);
  end Create;

  function Create ( a,b : double_float ) return Interval is

    res : Interval;
  
  begin
    res.a := a;
    res.b := b;
    return res;
  end Create;

-- COPY and EQUAL :

  function Equal ( x,y : Interval ) return boolean is
  begin
    if x.a /= y.a then
      return false;
    elsif x.b /= y.b then
      return false;
    else
      return true;
    end if;
  end Equal;

  procedure Copy ( x : in Interval; y : in out Interval ) is
  begin
    y.a := x.a;
    y.b := x.a;
  end Copy;

-- SELECTORS :

  function Left ( i : Interval ) return double_float is
  begin
    return i.a;
  end Left;

  function Right ( i : Interval ) return double_float is
  begin
    return i.b;
  end Right;

  function Width ( i : Interval ) return double_float is
  begin
    return i.b - i.a;
  end Width;

  function Width ( i : Interval ) return Interval is

    w : constant double_float := Width(i);

  begin
    return Create(w,w);
  end Width;

  function Middle ( i : Interval ) return double_float is
  begin
    return (i.a + i.b)/2.0;
  end Middle;

  function "<" ( x,y : Interval ) return boolean is
  begin
    return Middle(x) < Middle(y);
  end "<";

  function "<" ( x : Interval; y : double_float ) return boolean is
  begin
    return Middle(x) < y;
  end "<";

  function "<" ( x : double_float; y : Interval ) return boolean is
  begin
    return x < Middle(y);
  end "<";

  function ">" ( x,y : Interval ) return boolean is
  begin
    return Middle(x) > Middle(y);
  end ">";

  function ">" ( x : Interval; y : double_float ) return boolean is
  begin
    return Middle(x) > y;
  end ">";

  function ">" ( x : double_float; y : Interval ) return boolean is
  begin
    return x > Middle(y);
  end ">";

-- ARITHMETICAL OPERATIONS :

  function "+" ( x : Interval ) return Interval is
  begin
    return Create(x.a,x.b);
  end "+";

  function "+" ( x,y : Interval ) return Interval is
  begin
    return Create(x.a + y.a,x.b + y.b);
  end "+";

  function "+" ( x : Interval; y : double_float ) return Interval is
  begin
    return Create(x.a + y,x.b + y);
  end "+";

  function "+" ( x : double_float; y : Interval ) return Interval is
  begin
    return Create(x + y.a,x + y.b);
  end "+";

  procedure Add ( x : in out Interval; y : in Interval ) is
  begin
    x.a := x.a + y.a;
    x.b := x.b + y.b;
  end Add;

  procedure Add ( x : in out Interval; y : in double_float ) is
  begin
    x.a := x.a + y;
    x.b := x.b + y;
  end Add;

  function "-" ( x : Interval ) return Interval is
  begin
    return Create(-x.b,-x.a);
  end "-";

  procedure Min ( x : in out Interval ) is

    xa : constant double_float := x.a;

  begin
    x.a := -x.b; 
    x.b := -xa;
  end Min;

  function "-" ( x,y : Interval ) return Interval is
  begin
    return Create(x.a - y.b,x.b - y.a);
  end "-";

  function "-" ( x : Interval; y : double_float ) return Interval is
  begin
    return Create(x.a - y,x.b - y);
  end "-";

  function "-" ( x : double_float; y : Interval ) return Interval is
  begin
    return Create(x - y.b,x - y.a);
  end "-";

  procedure Sub ( x : in out Interval; y : in Interval ) is
  begin
    x.a := x.a - y.b;
    x.b := x.b - y.a;
  end Sub;

  procedure Sub ( x : in out Interval; y : in double_float ) is
  begin
    x.a := x.a - y;
    x.b := x.b - y;
  end Sub;

  function "*" ( x,y : Interval ) return Interval is

    ac : constant double_float := x.a*y.a;
    ad : constant double_float := x.a*y.b;
    bc : constant double_float := x.b*y.a;
    bd : constant double_float := x.b*y.b;

  begin
    return Create(min(ac,ad,bc,bd),max(ac,ad,bc,bd));
  end "*";

  function "*" ( x : Interval; y : double_float ) return Interval is
  begin
    if y >= 0.0 then
      return Create(y*x.a,y*x.b);
    else
      return Create(y*x.b,y*x.a);
    end if;
  end "*";

  function "*" ( x : double_float; y : Interval ) return Interval is
  begin
    if x >= 0.0 then
      return Create(x*y.a,x*y.b);
    else
      return Create(x*y.b,x*y.a);
    end if;
  end "*";

  procedure Mul ( x : in out Interval; y : in Interval ) is

    ac : constant double_float := x.a*y.a;
    ad : constant double_float := x.a*y.b;
    bc : constant double_float := x.b*y.a;
    bd : constant double_float := x.b*y.b;

  begin
    x.a := min(ac,ad,bc,bd);
    x.b := max(ac,ad,bc,bd);
  end Mul;

  procedure Mul ( x : in out Interval; y : in double_float ) is

    xa : double_float;
 
  begin
    if y >= 0.0 then
      x.a := y*x.a;
      x.b := y*x.b;
    else
      xa := x.a;
      x.a := y*x.b;
      x.b := y*xa;
    end if;
  end Mul;

  function "/" ( x,y : Interval ) return Interval is

    ac : constant double_float := x.a/y.a;
    ad : constant double_float := x.a/y.b;
    bc : constant double_float := x.b/y.a;
    bd : constant double_float := x.b/y.b;

  begin
    return Create(min(ac,ad,bc,bd),max(ac,ad,bc,bd));
  end "/";

  function "/" ( x : Interval; y : double_float ) return Interval is
  begin
    if y >= 0.0 then
      return Create(x.a/y,x.b/y);
    else
      return Create(x.b/y,x.a/y);
    end if;
  end "/";

  function "/" ( x : double_float; y : Interval ) return Interval is

    c : constant double_float := x/y.a;
    d : constant double_float := x/y.b;

  begin
    return Create(min(c,d),max(c,d));
  end "/";

  procedure Div ( x : in out Interval; y : in Interval ) is

    ac : constant double_float := x.a/y.a;
    ad : constant double_float := x.a/y.b;
    bc : constant double_float := x.b/y.a;
    bd : constant double_float := x.b/y.b;

  begin
    x.a := min(ac,ad,bc,bd);
    x.b := max(ac,ad,bc,bd);
  end Div;

  procedure Div ( x : in out Interval; y : in double_float ) is

    xa : double_float;

  begin
    if y >= 0.0 then
      x.a := x.a/y;
      x.b := x.b/y;
    else
      xa := x.a;
      x.a := x.b/y;
      x.b := xa/y;
    end if;
  end Div;

-- DESTRUCTOR :

  procedure Clear ( x : in out Interval ) is
  begin
    x.a := 0.0;
    x.b := 0.0;
  end Clear;

end Standard_Interval_Numbers;
