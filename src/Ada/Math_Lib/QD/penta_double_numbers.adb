with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with Double_Double_Basics;
with Double_Double_Numbers;
with Fast_Double_Renormalizations;       use Fast_Double_Renormalizations;

package body Penta_Double_Numbers is

-- CONSTRUCTORS :

  function create ( i : integer ) return penta_double is

    res : constant penta_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( n : natural32 ) return penta_double is

    res : constant penta_double := create(double_float(n));

  begin
    return res;
  end create;

  function create ( n : natural64 ) return penta_double is

    res : constant penta_double := create(double_float(n));

  begin
    return res;
  end create;

  function create ( i : integer32 ) return penta_double is

    res : constant penta_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( i : integer64 ) return penta_double is

    res : constant penta_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( x : double_float ) return penta_double is

    res : constant penta_double := create(x,0.0,0.0,0.0,0.0);

  begin
    return res;
  end create;

  function create ( thumb,index,middle,ring,pink : double_float )
                  return penta_double is

    res : penta_double;

  begin
    res.thumb := thumb;
    res.index := index;
    res.middle := middle;
    res.ring := ring;
    res.pink := pink;
    return res;
  end create;

  function "abs" ( x : penta_double ) return penta_double is

    res : penta_double;

  begin
    if x.thumb < 0.0 then
      res.thumb := -x.thumb; res.index := -x.index;
      res.middle := -x.middle; res.ring := -x.ring; res.pink := -x.pink;
    else
      res.thumb := x.thumb; res.index := x.index;
      res.middle := x.middle; res.ring := x.ring; res.pink := x.pink;
    end if;
    return res;
  end "abs";

  function AbsVal ( x : penta_double ) return penta_double is

    res : penta_double;

  begin
    if x.thumb < 0.0 then
      res.thumb := -x.thumb; res.index := -x.index;
      res.middle := -x.middle; res.ring := -x.ring; res.pink := -x.pink;
    else
      res.thumb := x.thumb; res.index := x.index;
      res.middle := x.middle; res.ring := x.ring; res.pink := x.pink;
    end if;
    return res;
  end AbsVal;

  function floor ( x : penta_double ) return penta_double is

    res : penta_double;
    f0,f1,f2,f3,f4,f5 : double_float;

  begin
    f0 := double_float'floor(x.thumb);
    f1 := 0.0; f2 := 0.0; f3 := 0.0; f4 := 0.0; f5 := 0.0;
    if f0 = x.thumb then
      f1 := double_float'floor(x.index);
      if f1 = x.index then
        f2 := double_float'floor(x.middle);
        if f2 = x.middle then
          f3 := double_float'floor(x.ring);
          if f3 = x.ring then
            f4 := double_float'floor(x.pink);
          end if;
        end if;
      end if;
    end if;
    fast_renorm(f0,f1,f2,f3,f4,f5,
                res.thumb,res.index,res.middle,res.ring,res.pink);
    return res;
  end floor;

  function nint ( x : penta_double ) return penta_double is

    res : penta_double;
    x0,x1,x2,x3,x4,x5 : double_float;

  begin
    x1 := 0.0; x2 := 0.0; x3 := 0.0; x4 := 0.0; x5 := 0.0;
    x0 := Double_Double_Numbers.nint(x.thumb);
    if x0 = x.thumb then      -- first double is already an integer
      x1 := Double_Double_Numbers.nint(x.index);
      if x1 = x.index then    -- second double is already an integer
        x2 := Double_Double_Numbers.nint(x.ring);
        if x2 = x.ring then  -- third double is already an integer
          x3 := Double_Double_Numbers.nint(x.middle);
          if x3 = x.middle then  -- 4-th double is an integer
            x4 := Double_Double_Numbers.nint(x.pink);
          else
            if abs(x3 - x.ring) = 0.5 and x.pink < 0.0
             then x3 := x3 - 1.0;
            end if;
          end if;
        else
          if abs(x2 - x.middle) = 0.5 and x.ring < 0.0
           then x2 := x2 - 1.0;
          end if;
        end if;
      else
        if abs(x1 - x.index) = 0.5 and x.middle < 0.0
         then x1 := x1 - 1.0;
        end if;
      end if;
    else                     -- first double is not an integer
      if abs(x0 - x.thumb) = 0.5 and x.index < 0.0
       then x0 := x0 - 1.0;
      end if;
    end if;
    fast_renorm(x0,x1,x2,x3,x4,x5,
                res.thumb,res.index,res.middle,res.ring,res.pink);
    return res;
  end nint;


-- SELECTORS :

  function thumb_part ( x : penta_double ) return double_float is
  begin
    return x.thumb;
  end thumb_part;

  function index_part ( x : penta_double ) return double_float is
  begin
    return x.index;
  end index_part;

  function middle_part ( x : penta_double ) return double_float is
  begin
    return x.middle;
  end middle_part;

  function ring_part ( x : penta_double ) return double_float is
  begin
    return x.ring;
  end ring_part;

  function pink_part ( x : penta_double ) return double_float is
  begin
    return x.pink;
  end pink_part;

-- TYPE CASTS :

  function to_int ( x : penta_double ) return integer32 is
  begin
    return integer32(x.thumb);
  end to_int;

  function to_double ( x : penta_double ) return double_float is
  begin
    return x.thumb;
  end to_double;

-- COMPARISON and COPYING :

  function is_zero ( x : penta_double ) return boolean is
  begin
    return ((x.thumb = 0.0) and (x.index = 0.0) and
            (x.middle = 0.0) and (x.ring = 0.0) and (x.pink = 0.0));
  end is_zero;

  function is_one ( x : penta_double ) return boolean is
  begin
    return ((x.thumb = 1.0) and (x.index = 0.0) and
            (x.middle = 0.0) and (x.ring = 0.0) and (x.pink = 0.0));
  end is_one;

  function is_positive ( x : penta_double ) return boolean is
  begin
    return (x.thumb > 0.0);
  end is_positive;

  function is_negative ( x : penta_double ) return boolean is
  begin
    return (x.thumb < 0.0);
  end is_negative;

  function equal ( x,y : penta_double ) return boolean is
  begin
    return ((x.thumb = y.thumb) and (x.index = y.index) and
            (x.middle = y.middle) and (x.ring = y.ring) and
            (x.pink = y.pink));
  end equal;

  function equal ( x : penta_double; y : double_float ) return boolean is
  begin
    return ((x.thumb = y) and (x.index = 0.0) and
            (x.middle = 0.0) and (x.ring = 0.0) and (x.pink = 0.0));
  end equal;

  function "<" ( x,y : penta_double ) return boolean is
  begin
    return ((x.thumb < y.thumb)
         or (x.thumb = y.thumb and x.index < y.index)
         or (x.thumb = y.thumb and x.index = y.index and x.middle < y.middle)
         or (x.thumb = y.thumb and x.index = y.index and x.middle = y.middle
             and x.ring < y.ring)
         or (x.thumb = y.thumb and x.index = y.index and x.middle = y.middle
             and x.ring = y.ring and x.pink < y.pink));
  end "<";

  function "<" ( x : penta_double; y : double_float ) return boolean is
  begin
    return ((x.thumb < y)
         or (x.thumb = y and x.index < 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle < 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle = 0.0 and
             x.ring < 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle = 0.0 and
             x.ring = 0.0 and x.pink < 0.0));
  end "<";

  function "<" ( x : double_float; y : penta_double ) return boolean is
  begin
    return (y > x);
  end "<";

  function "<=" ( x,y : penta_double ) return boolean is
  begin
    return x < y or equal(x,y);
  end "<=";

  function "<=" ( x : penta_double; y : double_float ) return boolean is
  begin
    return ((x.thumb < y)
         or (x.thumb = y and x.index < 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle < 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle = 0.0 and
             x.ring < 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle = 0.0 and
             x.ring = 0.0 and x.pink <= 0.0));
  end "<=";

  function "<=" ( x : double_float; y : penta_double ) return boolean is
  begin
    return (y >= x);
  end "<=";

  function ">" ( x,y : penta_double ) return boolean is
  begin
    return ((x.thumb > y.thumb)
         or (x.thumb = y.thumb and x.index > y.index)
         or (x.thumb = y.thumb and x.index = y.index and x.middle > y.middle)
         or (x.thumb = y.thumb and x.index = y.index and x.middle = y.middle
             and x.ring > y.ring)
         or (x.thumb = y.thumb and x.index = y.index and x.middle = y.middle
             and x.ring = y.ring and x.pink > y.pink));
  end ">";

  function ">" ( x : penta_double; y : double_float ) return boolean is
  begin
    return ((x.thumb > y)
         or (x.thumb = y and x.index > 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle > 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle = 0.0 and
             x.ring > 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle = 0.0 and
             x.ring = 0.0 and x.pink > 0.0));
  end ">";

  function ">" ( x : double_float; y : penta_double ) return boolean is
  begin
    return (y < x);
  end ">";

  function ">=" ( x,y : penta_double ) return boolean is
  begin
    return x > y or equal(x,y);
  end ">=";

  function ">=" ( x : penta_double; y : double_float ) return boolean is
  begin
    return ((x.thumb > y)
         or (x.thumb = y and x.index > 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle > 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle = 0.0 and
             x.ring > 0.0)
         or (x.thumb = y and x.index = 0.0 and x.middle = 0.0 and
             x.ring = 0.0 and x.pink >= 0.0));
  end ">=";

  function ">=" ( x : double_float; y : penta_double ) return boolean is
  begin
    return (y <= x);
  end ">=";

  procedure copy ( x : in penta_double; y : in out penta_double ) is
  begin
    y.thumb := x.thumb; y.index := x.index;
    y.middle := x.middle; y.ring := x.ring; y.pink := x.pink;
  end copy;

-- ARITHMETICAL FUNCTIONS :

  function "+" ( x,y : penta_double ) return penta_double is

  -- ALGORITHM : baileyAdd_fast<5,5,5> generated by CAMPARY.

    res : penta_double;
    f0,f1,f2,f3,f4,f5,e : double_float;

  begin
    f5 := 0.0;
    Double_Double_Basics.two_sum(x.pink,y.pink,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(x.ring,y.ring,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(x.middle,y.middle,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(x.index,y.index,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(x.thumb,y.thumb,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,
                res.thumb,res.index,res.middle,res.ring,res.pink);
    return res;
  end "+";

  function "+" ( x : penta_double; y : double_float ) return penta_double is

  -- ALGORITHM : baileyAdd_fast<5,1,5> generated by CAMPARY.

    res : penta_double;

  begin
    renorm_add1(x.thumb,x.index,x.middle,x.ring,x.pink,y,
                res.thumb,res.index,res.middle,res.ring,res.pink);
    return res;
  end "+";

  function "+" ( x : double_float; y : penta_double ) return penta_double is

    res : constant penta_double := y + x;

  begin
    return res;
  end "+";

  function "+" ( x : penta_double ) return penta_double is

    res : penta_double;

  begin
    res.thumb := x.thumb;
    res.index := x.index;
    res.middle := x.middle;
    res.ring := x.ring;
    res.pink := x.pink;
    return res;
  end "+";

  function "-" ( x : penta_double ) return penta_double is

    res : penta_double;

  begin
    res.thumb := -x.thumb;
    res.index := -x.index;
    res.middle := -x.middle;
    res.ring := -x.ring;
    res.pink := -x.pink;
    return res;
  end "-";

  function "-" ( x,y : penta_double ) return penta_double is

    mny : constant penta_double := -y;
    res : constant penta_double := x + mny;

  begin
    return res;
  end "-";

  function "-" ( x : penta_double; y : double_float ) return penta_double is

    mny : constant double_float := -y;
    res : constant penta_double := x + mny;

  begin
    return res;
  end "-";

  function "-" ( x : double_float; y : penta_double ) return penta_double is

    mny : constant penta_double := -y;
    res : constant penta_double := x + mny;

  begin
    return res;
  end "-";

  function "*" ( x,y : penta_double ) return penta_double is

  -- ALGORITHM : baileyMul_fast<5,5,5> generated by CAMPARY.

    res : penta_double;
    f0,f1,f2,f3,f4,f5,p,e : double_float;

  begin
    f5 := x.index*y.pink;
    f5 := f5 + x.middle*y.ring;
    f5 := f5 + x.ring*y.middle;
    f5 := f5 + x.pink*y.index;
    Double_Double_Basics.two_prod(x.thumb,y.pink,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.index,y.ring,p,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.middle,y.middle,p,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.ring,y.index,p,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.pink,y.thumb,p,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.thumb,y.ring,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.index,y.middle,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.middle,y.index,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.ring,y.thumb,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.thumb,y.middle,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.index,y.index,p,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(f2,p,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.middle,y.thumb,p,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(f2,p,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.thumb,y.index,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.index,y.thumb,p,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_sum(f1,p,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.thumb,y.thumb,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,
                res.thumb,res.index,res.middle,res.ring,res.pink);
    return res;
  end "*";

  function "*" ( x : penta_double; y : double_float ) return penta_double is

  -- ALGORITHM : baileyMul_fast<5,1,5>

    res : penta_double;
    f0,f1,f2,f3,f4,f5,e : double_float;

  begin
    f5 := 0.0;
    Double_Double_Basics.two_prod(x.pink,y,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.ring,y,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.middle,y,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.index,y,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    Double_Double_Basics.two_prod(x.thumb,y,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    f5 := f5 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,
                res.thumb,res.index,res.middle,res.ring,res.pink);
    return res;
  end "*";

  function "*" ( x : double_float; y : penta_double ) return penta_double is

  -- DESCRIPTION :
  --   Returns x*y, the product of the double x with the penta double y.

    res : constant penta_double := y*x;

  begin
    return res;
  end "*";

  function Mul_pwr2 ( x : penta_double; y : double_float )
                    return penta_double is

    res : penta_double;

  begin
    res.thumb := x.thumb*y; res.index := x.index*y;
    res.middle := x.middle*y; res.ring := x.ring*y; res.pink := x.pink*y;
    return res;
  end Mul_pwr2;

  procedure Mul_pwr2 ( x : in out penta_double; y : in double_float ) is
  begin
    x.thumb := x.thumb*y; x.index := x.index*y;
    x.middle := x.middle*y; x.ring := x.ring*y; x.pink := x.pink*y;
  end Mul_pwr2;

  function "/" ( x,y : penta_double ) return penta_double is

    res,acc : penta_double;
    q0,q1,q2,q3,q4,q5 : double_float;

  begin
    q0 := x.thumb/y.thumb;   acc := q0*y; res := x - acc;
    q1 := res.thumb/y.thumb; acc := q1*y; res := res - acc;
    q2 := res.thumb/y.thumb; acc := q2*y; res := res - acc;
    q3 := res.thumb/y.thumb; acc := q3*y; res := res - acc;
    q4 := res.thumb/y.thumb; acc := q4*y; res := res - acc;
    q5 := res.thumb/y.thumb;
    fast_renorm(q0,q1,q2,q3,q4,q5,
                res.thumb,res.index,res.middle,res.ring,res.pink);
    return res;
  end "/";

  function "/" ( x : penta_double; y : double_float ) return penta_double is

    pdy : constant penta_double := create(y);
    res : constant penta_double := x/pdy;

  begin
    return res;
  end "/";

  function "/" ( x : double_float; y : penta_double ) return penta_double is

    pdx : constant penta_double := create(x);
    res : constant penta_double := pdx/y;

  begin
    return res;
  end "/";

  function sqr ( x : penta_double ) return penta_double is
  begin
    return x*x;
  end sqr;

  function "**" ( x : penta_double; n : integer ) return penta_double is

    res,acc : penta_double;
    absn : natural;

  begin
    if n = 0 then
      res.thumb := 1.0; res.index := 0.0; res.middle := 0.0;
      res.ring := 0.0; res.pink := 0.0;
    else
      if n > 0
       then absn := n;
       else absn := -n;
      end if;
      res.thumb := x.thumb; res.index := x.index;
      res.middle := x.middle; res.ring := x.ring; res.pink := x.pink;
      acc.thumb := 1.0; acc.index := 0.0; acc.middle := 0.0;
      acc.ring := 0.0; acc.pink := 0.0;
      if absn > 1 then          -- use binary exponentiation
        while absn > 0 loop
          if absn mod 2 = 1
           then acc := acc*res;
          end if;
          absn := absn/2;
          if absn > 0
           then res := res*res;
          end if;
        end loop;
      else
        acc.thumb := res.thumb; acc.index := res.index;
        acc.middle := res.middle; acc.ring := res.ring;
        acc.pink := res.pink;
      end if;
      if n < 0 then
        res := 1.0/acc;          -- compute reciprocal
      else 
        res.thumb := acc.thumb; res.index := acc.index;
        res.middle := acc.middle; res.ring := acc.ring;
        res.pink := acc.pink;
      end if;
    end if;
    return res;
  end "**";

  function ldexp ( x : penta_double; n : integer ) return penta_double is

    res : penta_double;
    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

  begin
    res.thumb := C_ldexp(x.thumb,n);
    res.index := C_ldexp(x.index,n);
    res.middle := C_ldexp(x.middle,n);
    res.ring := C_ldexp(x.ring,n);
    res.pink := C_ldexp(x.pink,n);
    return res;
  end ldexp;

  function "**" ( x,y : penta_double ) return penta_double is
  begin
    return exp(y*log(x));
  end "**";

  function "**" ( x : penta_double; y : double_float ) return penta_double is

    pd_y : constant penta_double := create(y);

  begin
    return x**pd_y;
  end "**";

  function exp ( x : penta_double ) return penta_double is

  -- Strategy:  We first reduce the size of x by noting that
  --    exp(kr + m * log(2)) = 2^m * exp(r)^k
  -- where m and k are integers.  By choosing m appropriately
  -- we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
  -- evaluated using the familiar Taylor series.  Reducing the
  -- argument substantially speeds up the convergence.

    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

    res : penta_double;
    k : constant double_float := C_ldexp(1.0,20); -- +4 in qd k
    inv_k : constant double_float := 1.0/k;
    e0 : constant double_float :=  2.71828182845904509E+00;
    e1 : constant double_float :=  1.44564689172925016E-16;
    e2 : constant double_float := -2.12771710803817676E-33;
    e3 : constant double_float :=  1.51563015984121914E-49;
    e4 : constant double_float := -9.33538137882084738E-66;
    exp1 : constant penta_double := Create(e0,e1,e2,e3,e4);
    L0 : constant double_float :=  6.93147180559945286E-01;
    L1 : constant double_float :=  2.31904681384629956E-17;
    L2 : constant double_float :=  5.70770843841621207E-34;
    L3 : constant double_float := -3.58243221060181142E-50;
    L4 : constant double_float := -1.35216967579886296E-66;
    log2 : constant penta_double := Create(L0,L1,L2,L3,L4);
    pd_eps : constant double_float := 6.747006683667535e-80; -- 2^-263
    tol : constant double_float := inv_k*pd_eps;
    m : constant double_float := double_float'floor(x.thumb/L0 + 0.5);
    i_fac : array(0..14) of penta_double;
      -- inverse factorials for Taylor expansion
    p,s,t : penta_double;
    cnt : integer;

  begin
    if x.thumb <= -709.0 then
      res := Create(0.0);
    elsif x.thumb >= 709.0 then
      res := Create(-1.0);
    elsif is_zero(x) then
      res := Create(1.0);
    elsif is_one(x) then
      res := exp1;
    else
      i_fac(0) := Create(1.66666666666666657E-01, 9.25185853854297066E-18,
                         5.13581318503262866E-34, 2.85094902409834186E-50,
                         1.58259462429329970E-66);
      i_fac(1) := Create(4.16666666666666644E-02, 2.31296463463574266E-18,
                         1.28395329625815716E-34, 7.12737256024585466E-51,
                         3.95648656073324926E-67);
      i_fac(2) := Create(8.33333333333333322E-03, 1.15648231731787138E-19,
                         1.60494162032269652E-36, 2.22730392507682967E-53,
                         3.09100512557285111E-70);
      i_fac(3) := Create(1.38888888888888894E-03, -5.30054395437357706E-20,
                        -1.73868675534958776E-36, -1.63335621172300840E-52,
                        -5.35774221765960817E-69);
      i_fac(4) := Create(1.98412698412698413E-04, 1.72095582934207053E-22,
                         1.49269123913941271E-40, 1.29470326746002471E-58,
                         1.12297607624339191E-76);
      i_fac(5) := Create(2.48015873015873016E-05, 2.15119478667758816E-23,
                         1.86586404892426588E-41, 1.61837908432503088E-59,
                         1.40372009530423989E-77);
      i_fac(6) := Create(2.75573192239858925E-06, -1.85839327404647208E-22,
                         8.49175460488199287E-39, -5.72661640789429621E-55,
                         2.61672391582886888E-71);
      i_fac(7) := Create(2.75573192239858883E-07, 2.37677146222502973E-23,
                        -3.26318890334088294E-40, 1.61435111860404415E-56,
                        -9.99798078191450073E-74);
      i_fac(8) := Create(2.50521083854417202E-08, -1.44881407093591197E-24,
                         2.04267351467144546E-41, -8.49632672007163175E-58,
                        -9.08907343810409139E-75);
      i_fac(9) := Create(2.08767569878681002E-09, -1.20734505911325997E-25,
                         1.70222792889287100E-42, 1.41609532150396700E-58,
                         5.13820161376913309E-75);
      i_fac(10) := Create(1.60590438368216133E-10, 1.25852945887520981E-26,
                         -5.31334602762985031E-43, 3.54021472597605528E-59,
                         -1.98567896059148884E-75);
      i_fac(11) := Create(1.14707455977297245E-11, 2.06555127528307454E-28,
                          6.88907923246664603E-45, 5.72920002655109095E-61,
                         -3.65551458930952487E-78);
      i_fac(12) := Create(7.64716373181981641E-13, 7.03872877733453001E-30,
                         -7.82753927716258345E-48, 1.92138649443790242E-64,
                         -1.43027453759388186E-80);
      i_fac(13) := Create(4.77947733238738525E-14, 4.39920548583408126E-31,
                         -4.89221204822661465E-49, 1.20086655902368901E-65,
                         -8.93921585996176165E-82);
      i_fac(14) := Create(2.81145725434552060E-15, 1.65088427308614326E-31,
                         -2.87777179307447918E-50, 4.27110689256293549E-67,
                         -2.93287743014724397E-83);
      res := mul_pwr2(x - m*log2,inv_k);
      p := res*res;
      s := res + mul_pwr2(p,0.5);
      cnt := 0;
      loop
        p := p*res;
        t := i_fac(cnt)*p;
        s := s + t;
        exit when abs(t.thumb) <= tol;
        cnt := cnt + 1;
        exit when (cnt >= 15); -- use all of the expansion ...
      end loop;
      for i in 1..20 loop -- 20 times s = mul_pwr2(s,2.0) + sqr(s);
        p := Mul_pwr2(s,2.0);
        t := s*s;
        s := p + t;
      end loop;
      res := s + 1.0;
      cnt := integer(m);
      res := ldexp(res,cnt);
    end if;
    return res;
  end exp;

  function log ( x : penta_double ) return penta_double is

  -- Strategy.  The Taylor series for log converges much more
  --   slowly than that of exp, due to the lack of the factorial
  --   term in the denominator.  Hence this routine instead tries
  --   to determine the root of the function f(x) = exp(x) - a
  --   using Newton iteration.  The iteration is given by
  --       x' = x - f(x)/f'(x)
  --          = x - (1 - a * exp(-x))
  --          = x + a * exp(-x) - 1.
  --   Four iterations are needed, since Newton's iteration
  --   approximately doubles the number of digits per iteration.

    res,acc : penta_double;

  begin
    if is_one(x) then
      res := Create(0.0);
    elsif x.thumb <= 0.0 then
      put_line("td_log: argument is not positive");
      res := Create(-1.0);
    else
      res := Create(Standard_Mathematical_Functions.LN(x.thumb));
      for i in 1..4 loop       -- r = r + x*exp(-r) - 1
        acc := x*exp(-res);    -- acc = x*exp(-res)      
        res := res + acc;      -- res = res + x*exp(-res)
        res := res - 1.0;
      end loop;
    end if;
    return res;
  end log;

  function log10 ( x : penta_double ) return penta_double is

    res : penta_double;
    log10_0 : constant double_float :=  2.30258509299404590E+00;
    log10_1 : constant double_float := -2.17075622338224935E-16;
    log10_2 : constant double_float := -9.98426245446577657E-33;
    log10_3 : constant double_float := -4.02335745445020638E-49;
    log10_4 : constant double_float :=  1.92889952896933719E-65;
    logten : constant penta_double
           := create(log10_0,log10_1,log10_2,log10_3,log10_4);

  begin
    res := log(x)/logten;
    return res;
  end log10;

-- ARITHMETICAL OPERATIONS AS PROCEDURES :

  procedure Add ( x : in out penta_double; y : in penta_double ) is
  begin
    x := x + y;
  end Add;

  procedure Sub ( x : in out penta_double; y : in penta_double ) is
  begin
    x := x - y;
  end Sub;

  procedure Min ( x : in out penta_double ) is
  begin
    x := -x;
  end Min;

  procedure Mul ( x : in out penta_double; y : in penta_double ) is
  begin
    x := x*y;
  end Mul;

  procedure Div ( x : in out penta_double; y : in penta_double ) is
  begin
    x := x/y;
  end Div;

-- DESTRUCTOR :

  procedure clear ( x : in out penta_double ) is
  begin
    x.thumb := 0.0; x.index := 0.0;
    x.middle := 0.0; x.ring := 0.0; x.pink := 0.0;
  end clear;

end Penta_Double_Numbers;
