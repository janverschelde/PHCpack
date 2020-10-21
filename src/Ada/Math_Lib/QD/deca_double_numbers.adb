with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with Double_Double_Basics;
with Fast_Double_Renormalizations;       use Fast_Double_Renormalizations;

package body Deca_Double_Numbers is

-- CONSTRUCTORS :

  function create ( i : integer ) return deca_double is

    res : constant deca_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( n : natural32 ) return deca_double is

    res : constant deca_double := create(double_float(n));

  begin
    return res;
  end create;

  function create ( n : natural64 ) return deca_double is

    res : constant deca_double := create(double_float(n));

  begin
    return res;
  end create;

  function create ( i : integer32 ) return deca_double is

    res : constant deca_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( i : integer64 ) return deca_double is

    res : constant deca_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( x : double_float ) return deca_double is

    res : constant deca_double
        := create(x,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

  begin
    return res;
  end create;

  function create ( right_thumb,right_index,right_middle : double_float;
                    right_ring,right_pink : double_float;
                    left_thumb,left_index,left_middle : double_float;
                    left_ring,left_pink : double_float ) return deca_double is

    res : deca_double;

  begin
    res.right_thumb := right_thumb;
    res.right_index := right_index;
    res.right_middle := right_middle; 
    res.right_ring := right_ring;
    res.right_pink := right_pink;
    res.left_thumb := left_thumb;
    res.left_index := left_index;
    res.left_middle := left_middle;
    res.left_ring := left_ring;
    res.left_pink := left_pink;
    return res;
  end create;

  function "abs" ( x : deca_double ) return deca_double is

    res : deca_double;

  begin
    if x.right_thumb < 0.0 then
      res.right_thumb := -x.right_thumb; res.right_index := -x.right_index;
      res.right_middle := -x.right_middle; res.right_ring := -x.right_ring;
      res.right_pink := -x.right_pink;
      res.left_thumb := -x.left_thumb; res.left_index := -x.left_index;
      res.left_middle := -x.left_middle; res.left_ring := -x.left_ring;
      res.left_pink := -x.left_pink;
    else
      res.right_thumb := x.right_thumb; res.right_index := x.right_index;
      res.right_middle := x.right_middle; res.right_ring := x.right_ring;
      res.right_pink := x.right_pink;
      res.left_thumb := x.left_thumb; res.left_index := x.left_index;
      res.left_middle := x.left_middle; res.left_ring := x.left_ring;
      res.left_pink := x.left_pink;
    end if;
    return res;
  end "abs";

  function AbsVal ( x : deca_double ) return deca_double is

    res : deca_double;

  begin
    if x.right_thumb < 0.0 then
      res.right_thumb := -x.right_thumb; res.right_index := -x.right_index;
      res.right_middle := -x.right_middle; res.right_ring := -x.right_ring;
      res.right_pink := -x.right_pink;
      res.left_thumb := -x.left_thumb; res.left_index := -x.left_index;
      res.left_middle := -x.left_middle; res.left_ring := -x.left_ring;
      res.left_pink := -x.left_pink;
    else
      res.right_thumb := x.right_thumb; res.right_index := x.right_index;
      res.right_middle := x.right_middle; res.right_ring := x.right_ring;
      res.right_pink := x.right_pink;
      res.left_thumb := x.left_thumb; res.left_index := x.left_index;
      res.left_middle := x.left_middle; res.left_ring := x.left_ring;
      res.left_pink := x.left_pink;
    end if;
    return res;
  end AbsVal;

  function floor ( x : deca_double ) return deca_double is

    res : deca_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10 : double_float;

  begin
    f0 := double_float'floor(x.right_thumb);
    f1 := 0.0; f2 := 0.0; f3 := 0.0; f4 := 0.0; f5 := 0.0;
    f6 := 0.0; f7 := 0.0; f8 := 0.0; f9 := 0.0; f10 := 0.0;
    if f0 = x.right_thumb then
      f1 := double_float'floor(x.right_index);
      if f1 = x.right_index then
        f2 := double_float'floor(x.right_middle);
        if f2 = x.right_middle then
          f3 := double_float'floor(x.right_ring);
          if f3 = x.right_ring then
            f4 := double_float'floor(x.right_pink);
            if f4 = x.right_pink then
              f5 := double_float'floor(x.left_thumb);
              if f5 = x.left_thumb then
                f6 := double_float'floor(x.left_index);
                if f6 = x.left_index then
                  f7 := double_float'floor(x.left_middle);
                  if f7 = x.left_middle then
                    f8 := double_float'floor(x.left_ring);
                    if f8 = x.left_ring then
                      f9 := double_float'floor(x.left_pink);
                    end if;
                  end if;
                end if;
              end if;
            end if;
          end if;
        end if;
      end if;
    end if;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                res.right_thumb,res.right_index,res.right_middle,
                res.right_ring,res.right_pink,
                res.left_thumb,res.left_index,res.left_middle,
                res.left_ring,res.left_pink);
    return res;
  end floor;

  function nint ( x : deca_double ) return deca_double is

    res : deca_double;
    x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 : double_float;

  begin
    x1 := 0.0; x2 := 0.0; x3 := 0.0; x4 := 0.0; x5 := 0.0;
    x6 := 0.0; x7 := 0.0; x8 := 0.0; x9 := 0.0; x10 := 0.0;
    x0 := Double_Double_Numbers.nint(x.right_thumb);
    if x0 = x.right_thumb then      -- first double is already an integer
      x1 := Double_Double_Numbers.nint(x.right_index);
      if x1 = x.right_index then    -- second double is already an integer
        x2 := Double_Double_Numbers.nint(x.right_middle);
        if x2 = x.right_middle then  -- third double is already an integer
          x3 := Double_Double_Numbers.nint(x.right_ring);
          if x3 = x.right_ring then  -- 4-th double is an integer
            x4 := Double_Double_Numbers.nint(x.right_pink);
            if x4 = x.right_pink then  -- 5-th double is an integer
              x5 := Double_Double_Numbers.nint(x.left_thumb);
              if x5 = x.left_thumb then  -- 6-th double is an integer
                x6 := Double_Double_Numbers.nint(x.left_index);
                if x6 = x.left_index then  -- 7-th double is an integer
                  x7 := Double_Double_Numbers.nint(x.left_middle);
                  if x7 = x.left_middle then  -- 8-th double is an integer
                    x8 := Double_Double_Numbers.nint(x.left_ring);
                    if x8 = x.left_ring then  -- 9-th double is an integer
                      x9 := Double_Double_Numbers.nint(x.left_ring);
                    else
                      if abs(x8 - x.left_ring) = 0.5 and x.left_pink < 0.0
                       then x8 := x8 - 1.0;
                      end if;
                    end if;
                  else
                    if abs(x7 - x.left_middle) = 0.5 and x.left_ring < 0.0
                     then x7 := x7 - 1.0;
                    end if;
                  end if;
                else
                  if abs(x6 - x.left_index) = 0.5 and x.left_middle < 0.0
                   then x6 := x6 - 1.0;
                  end if;
                end if;
              else
                if abs(x5 - x.left_thumb) = 0.5 and x.left_index < 0.0
                 then x5 := x5 - 1.0;
                end if;
              end if;
            else
              if abs(x4 - x.right_pink) = 0.5 and x.left_thumb < 0.0
               then x4 := x4 - 1.0;
              end if;
            end if;
          else
            if abs(x3 - x.right_ring) = 0.5 and x.right_pink < 0.0
             then x3 := x3 - 1.0;
            end if;
          end if;
        else
          if abs(x2 - x.right_middle) = 0.5 and x.right_ring < 0.0
           then x2 := x2 - 1.0;
          end if;
        end if;
      else
        if abs(x1 - x.right_index) = 0.5 and x.right_middle < 0.0
         then x1 := x1 - 1.0;
        end if;
      end if;
    else                     -- first double is not an integer
      if abs(x0 - x.right_thumb) = 0.5 and x.right_index < 0.0
       then x0 := x0 - 1.0;
      end if;
    end if;
    fast_renorm(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,
                res.right_thumb,res.right_index,res.right_middle,
                res.right_ring,res.right_pink,
                res.left_thumb,res.left_index,res.left_middle,
                res.left_ring,res.left_pink);
    return res;
  end nint;

-- SELECTORS :

  function thumb_right ( x : deca_double ) return double_float is
  begin
    return x.right_thumb;
  end thumb_right;

  function index_right ( x : deca_double ) return double_float is
  begin
    return x.right_index;
  end index_right;

  function middle_right ( x : deca_double ) return double_float is
  begin
    return x.right_middle;
  end middle_right;

  function ring_right ( x : deca_double ) return double_float is
  begin
    return x.right_ring;
  end ring_right;

  function pink_right ( x : deca_double ) return double_float is
  begin
    return x.right_pink;
  end pink_right;

  function thumb_left ( x : deca_double ) return double_float is
  begin
    return x.left_thumb;
  end thumb_left;

  function index_left ( x : deca_double ) return double_float is
  begin
    return x.left_index;
  end index_left;

  function middle_left ( x : deca_double ) return double_float is
  begin
    return x.left_middle;
  end middle_left;

  function ring_left ( x : deca_double ) return double_float is
  begin
    return x.left_ring;
  end ring_left;

  function pink_left ( x : deca_double ) return double_float is
  begin
    return x.left_pink;
  end pink_left;

-- TYPE CASTS :

  function to_int ( x : deca_double ) return integer32 is
  begin
    return integer32(x.right_thumb);
  end to_int;

  function to_double ( x : deca_double ) return double_float is
  begin
    return x.right_thumb;
  end to_double;

  function to_double_double ( x : deca_double ) return double_double is

    res : constant double_double := create(x.right_thumb,x.right_index);

  begin
    return res;
  end to_double_double;

  function to_triple_double ( x : deca_double ) return triple_double is

    res : constant triple_double
        := create(x.right_thumb,x.right_index,x.right_middle);

  begin
    return res;
  end to_triple_double;

  function to_quad_double ( x : deca_double ) return quad_double is

    res : constant quad_double
        := create(x.right_thumb,x.right_index,x.right_middle,x.right_ring);

  begin
    return res;
  end to_quad_double;

  function to_penta_double ( x : deca_double ) return penta_double is

    res : constant penta_double
        := create(x.right_thumb,x.right_index,x.right_middle,
                  x.right_ring,x.right_pink);

  begin
    return res;
  end to_penta_double;

  function to_octo_double ( x : deca_double ) return octo_double is

    res : constant octo_double
        := create(x.right_thumb,x.right_index,x.right_middle,x.right_ring,
                  x.right_pink,x.left_thumb,x.left_index,x.left_middle);

  begin
    return res;
  end to_octo_double;

-- COMPARISON and COPYING :

  function is_zero ( x : deca_double ) return boolean is
  begin
    return ((x.right_thumb = 0.0) and
            (x.right_index = 0.0) and
            (x.right_middle = 0.0) and
            (x.right_ring = 0.0) and
            (x.right_pink = 0.0) and
            (x.left_thumb = 0.0) and
            (x.left_index = 0.0) and
            (x.left_middle = 0.0) and
            (x.left_ring = 0.0) and
            (x.left_pink = 0.0));
  end is_zero;

  function is_one ( x : deca_double ) return boolean is
  begin
    return ((x.right_thumb = 1.0) and
            (x.right_index = 0.0) and
            (x.right_middle = 0.0) and
            (x.right_ring = 0.0) and
            (x.right_pink = 0.0) and
            (x.left_thumb = 0.0) and
            (x.left_index = 0.0) and
            (x.left_middle = 0.0) and
            (x.left_ring = 0.0) and
            (x.left_pink = 0.0));
  end is_one;

  function is_positive ( x : deca_double ) return boolean is
  begin
    return (x.right_thumb > 0.0);
  end is_positive;

  function is_negative ( x : deca_double ) return boolean is
  begin
    return (x.right_thumb < 0.0);
  end is_negative;

  function equal ( x,y : deca_double ) return boolean is
  begin
    return ((x.right_thumb = y.right_thumb) and
            (x.right_index = y.right_index) and
            (x.right_middle = y.right_middle) and
            (x.right_ring = y.right_ring) and
            (x.right_pink = y.right_pink) and
            (x.left_thumb = y.left_thumb) and
            (x.left_index = y.left_index) and
            (x.left_middle = y.left_middle) and
            (x.left_ring = y.left_ring) and
            (x.left_pink = y.left_pink));
  end equal;

  function equal ( x : deca_double; y : double_float ) return boolean is
  begin
    return ((x.right_thumb = y) and
            (x.right_index = 0.0) and
            (x.right_middle = 0.0) and
            (x.right_ring = 0.0) and
            (x.right_pink = 0.0) and
            (x.left_thumb = 0.0) and
            (x.left_index = 0.0) and
            (x.left_middle = 0.0) and
            (x.left_ring = 0.0) and
            (x.left_pink = 0.0));
  end equal;

  function "<" ( x,y : deca_double ) return boolean is
  begin
    return ((x.right_thumb < y.right_thumb)
         or (x.right_thumb = y.right_thumb and x.right_index < y.right_index)
         or (x.right_thumb = y.right_thumb and x.right_index = y.right_index
             and x.right_middle < y.right_middle)
         or (x.right_thumb = y.right_thumb and x.right_index = y.right_index
             and x.right_middle = y.right_middle
             and x.right_ring < y.right_ring)
         or (x.right_thumb = y.right_thumb and x.right_index = y.right_index
             and x.right_middle = y.right_middle
             and x.right_ring = y.right_ring and x.right_pink < y.right_pink)
         or (x.right_thumb = y.right_thumb and x.right_index = y.right_index
             and x.right_middle = y.right_middle
             and x.right_ring = y.right_ring and x.right_pink = y.right_pink
             and x.left_thumb < y.left_thumb)
         or (x.right_thumb = y.right_thumb and x.right_index = y.right_index
             and x.right_middle = y.right_middle
             and x.right_ring = y.right_ring and x.right_pink = y.right_pink
             and x.left_thumb = y.left_thumb and x.left_index < y.left_index)
         or (x.right_thumb = y.right_thumb and x.right_index = y.right_index
             and x.right_middle = y.right_middle
             and x.right_ring = y.right_ring and x.right_pink = y.right_pink
             and x.left_thumb = y.left_thumb and x.left_index = y.left_index
             and x.left_middle < y.left_middle)
         or (x.right_thumb = y.right_thumb and x.right_index = y.right_index
             and x.right_middle = y.right_middle
             and x.right_ring = y.right_ring and x.right_pink = y.right_pink
             and x.left_thumb = y.left_thumb and x.left_index = y.left_index
             and x.left_middle = y.left_middle
             and x.left_ring < y.left_ring)
         or (x.right_thumb = y.right_thumb and x.right_index = y.right_index
             and x.right_middle = y.right_middle
             and x.right_ring = y.right_ring and x.right_pink = y.right_pink
             and x.left_thumb = y.left_thumb and x.left_index = y.left_index
             and x.left_middle = y.left_middle
             and x.left_ring = y.left_ring and x.left_pink = y.left_pink));
  end "<";

  function "<" ( x : deca_double; y : double_float ) return boolean is
  begin
   return ((x.right_thumb < y)
        or (x.right_thumb = y and x.right_index < 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle < 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring < 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink < 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink = 0.0 and x.left_thumb < 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink = 0.0 and x.left_thumb = 0.0 and
            x.left_index < 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink = 0.0 and x.left_thumb = 0.0 and
            x.left_index = 0.0 and x.left_middle < 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink = 0.0 and x.left_thumb = 0.0 and
            x.left_index = 0.0 and x.left_middle = 0.0 and
            x.left_ring < 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink = 0.0 and x.left_thumb = 0.0 and
            x.left_index = 0.0 and x.left_middle = 0.0 and
            x.left_ring = 0.0 and x.left_pink < 0.0));
  end "<";

  function "<" ( x : double_float; y : deca_double ) return boolean is
  begin
    return (y > x);
  end "<";

  function "<=" ( x,y : deca_double ) return boolean is
  begin
    return (x < y) or equal(x,y);
  end "<=";

  function "<=" ( x : deca_double; y : double_float ) return boolean is
  begin
    return (x < y) or equal(x,y);
  end "<=";

  function "<=" ( x : double_float; y : deca_double ) return boolean is
  begin
    return (x < y) or equal(y,x);
  end "<=";

  function ">" ( x,y : deca_double ) return boolean is
  begin
    return (y < x);
  end ">";

  function ">" ( x : deca_double; y : double_float ) return boolean is
  begin
   return ((x.right_thumb > y)
        or (x.right_thumb = y and x.right_index > 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle > 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring > 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink > 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink = 0.0 and x.left_thumb > 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink = 0.0 and x.left_thumb = 0.0 and
            x.left_index > 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink = 0.0 and x.left_thumb = 0.0 and
            x.left_index = 0.0 and x.left_middle > 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink = 0.0 and x.left_thumb = 0.0 and
            x.left_index = 0.0 and x.left_middle = 0.0 and
            x.left_ring > 0.0)
        or (x.right_thumb = y and x.right_index = 0.0 and
            x.right_middle = 0.0 and x.right_ring = 0.0 and
            x.right_pink = 0.0 and x.left_thumb = 0.0 and
            x.left_index = 0.0 and x.left_middle = 0.0 and
            x.left_ring = 0.0 and x.left_pink > 0.0));
  end ">";

  function ">" ( x : double_float; y : deca_double ) return boolean is
  begin
    return (y < x);
  end ">";

  function ">=" ( x,y : deca_double ) return boolean is
  begin
    return (x > y) or equal(x,y);
  end ">=";

  function ">=" ( x : deca_double; y : double_float ) return boolean is
  begin
    return (x > y) or equal(x,y);
  end ">=";

  function ">=" ( x : double_float; y : deca_double ) return boolean is
  begin
    return (x > y) or equal(y,x);
  end ">=";

  procedure copy ( x : in deca_double; y : in out deca_double ) is
  begin
    y.right_thumb := x.right_thumb;
    y.right_index := x.right_index;
    y.right_middle := x.right_middle;
    y.right_ring := x.right_ring;
    y.right_pink := x.right_pink;
    y.left_thumb := x.left_thumb;
    y.left_index := x.left_index;
    y.left_middle := x.left_middle;
    y.left_ring := x.left_ring;
    y.left_pink := x.left_pink;
  end copy;

-- ARITHMETICAL FUNCTIONS :

  function "+" ( x,y : deca_double ) return deca_double is

  -- ALGORITHM : baileyAdd_fast<10,10,10> generated by CAMPARY.
   
    res : deca_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,e : double_float;

  begin
    f10 := 0.0;
    Double_Double_Basics.two_sum(x.left_pink,y.left_pink,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.left_ring,y.left_ring,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.left_middle,y.left_middle,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.left_index,y.left_index,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.left_thumb,y.left_thumb,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.right_pink,y.right_pink,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.right_ring,y.right_ring,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.right_middle,y.right_middle,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.right_index,y.right_index,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.right_thumb,y.right_thumb,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                res.right_thumb,res.right_index,res.right_middle,
                res.right_ring,res.right_pink,
                res.left_thumb,res.left_index,res.left_middle,
                res.left_ring,res.left_pink);
    return res;
  end "+";

  function "+" ( x : deca_double; y : double_float ) return deca_double is

    res : deca_double;

  begin
    renorm_add1(x.right_thumb,x.right_index,x.right_middle,
                x.right_ring,x.right_pink,
                x.left_thumb,x.left_index,x.left_middle,
                x.left_ring,x.left_pink,y,
                res.right_thumb,res.right_index,res.right_middle,
                res.right_ring,res.right_pink,
                res.left_thumb,res.left_index,res.left_middle,
                res.left_ring,res.left_pink);
    return res;
  end "+";

  function "+" ( x : double_float; y : deca_double ) return deca_double is

    res : constant deca_double := y + x;

  begin
    return res;
  end "+";

  function "+" ( x : deca_double ) return deca_double is

    res : deca_double;

  begin
    res.right_thumb := x.right_thumb;
    res.right_index := x.right_index;
    res.right_middle := x.right_middle;
    res.right_ring := x.right_ring;
    res.right_pink := x.right_pink;
    res.left_thumb := x.left_thumb;
    res.left_index := x.left_index;
    res.left_middle := x.left_middle;
    res.left_ring := x.left_ring;
    res.left_pink := x.left_pink;
    return res;
  end "+";

  function "-" ( x : deca_double ) return deca_double is

    res : deca_double;

  begin
    res.right_thumb := -x.right_thumb;
    res.right_index := -x.right_index;
    res.right_middle := -x.right_middle;
    res.right_ring := -x.right_ring;
    res.right_pink := -x.right_pink;
    res.left_thumb := -x.left_thumb;
    res.left_index := -x.left_index;
    res.left_middle := -x.left_middle;
    res.left_ring := -x.left_ring;
    res.left_pink := -x.left_pink;
    return res;
  end "-";

  function "-" ( x,y : deca_double ) return deca_double is

    mny : constant deca_double := -y;
    res : constant deca_double := x + mny;

  begin
    return res;
  end "-";

  function "-" ( x : deca_double; y : double_float ) return deca_double is

    mny : constant double_float := -y;
    res : constant deca_double := x + mny;

  begin
    return res;
  end "-";

  function "-" ( x : double_float; y : deca_double ) return deca_double is

    mny : constant deca_double := -y;
    res : constant deca_double := x + mny;

  begin
    return res;
  end "-";

  function "*" ( x,y : deca_double ) return deca_double is

  -- ALGORITHM :baileyMul_fast<10,10,10> generated by CAMPARY.

    res : deca_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,p,e : double_float;

  begin
    f10 := x.right_index*y.left_pink;
    f10 := f10 + x.right_middle*y.left_ring;
    f10 := f10 + x.right_ring*y.left_middle;
    f10 := f10 + x.right_pink*y.left_index;
    f10 := f10 + x.left_thumb*y.left_thumb;
    f10 := f10 + x.left_index*y.right_pink;
    f10 := f10 + x.left_middle*y.right_ring;
    f10 := f10 + x.left_ring*y.right_middle;
    f10 := f10 + x.left_pink*y.right_index;
    Double_Double_Basics.two_prod(x.right_thumb,y.left_pink,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_index,y.left_ring,p,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_middle,y.left_middle,p,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_ring,y.left_index,p,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_pink,y.left_thumb,p,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_thumb,y.right_pink,p,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_index,y.right_ring,p,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_middle,y.right_middle,p,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_ring,y.right_index,p,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_pink,y.right_thumb,p,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_thumb,y.left_ring,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_index,y.left_middle,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_middle,y.left_index,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_ring,y.left_thumb,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_pink,y.right_pink,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_thumb,y.right_ring,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_index,y.right_middle,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_middle,y.right_index,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_ring,y.right_thumb,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_thumb,y.left_middle,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_index,y.left_index,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_middle,y.left_thumb,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_ring,y.right_pink,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_pink,y.right_ring,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_thumb,y.right_middle,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_index,y.right_index,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_middle,y.right_thumb,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_thumb,y.left_index,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_index,y.left_thumb,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_middle,y.right_pink,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_ring,y.right_ring,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_pink,y.right_middle,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_thumb,y.right_index,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_index,y.right_thumb,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_thumb,y.left_thumb,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_index,y.right_pink,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_middle,y.right_ring,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_ring,y.right_middle,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_pink,y.right_index,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_thumb,y.right_thumb,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_thumb,y.right_pink,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_index,y.right_ring,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_middle,y.right_middle,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_ring,y.right_index,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_pink,y.right_thumb,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_thumb,y.right_ring,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_index,y.right_middle,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_middle,y.right_index,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_ring,y.right_thumb,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_thumb,y.right_middle,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_index,y.right_index,p,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f2,p,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_middle,y.right_thumb,p,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f2,p,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_thumb,y.right_index,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_index,y.right_thumb,p,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(f1,p,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_thumb,y.right_thumb,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                res.right_thumb,res.right_index,res.right_middle,
                res.right_ring,res.right_pink,
                res.left_thumb,res.left_index,res.left_middle,
                res.left_ring,res.left_pink);
    return res;
  end "*";

  function "*" ( x : deca_double; y : double_float ) return deca_double is

  -- ALGORITHM : baileyMul_fast<10,1,10> generated by CAMPARY.

    res : deca_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,e : double_float;

  begin
    f10 := 0.0;
    Double_Double_Basics.two_prod(x.left_pink,y,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_ring,y,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_middle,y,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_index,y,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.left_thumb,y,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_pink,y,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_ring,y,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_middle,y,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_index,y,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_prod(x.right_thumb,y,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                res.right_thumb,res.right_index,res.right_middle,
                res.right_ring,res.right_pink,
                res.left_thumb,res.left_index,res.left_middle,
                res.left_ring,res.left_pink);
    return res;
  end "*";

  function "*" ( x : double_float; y : deca_double ) return deca_double is

    res : constant deca_double := y*x;

  begin
    return res;
  end "*";

  function Mul_pwr2 ( x : deca_double; y : double_float )
                    return deca_double is

    res : deca_double;

  begin
    res.right_thumb := y*x.right_thumb;
    res.right_index := y*x.right_index;
    res.right_middle := y*x.right_middle;
    res.right_ring := y*x.right_ring;
    res.right_pink := y*x.right_pink;
    res.left_thumb := y*x.left_thumb;
    res.left_index := y*x.left_index;
    res.left_middle := y*x.left_middle;
    res.left_ring := y*x.left_ring;
    res.left_pink := y*x.left_pink;
    return res;
  end Mul_pwr2;

  procedure Mul_pwr2 ( x : in out deca_double; y : in double_float ) is
  begin
    x.right_thumb := y*x.right_thumb;
    x.right_index := y*x.right_index;
    x.right_middle := y*x.right_middle;
    x.right_ring := y*x.right_ring;
    x.right_pink := y*x.right_pink;
    x.left_thumb := y*x.left_thumb;
    x.left_index := y*x.left_index;
    x.left_middle := y*x.left_middle;
    x.left_ring := y*x.left_ring;
    x.left_pink := y*x.left_pink;
  end Mul_pwr2;

  function "/" ( x,y : deca_double ) return deca_double is

    res,acc : deca_double;
    q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10 : double_float;

  begin
    q0 := x.right_thumb/y.right_thumb;   acc := q0*y; res := x - acc;
    q1 := res.right_thumb/y.right_thumb; acc := q1*y; res := res - acc;
    q2 := res.right_thumb/y.right_thumb; acc := q2*y; res := res - acc;
    q3 := res.right_thumb/y.right_thumb; acc := q3*y; res := res - acc;
    q4 := res.right_thumb/y.right_thumb; acc := q4*y; res := res - acc;
    q5 := res.right_thumb/y.right_thumb; acc := q5*y; res := res - acc;
    q6 := res.right_thumb/y.right_thumb; acc := q6*y; res := res - acc;
    q7 := res.right_thumb/y.right_thumb; acc := q7*y; res := res - acc;
    q8 := res.right_thumb/y.right_thumb; acc := q8*y; res := res - acc;
    q9 := res.right_thumb/y.right_thumb; acc := q9*y; res := res - acc;
    q10 := res.right_thumb/y.right_thumb;
    fast_renorm(q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,
                res.right_thumb,res.right_index,res.right_middle,
                res.right_ring,res.right_pink,
                res.left_thumb,res.left_index,res.left_middle,
                res.left_ring,res.left_pink);
    return res;
  end "/";

  function "/" ( x : deca_double; y : double_float ) return deca_double is

    ddy : constant deca_double := create(y);
    res : constant deca_double := x/ddy;

  begin
    return res;
  end "/";

  function "/" ( x : double_float; y : deca_double ) return deca_double is

    ddx : constant deca_double := create(x);
    res : constant deca_double := ddx/y;

  begin
    return res;
  end "/";

  function sqr ( x : deca_double ) return deca_double is
  begin
    return x*x;
  end sqr;

  function "**" ( x : deca_double; n : integer ) return deca_double is

    res,acc : deca_double;
    absn : natural;

  begin
    if n = 0 then
      res.right_thumb := 1.0; res.right_index := 0.0;
      res.right_middle := 0.0; res.right_ring := 0.0; res.right_pink := 0.0;
      res.left_thumb := 0.0; res.left_index := 0.0;
      res.left_middle := 0.0; res.left_ring := 0.0; res.left_pink := 0.0;
    else
      if n > 0
       then absn := n;
       else absn := -n;
      end if;
      res.right_thumb := x.right_thumb; res.right_index := x.right_index;
      res.right_middle := x.right_middle; res.right_ring := x.right_ring;
      res.right_pink := x.right_pink;
      res.left_thumb := x.left_thumb; res.left_index := x.left_index;
      res.left_middle := x.left_middle; res.left_ring := x.left_ring;
      res.left_pink := x.left_pink;
      acc.right_thumb := 1.0; acc.right_index := 0.0;
      acc.right_middle := 0.0; acc.right_ring := 0.0; acc.right_pink := 0.0;
      acc.left_thumb := 0.0; acc.left_index := 0.0;
      acc.left_middle := 0.0; acc.left_ring := 0.0; acc.left_pink := 0.0;
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
        acc.right_thumb := res.right_thumb;
        acc.right_index := res.right_index;
        acc.right_middle := res.right_middle;
        acc.right_ring := res.right_ring;
        acc.right_pink := res.right_pink;
        acc.left_thumb := res.left_thumb;
        acc.left_index := res.left_index;
        acc.left_middle := res.left_middle;
        acc.left_ring := res.left_ring;
        acc.left_pink := res.left_pink;
      end if;
      if n < 0 then
        res := 1.0/acc;          -- compute reciprocal
      else 
        res.right_thumb := acc.right_thumb;
        res.right_index := acc.right_index;
        res.right_middle := acc.right_middle;
        res.right_ring := acc.right_ring;
        res.right_pink := acc.right_pink;
        res.left_thumb := acc.left_thumb;
        res.left_index := acc.left_index;
        res.left_middle := acc.left_middle;
        res.left_ring := acc.left_ring;
        res.left_pink := acc.left_pink;
      end if;
    end if;
    return res;
  end "**";

  function ldexp ( x : deca_double; n : integer ) return deca_double is

    res : deca_double;
    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

  begin
    res.right_thumb := C_ldexp(x.right_thumb,n);
    res.right_index := C_ldexp(x.right_index,n);
    res.right_middle := C_ldexp(x.right_middle,n);
    res.right_ring := C_ldexp(x.right_ring,n);
    res.right_pink := C_ldexp(x.right_pink,n);
    res.left_thumb := C_ldexp(x.left_thumb,n);
    res.left_index := C_ldexp(x.left_index,n);
    res.left_middle := C_ldexp(x.left_middle,n);
    res.left_ring := C_ldexp(x.left_ring,n);
    res.left_pink := C_ldexp(x.left_pink,n);
    return res;
  end ldexp;

  function "**" ( x,y : deca_double ) return deca_double is
  begin
    return exp(y*log(x));
  end "**";

  function "**" ( x : deca_double; y : double_float ) return deca_double is

    da_y : constant deca_double := create(y);

  begin
    return x**da_y;
  end "**";

  function exp ( x : deca_double ) return deca_double is

  -- Strategy:  We first reduce the size of x by noting that
  --    exp(kr + m * log(2)) = 2^m * exp(r)^k
  -- where m and k are integers.  By choosing m appropriately
  -- we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
  -- evaluated using the familiar Taylor series.  Reducing the
  -- argument substantially speeds up the convergence.

    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

    res : deca_double;
    k : constant double_float := C_ldexp(1.0,26);
    inv_k : constant double_float := 1.0/k;
    e0 : constant double_float :=  2.71828182845904509E+00;
    e1 : constant double_float :=  1.44564689172925016E-16;
    e2 : constant double_float := -2.12771710803817676E-33;
    e3 : constant double_float :=  1.51563015984121914E-49;
    e4 : constant double_float := -9.33538137882084738E-66;
    e5 : constant double_float := -2.02185260335762112E-82;
    e6 : constant double_float :=  8.68351053192660636E-99;
    e7 : constant double_float := -7.14388345698345838E-115;
    e8 : constant double_float := -4.29981073684482334E-131;
    e9 : constant double_float := -5.09120833097630782E-149;
    exp1 : constant deca_double := Create(e0,e1,e2,e3,e4,e5,e6,e7,e8,e9);
    L0 : constant double_float :=  6.93147180559945286E-01;
    L1 : constant double_float :=  2.31904681384629956E-17;
    L2 : constant double_float :=  5.70770843841621207E-34;
    L3 : constant double_float := -3.58243221060181142E-50;
    L4 : constant double_float := -1.35216967579886296E-66;
    L5 : constant double_float :=  6.08063874024081391E-83;
    L6 : constant double_float :=  2.89550243323471469E-99;
    L7 : constant double_float :=  2.35138671214564106E-116;
    L8 : constant double_float :=  4.45977441701428101E-133;
    L9 : constant double_float := -3.06993326323528609E-149;
    log2 : constant deca_double := Create(L0,L1,L2,L3,L4,L5,L6,L7,L8,L9);
    da_eps : constant double_float := 5.6902623986817984e-160; -- 2^-529
    tol : constant double_float := inv_k*da_eps;
    m : constant double_float := double_float'floor(x.right_thumb/L0 + 0.5);
    i_fac : array(0..14) of deca_double;
      -- inverse factorials for Taylor expansion
    p,s,t : deca_double;
    cnt : integer;

  begin
    if x.right_thumb <= -709.0 then
      res := Create(0.0);
    elsif x.right_thumb >= 709.0 then
      res := Create(-1.0);
    elsif is_zero(x) then
      res := Create(1.0);
    elsif is_one(x) then
      res := exp1;
    else
      i_fac(0) := Create(1.66666666666666657E-01, 9.25185853854297066E-18,
                         5.13581318503262866E-34, 2.85094902409834186E-50,
                         1.58259462429329970E-66, 8.78516495269210220E-83,
                         4.87674620280437299E-99, 2.70713795980335902E-115,
                         1.50276344690523035E-131, 8.34201289659658680E-148);
      i_fac(1) := Create(4.16666666666666644E-02, 2.31296463463574266E-18,
                         1.28395329625815716E-34, 7.12737256024585466E-51,
                         3.95648656073324926E-67, 2.19629123817302555E-83,
                         1.21918655070109325E-99, 6.76784489950839755E-116,
                         3.75690861726307587E-132, 2.08550322414909669E-148);
      i_fac(2) := Create(8.33333333333333322E-03, 1.15648231731787138E-19,
                         1.60494162032269652E-36, 2.22730392507682967E-53,
                         3.09100512557285111E-70, 4.28963132455669071E-87,
                         5.95305932959518212E-104, 8.26152941834521220E-121,
                         1.14651752235811646E-137, 1.59111268581469358E-154);
      i_fac(3) := Create(1.38888888888888894E-03, -5.30054395437357706E-20,
                        -1.73868675534958776E-36, -1.63335621172300840E-52,
                        -5.35774221765960817E-69, -5.03316742081318367E-85,
                        -1.65098178740773037E-101, -1.55096445613734114E-117,
                        -5.08748041921041505E-134, -4.77927822200816229E-150);
      i_fac(4) := Create(1.98412698412698413E-04, 1.72095582934207053E-22,
                         1.49269123913941271E-40, 1.29470326746002471E-58,
                         1.12297607624339191E-76, 9.74026481209866379E-95,
                         8.44833301588918793E-113, 7.32776080776645696E-131,
                         6.35581934978762111E-149, -1.43568009093965205E-165);
      i_fac(5) := Create(2.48015873015873016E-05, 2.15119478667758816E-23,
                         1.86586404892426588E-41, 1.61837908432503088E-59,
                         1.40372009530423989E-77, 1.21753310151233297E-95,
                         1.05604162698614849E-113, 9.15970100970807120E-132,
                         7.94477418723452422E-150, -4.59519619682793404E-166);
      i_fac(6) := Create(2.75573192239858925E-06, -1.85839327404647208E-22,
                         8.49175460488199287E-39, -5.72661640789429621E-55,
                         2.61672391582886888E-71, -1.76464992319726308E-87,
                         8.06340311310246955E-104, -5.43774740551399185E-120,
                         2.48472792147028858E-136, -1.67563528932629318E-152);
      i_fac(7) := Create(2.75573192239858883E-07, 2.37677146222502973E-23,
                        -3.26318890334088294E-40, 1.61435111860404415E-56,
                        -9.99798078191450073E-74, -5.23082523455780306E-91,
                        -2.46767621898134459E-107, 8.08497474286222354E-124,
                        -4.24517638331796814E-140, 6.26790008230286277E-157);
      i_fac(8) := Create(2.50521083854417202E-08, -1.44881407093591197E-24,
                         2.04267351467144546E-41, -8.49632672007163175E-58,
                        -9.08907343810409139E-75, -2.26065446324988261E-91,
                        -7.19805892199252488E-108, -6.40212493137993289E-125,
                        -1.95076372700295719E-141, 5.69809073654579017E-158);
      i_fac(9) := Create(2.08767569878681002E-09, -1.20734505911325997E-25,
                         1.70222792889287100E-42, 1.41609532150396700E-58,
                         5.13820161376913309E-75, 1.44797661649508516E-91,
                         1.30256332445267049E-107, 3.72847700075020451E-124,
                         8.58467087114439999E-141, 3.93203656395581129E-157);
      i_fac(10) := Create(1.60590438368216133E-10, 1.25852945887520981E-26,
                         -5.31334602762985031E-43, 3.54021472597605528E-59,
                         -1.98567896059148884E-75, -1.02148490610754577E-91,
                          5.19442455358700245E-108, 2.03226501937080225E-124,
                         -2.56938883085768487E-141, -8.18079628975045645E-158);
      i_fac(11) := Create(1.14707455977297245E-11, 2.06555127528307454E-28,
                          6.88907923246664603E-45, 5.72920002655109095E-61,
                         -3.65551458930952487E-78, -1.73752114063892946E-94,
                         -1.29464103159491163E-110, -4.67626589929406776E-127,
                         -7.80208024976360939E-144, 9.41483637653369274E-162);
      i_fac(12) := Create(7.64716373181981641E-13, 7.03872877733453001E-30,
                         -7.82753927716258345E-48, 1.92138649443790242E-64,
                         -1.43027453759388186E-80, 7.76151305461206293E-97,
                          3.19012607380048399E-114, 1.78550560936942146E-130,
                          4.02921430032281832E-148, 2.66951511873000631E-164);
      i_fac(13) := Create(4.77947733238738525E-14, 4.39920548583408126E-31,
                         -4.89221204822661465E-49, 1.20086655902368901E-65,
                         -8.93921585996176165E-82, 4.85094565913253933E-98,
                          1.99382879612530249E-115, 1.11594100585588841E-131,
                          2.51825893770176058E-149, 9.05596094344050963E-166);
      i_fac(14) := Create(2.81145725434552060E-15, 1.65088427308614326E-31,
                         -2.87777179307447918E-50, 4.27110689256293549E-67,
                         -2.93287743014724397E-83, -1.23436334109628973E-99,
                         -2.14851021924804564E-118, 2.19396004831602516E-134,
                          4.74996808318817135E-151, 2.77386830223949153E-168);
      res := mul_pwr2(x - m*log2,inv_k);
      p := res*res;
      s := res + mul_pwr2(p,0.5);
      cnt := 0;
      loop
        p := p*res;
        t := i_fac(cnt)*p;
        s := s + t;
        exit when abs(t.right_thumb) <= tol;
        cnt := cnt + 1;
        exit when (cnt >= 15);
      end loop;
      for i in 1..26 loop -- 26 times s = mul_pwr2(s,2.0) + sqr(s);
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

  function log ( x : deca_double ) return deca_double is

  -- Strategy.  The Taylor series for log converges much more
  --   slowly than that of exp, due to the lack of the factorial
  --   term in the denominator.  Hence this routine instead tries
  --   to determine the root of the function f(x) = exp(x) - a
  --   using Newton iteration.  The iteration is given by
  --       x' = x - f(x)/f'(x)
  --          = x - (1 - a * exp(-x))
  --          = x + a * exp(-x) - 1.
  --   Six iterations are needed, since Newton's iteration
  --   approximately doubles the number of digits per iteration.

    res,acc : deca_double;

  begin
    if is_one(x) then
      res := Create(0.0);
    elsif x.right_thumb <= 0.0 then
      put_line("Dd_log: argument is not positive");
      res := Create(-1.0);
    else
      res := Create(Standard_Mathematical_Functions.LN(x.right_thumb));
      for i in 1..6 loop       -- r = r + x*exp(-r) - 1
        acc := x*exp(-res);    -- acc = x*exp(-res)      
        res := res + acc;      -- res = res + x*exp(-res)
        res := res - 1.0;
      end loop;
    end if;
    return res;
  end log;

  function log10 ( x : deca_double ) return deca_double is

    res : deca_double;
    log10_0 : constant double_float :=  2.30258509299404590E+00;
    log10_1 : constant double_float := -2.17075622338224935E-16;
    log10_2 : constant double_float := -9.98426245446577657E-33;
    log10_3 : constant double_float := -4.02335745445020638E-49;
    log10_4 : constant double_float :=  1.92889952896933719E-65;
    log10_5 : constant double_float := -5.21257011815125513E-82;
    log10_6 : constant double_float := -2.60373698986932938E-98;
    log10_7 : constant double_float :=  8.29741762082190114E-115;
    log10_8 : constant double_float := -4.10600981629892264E-131;
    log10_9 : constant double_float := -5.26148805167447880E-148;
    logten : constant deca_double
           := create(log10_0,log10_1,log10_2,log10_3,log10_4,
                     log10_5,log10_6,log10_7,log10_8,log10_9);

  begin
    res := log(x)/logten;
    return res;
  end log10;

-- ARITHMETICAL OPERATIONS AS PROCEDURES :

  procedure Add ( x : in out deca_double; y : in deca_double ) is
  begin
    x := x + y;
  end Add;

  procedure Sub ( x : in out deca_double; y : in deca_double ) is
  begin
    x := x - y;
  end Sub;

  procedure Min ( x : in out deca_double ) is
  begin
    x := -x;
  end Min;

  procedure Mul ( x : in out deca_double; y : in deca_double ) is
  begin
    x := x*y;
  end Mul;

  procedure Div ( x : in out deca_double; y : in deca_double ) is
  begin
    x := x/y;
  end Div;

-- DESTRUCTOR :

  procedure clear ( x : in out deca_double ) is
  begin
    x.right_thumb := 0.0; x.right_index := 0.0;
    x.right_middle := 0.0; x.right_ring := 0.0; x.right_pink := 0.0;
    x.left_thumb := 0.0; x.left_index := 0.0;
    x.left_middle := 0.0; x.left_ring := 0.0; x.left_pink := 0.0;
  end clear;

end Deca_Double_Numbers;
