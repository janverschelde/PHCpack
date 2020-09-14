with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with Double_Double_Basics;
with Double_Double_Numbers;
with Fast_Double_Renormalizations;       use Fast_Double_Renormalizations;

package body Octo_Double_Numbers is

-- CONSTRUCTORS :

  function create ( i : integer ) return octo_double is

    res : constant octo_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( n : natural32 ) return octo_double is

    res : constant octo_double := create(double_float(n));

  begin
    return res;
  end create;

  function create ( n : natural64 ) return octo_double is

    res : constant octo_double := create(double_float(n));

  begin
    return res;
  end create;

  function create ( i : integer32 ) return octo_double is

    res : constant octo_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( i : integer64 ) return octo_double is

    res : constant octo_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( x : double_float ) return octo_double is

    res : octo_double;

  begin
    res.hihihi := x;
    res.lohihi := 0.0;
    res.hilohi := 0.0;
    res.lolohi := 0.0;
    res.hihilo := 0.0;
    res.lohilo := 0.0;
    res.hilolo := 0.0;
    res.lololo := 0.0;
    return res;
  end Create;

  function create ( hihihi,lohihi,hilohi,lolohi : double_float;
                    hihilo,lohilo,hilolo,lololo : double_float )
                  return octo_double is

    res : octo_double;

  begin
    res.hihihi := hihihi;
    res.lohihi := lohihi;
    res.hilohi := hilohi;
    res.lolohi := lolohi;
    res.hihilo := hihilo;
    res.lohilo := lohilo;
    res.hilolo := hilolo;
    res.lololo := lololo;
    return res;
  end create;

  function "abs" ( x : octo_double ) return octo_double is

    res : octo_double;

  begin
    if x.hihihi < 0.0 then
      res.hihihi := -x.hihihi; res.lohihi := -x.lohihi;
      res.hilohi := -x.hilohi; res.lolohi := -x.lolohi;
      res.hihilo := -x.hihilo; res.lohilo := -x.lohilo;
      res.hilolo := -x.hilolo; res.lololo := -x.lololo;
    else
      res.hihihi := x.hihihi; res.lohihi := x.lohihi;
      res.hilohi := x.hilohi; res.lolohi := x.lolohi;
      res.hihilo := x.hihilo; res.lohilo := x.lohilo;
      res.hilolo := x.hilolo; res.lololo := x.lololo;
    end if;
    return res;
  end "abs";

  function AbsVal ( x : octo_double ) return octo_double is

    res : octo_double;

  begin
    if x.hihihi < 0.0 then
      res.hihihi := -x.hihihi; res.lohihi := -x.lohihi;
      res.hilohi := -x.hilohi; res.lolohi := -x.lolohi;
      res.hihilo := -x.hihilo; res.lohilo := -x.lohilo;
      res.hilolo := -x.hilolo; res.lololo := -x.lololo;
    else
      res.hihihi := x.hihihi; res.lohihi := x.lohihi;
      res.hilohi := x.hilohi; res.lolohi := x.lolohi;
      res.hihilo := x.hihilo; res.lohilo := x.lohilo;
      res.hilolo := x.hilolo; res.lololo := x.lololo;
    end if;
    return res;
  end AbsVal;

  function floor ( x : octo_double ) return octo_double is

    res : octo_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8 : double_float;

  begin
    f0 := double_float'floor(x.hihihi);
    f1 := 0.0; f2 := 0.0; f3 := 0.0; f4 := 0.0;
    f5 := 0.0; f6 := 0.0; f7 := 0.0; f8 := 0.0;
    if f0 = x.hihihi then
      f1 := double_float'floor(x.lohihi);
      if f1 = x.lohihi then
        f2 := double_float'floor(x.hilohi);
        if f2 = x.hilohi then
          f3 := double_float'floor(x.lolohi);
          if f3 = x.lolohi then
            f4 := double_float'floor(x.hihilo);
            if f4 = x.hihilo then
              f5 := double_float'floor(x.lohilo);
              if f5 = x.lohilo then
                f6 := double_float'floor(x.hilolo);
                if f6 = x.hilolo then
                  f7 := double_float'floor(x.lololo);
                end if;
              end if;
            end if;
          end if;
        end if;
      end if;
    end if;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,
                res.hihihi,res.lohihi,res.hilohi,res.lolohi,
                res.hihilo,res.lohilo,res.hilolo,res.lololo);
    return res;
  end floor;

  function nint ( x : octo_double ) return octo_double is

    res : octo_double;
    x0,x1,x2,x3,x4,x5,x6,x7,x8 : double_float;

  begin
    x1 := 0.0; x2 := 0.0; x3 := 0.0; x4 := 0.0; x5 := 0.0;
    x6 := 0.0; x7 := 0.0; x8 := 0.0;
    x0 := Double_Double_Numbers.nint(x.hihihi);
    if x0 = x.hihihi then      -- first double is already an integer
      x1 := Double_Double_Numbers.nint(x.lohihi);
      if x1 = x.lohihi then    -- second double is already an integer
        x2 := Double_Double_Numbers.nint(x.hilohi);
        if x2 = x.hilohi then  -- third double is already an integer
          x3 := Double_Double_Numbers.nint(x.lolohi);
          if x3 = x.lolohi then  -- 4-th double is an integer
            x4 := Double_Double_Numbers.nint(x.hihilo);
            if x4 = x.hihilo then  -- 5-th double is an integer
              x5 := Double_Double_Numbers.nint(x.lohilo);
              if x5 = x.lohilo then  -- 6-th double is an integer
                x6 := Double_Double_Numbers.nint(x.hilolo);
                if x6 = x.hilolo then  -- 7-th double is an integer
                  x7 := Double_Double_Numbers.nint(x.lololo);
                else
                  if abs(x6 - x.hilolo) = 0.5 and x.lololo < 0.0
                   then x6 := x6 - 1.0;
                  end if;
                end if;
              else
                if abs(x5 - x.lohilo) = 0.5 and x.hilolo < 0.0
                 then x5 := x5 - 1.0;
                end if;
              end if;
            else
              if abs(x4 - x.hihilo) = 0.5 and x.lohilo < 0.0
               then x4 := x4 - 1.0;
              end if;
            end if;
          else
            if abs(x3 - x.lolohi) = 0.5 and x.hihilo < 0.0
             then x3 := x3 - 1.0;
            end if;
          end if;
        else
          if abs(x2 - x.hilohi) = 0.5 and x.lolohi < 0.0
           then x2 := x2 - 1.0;
          end if;
        end if;
      else
        if abs(x1 - x.lohihi) = 0.5 and x.hilohi < 0.0
         then x1 := x1 - 1.0;
        end if;
      end if;
    else                     -- first double is not an integer
      if abs(x0 - x.hihihi) = 0.5 and x.lohihi < 0.0
       then x0 := x0 - 1.0;
      end if;
    end if;
    fast_renorm(x0,x1,x2,x3,x4,x5,x6,x7,x8,
                res.hihihi,res.lohihi,res.hilohi,res.lolohi,
                res.hihilo,res.lohilo,res.hilolo,res.lololo);
    return res;
  end nint;

-- SELECTORS :

  function hihihi_part ( x : octo_double ) return double_float is
  begin
    return x.hihihi;
  end hihihi_part;

  function lohihi_part ( x : octo_double ) return double_float is
  begin
    return x.lohihi;
  end lohihi_part;

  function hilohi_part ( x : octo_double ) return double_float is
  begin
    return x.hilohi;
  end hilohi_part;

  function lolohi_part ( x : octo_double ) return double_float is
  begin
    return x.lolohi;
  end lolohi_part;

  function hihilo_part ( x : octo_double ) return double_float is
  begin
    return x.hihilo;
  end hihilo_part;

  function lohilo_part ( x : octo_double ) return double_float is
  begin
    return x.lohilo;
  end lohilo_part;

  function hilolo_part ( x : octo_double ) return double_float is
  begin
    return x.hilolo;
  end hilolo_part;

  function lololo_part ( x : octo_double ) return double_float is
  begin
    return x.lololo;
  end lololo_part;

-- TYPE CASTS :

  function to_int ( x : octo_double ) return integer32 is
  begin
    return integer32(x.hihihi);
  end to_int;

  function to_double ( x : octo_double ) return double_float is
  begin
    return x.hihihi;
  end to_double;

-- COMPARISON and COPYING :

  function is_zero ( x : octo_double ) return boolean is
  begin
    return ((x.hihihi = 0.0) and (x.lohihi = 0.0) and
            (x.hilohi = 0.0) and (x.lolohi = 0.0) and
            (x.hihilo = 0.0) and (x.lohilo = 0.0) and
            (x.hilolo = 0.0) and (x.lololo = 0.0));
  end is_zero;

  function is_one ( x : octo_double ) return boolean is
  begin
    return ((x.hihihi = 1.0) and (x.lohihi = 0.0) and
            (x.hilohi = 0.0) and (x.lolohi = 0.0) and
            (x.hihilo = 0.0) and (x.lohilo = 0.0) and
            (x.hilolo = 0.0) and (x.lololo = 0.0));
  end is_one;

  function is_positive ( x : octo_double ) return boolean is
  begin
    return (x.hihihi > 0.0);
  end is_positive;

  function is_negative ( x : octo_double ) return boolean is
  begin
    return (x.hihihi < 0.0);
  end is_negative;

  function equal ( x,y : octo_double ) return boolean is
  begin
    return ((x.hihihi = y.hihihi) and (x.lohihi = y.lohihi) and
            (x.hilohi = y.hilohi) and (x.lolohi = y.lolohi) and
            (x.hihilo = y.hihilo) and (x.lohilo = y.lohilo) and
            (x.hilolo = y.hilolo) and (x.lololo = y.lololo));
  end equal;

  function equal ( x : octo_double; y : double_float ) return boolean is
  begin
    return ((x.hihihi = y) and (x.lohihi = 0.0) and
            (x.hilohi = 0.0) and (x.lolohi = 0.0) and
            (x.hihilo = 0.0) and (x.lohilo = 0.0) and
            (x.hilolo = 0.0) and (x.lololo = 0.0));
  end equal;

  function "<" ( x,y : octo_double ) return boolean is
  begin
    return ((x.hihihi < y.hihihi)
         or (x.hihihi = y.hihihi and x.lohihi < y.lohihi)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi < y.hilohi)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi = y.hilohi and x.lolohi < y.lolohi)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi = y.hilohi and x.lolohi = y.lolohi and
             x.hihilo < y.hihilo)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi = y.hilohi and x.lolohi = y.lolohi and
             x.hihilo = y.hihilo and x.lohilo < y.lohilo)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi = y.hilohi and x.lolohi = y.lolohi and
             x.hihilo = y.hihilo and x.lohilo = y.lohilo and
             x.hilolo < y.hilolo)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi = y.hilohi and x.lolohi = y.lolohi and
             x.hihilo = y.hihilo and x.lohilo = y.lohilo and
             x.hilolo = y.hilolo and x.lololo < y.lololo));
  end "<";

  function "<" ( x : octo_double; y : double_float ) return boolean is
  begin
    return ((x.hihihi < y)
         or (x.hihihi = y and x.lohihi < 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi < 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi = 0.0 and
             x.lolohi < 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi = 0.0 and
             x.lolohi = 0.0 and x.hihilo < 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi = 0.0 and
             x.lolohi = 0.0 and x.hihilo = 0.0 and x.lohilo < 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi = 0.0 and
             x.lolohi = 0.0 and x.hihilo = 0.0 and x.lohilo = 0.0 and
             x.hilolo < 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi = 0.0 and
             x.lolohi = 0.0 and x.hihilo = 0.0 and x.lohilo = 0.0 and
             x.hilolo = 0.0 and x.lololo < 0.0));
  end "<";

  function "<" ( x : double_float; y : octo_double ) return boolean is
  begin
    return (y > x);
  end "<";

  function "<=" ( x,y : octo_double ) return boolean is
  begin
    return (x < y) or equal(x,y);
  end "<=";

  function "<=" ( x : octo_double; y : double_float ) return boolean is
  begin
    return (x < y) or equal(x,y);
  end "<=";

  function "<=" ( x : double_float; y : octo_double ) return boolean is
  begin
    return (y >= x);
  end "<=";

  function ">" ( x,y : octo_double ) return boolean is
  begin
    return ((x.hihihi > y.hihihi)
         or (x.hihihi = y.hihihi and x.lohihi > y.lohihi)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi > y.hilohi)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi = y.hilohi and x.lolohi > y.lolohi)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi = y.hilohi and x.lolohi = y.lolohi and
             x.hihilo > y.hihilo)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi = y.hilohi and x.lolohi = y.lolohi and
             x.hihilo = y.hihilo and x.lohilo > y.lohilo)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi = y.hilohi and x.lolohi = y.lolohi and
             x.hihilo = y.hihilo and x.lohilo = y.lohilo and
             x.hilolo > y.hilolo)
         or (x.hihihi = y.hihihi and x.lohihi = y.lohihi and
             x.hilohi = y.hilohi and x.lolohi = y.lolohi and
             x.hihilo = y.hihilo and x.lohilo = y.lohilo and
             x.hilolo = y.hilolo and x.lololo > y.lololo));
  end ">";

  function ">" ( x : octo_double; y : double_float ) return boolean is
  begin
    return ((x.hihihi > y)
         or (x.hihihi = y and x.lohihi > 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi > 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi = 0.0 and
             x.lolohi > 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi = 0.0 and
             x.lolohi = 0.0 and x.hihilo > 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi = 0.0 and
             x.lolohi = 0.0 and x.hihilo = 0.0 and x.lohilo > 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi = 0.0 and
             x.lolohi = 0.0 and x.hihilo = 0.0 and x.lohilo = 0.0 and
             x.hilolo > 0.0)
         or (x.hihihi = y and x.lohihi = 0.0 and x.hilohi = 0.0 and
             x.lolohi = 0.0 and x.hihilo = 0.0 and x.lohilo = 0.0 and
             x.hilolo = 0.0 and x.lololo > 0.0));
  end ">";

  function ">" ( x : double_float; y : octo_double ) return boolean is
  begin
    return (y < x);
  end ">";

  function ">=" ( x,y : octo_double ) return boolean is
  begin
    return (x > y) or equal(x,y);
  end ">=";

  function ">=" ( x : octo_double; y : double_float ) return boolean is
  begin
    return (x > y) or equal(x,y);
  end ">=";

  function ">=" ( x : double_float; y : octo_double ) return boolean is
  begin
    return (y <= x);
  end ">=";

  procedure copy ( x : in octo_double; y : in out octo_double ) is
  begin
    y.hihihi := x.hihihi; y.lohihi := x.lohihi;
    y.hilohi := x.hilohi; y.lolohi := x.lolohi;
    y.hihilo := x.hihilo; y.lohilo := x.lohilo;
    y.hilolo := x.hilolo; y.lololo := x.lololo;
  end copy;

-- ARITHMETICAL OPERATIONS :

  function "+" ( x,y : octo_double ) return octo_double is

  -- ALGORITHM : baileyAdd_fast<8,8,8> generated by CAMPARY.

    res : octo_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8,e : double_float;

  begin
    f8 := 0.0;
    Double_Double_Basics.two_sum(x.lololo,y.lololo,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(x.hilolo,y.hilolo,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(x.lohilo,y.lohilo,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(x.hihilo,y.hihilo,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(x.lolohi,y.lolohi,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(x.hilohi,y.hilohi,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(x.lohihi,y.lohihi,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(x.hihihi,y.hihihi,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,
                res.hihihi,res.lohihi,res.hilohi,res.lolohi,
                res.hihilo,res.lohilo,res.hilolo,res.lololo);
    return res;
  end "+";

  function "+" ( x : octo_double; y : double_float ) return octo_double is

    res : octo_double;

  begin
    renorm_add1(x.hihihi,x.lohihi,x.hilohi,x.lolohi,
                x.hihilo,x.lohilo,x.hilolo,x.lololo,y,
                res.hihihi,res.lohihi,res.hilohi,res.lolohi,
                res.hihilo,res.lohilo,res.hilolo,res.lololo);
    return res;
  end "+";

  function "+" ( x : double_float; y : octo_double ) return octo_double is

    res : constant octo_double := y + x;

  begin
    return res;
  end "+";

  function "+" ( x : octo_double ) return octo_double is

    res : octo_double;

  begin
    res.hihihi := x.hihihi;
    res.lohihi := x.lohihi;
    res.hilohi := x.hilohi;
    res.lolohi := x.lolohi;
    res.hihilo := x.hihilo;
    res.lohilo := x.lohilo;
    res.hilolo := x.hilolo;
    res.lololo := x.lololo;
    return res;
  end "+";

  function "-" ( x : octo_double ) return octo_double is

    res : octo_double;

  begin
    res.hihihi := -x.hihihi;
    res.lohihi := -x.lohihi;
    res.hilohi := -x.hilohi;
    res.lolohi := -x.lolohi;
    res.hihilo := -x.hihilo;
    res.lohilo := -x.lohilo;
    res.hilolo := -x.hilolo;
    res.lololo := -x.lololo;
    return res;
  end "-";

  function "-" ( x,y : octo_double ) return octo_double is

    mny : constant octo_double := -y;
    res : constant octo_double := x + mny;

  begin
    return res;
  end "-";

  function "-" ( x : octo_double; y : double_float ) return octo_double is

    mny : constant double_float := -y;
    res : constant octo_double := x + mny;

  begin
    return res;
  end "-";

  function "-" ( x : double_float; y : octo_double ) return octo_double is

    mny : constant octo_double := -y;
    res : constant octo_double := x + mny;

  begin
    return res;
  end "-";

  function "*" ( x,y : octo_double ) return octo_double is

  -- ALGORITHM : baileyMul_fast<8,8,8> generated by CAMPARY.

    res : octo_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8,p,e : double_float;

  begin
    f8 := x.lohihi*y.lololo;
    f8 := f8 + x.hilohi*y.hilolo;
    f8 := f8 + x.lolohi*y.lohilo;
    f8 := f8 + x.hihilo*y.hihilo;
    f8 := f8 + x.lohilo*y.lolohi;
    f8 := f8 + x.hilolo*y.hilohi;
    f8 := f8 + x.lololo*y.lohihi;
    Double_Double_Basics.two_prod(x.hihihi,y.lololo,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohihi,y.hilolo,p,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hilohi,y.lohilo,p,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lolohi,y.hihilo,p,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihilo,y.lolohi,p,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohilo,y.hilohi,p,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hilolo, y.lohihi,p,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lololo,y.hihihi,p,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihihi,y.hilolo,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohihi,y.lohilo,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hilohi,y.hihilo,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lolohi,y.lolohi,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihilo,y.hilohi,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohilo,y.lohihi,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hilolo,y.hihihi,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihihi,y.lohilo,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohihi,y.hihilo,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hilohi,y.lolohi,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lolohi,y.hilohi,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihilo,y.lohihi,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohilo,y.hihihi,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihihi,y.hihilo,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohihi,y.lolohi,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hilohi,y.hilohi,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lolohi,y.lohihi,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihilo,y.hihihi,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihihi,y.lolohi,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohihi,y.hilohi,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hilohi,y.lohihi,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lolohi,y.hihihi,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihihi,y.hilohi,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohihi,y.lohihi,p,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f2,p,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hilohi,y.hihihi,p,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f2,p,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihihi,y.lohihi,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohihi,y.hihihi,p,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_sum(f1,p,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihihi,y.hihihi,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,
                res.hihihi,res.lohihi,res.hilohi,res.lolohi,
                res.hihilo,res.lohilo,res.hilolo,res.lololo);
    return res;
  end "*";

  function "*" ( x : octo_double; y : double_float ) return octo_double is

  -- ALGORITHM : baileyMul_fast<8,1,8> generated by CAMPARY.

    res : octo_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8,e : double_float;

  begin
    f8 := 0.0;
    Double_Double_Basics.two_prod(x.lololo,y,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hilolo,y,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohilo,y,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihilo,y,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lolohi,y,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hilohi,y,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.lohihi,y,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    Double_Double_Basics.two_prod(x.hihihi,y,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    f8 := f8 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,
                res.hihihi,res.lohihi,res.hilohi,res.lolohi,
                res.hihilo,res.lohilo,res.hilolo,res.lololo);
    return res;
  end "*";

  function "*" ( x : double_float; y : octo_double ) return octo_double is

    res : constant octo_double := y*x;

  begin
    return res;
  end "*";

  function Mul_pwr2 ( x : octo_double; y : double_float )
                    return octo_double is

    res : octo_double;

  begin
    res.hihihi := x.hihihi*y; res.lohihi := x.lohihi*y;
    res.hilohi := x.hilohi*y; res.lolohi := x.lolohi*y;
    res.hihilo := x.hihilo*y; res.lohilo := x.lohilo*y;
    res.hilolo := x.hilolo*y; res.lololo := x.lololo*y;
    return res;
  end Mul_pwr2;

  procedure Mul_pwr2 ( x : in out octo_double; y : in double_float ) is
  begin
    x.hihihi := x.hihihi*y; x.lohihi := x.lohihi*y;
    x.hilohi := x.hilohi*y; x.lolohi := x.lolohi*y;
    x.hihilo := x.hihilo*y; x.lohilo := x.lohilo*y;
    x.hilolo := x.hilolo*y; x.lololo := x.lololo*y;
  end Mul_pwr2;

  function "/" ( x,y : octo_double ) return octo_double is

    res,acc : octo_double;
    q0,q1,q2,q3,q4,q5,q6,q7,q8 : double_float;

  begin
    q0 := x.hihihi/y.hihihi;   acc := q0*y; res := x - acc;
    q1 := res.hihihi/y.hihihi; acc := q1*y; res := res - acc;
    q2 := res.hihihi/y.hihihi; acc := q2*y; res := res - acc;
    q3 := res.hihihi/y.hihihi; acc := q3*y; res := res - acc;
    q4 := res.hihihi/y.hihihi; acc := q4*y; res := res - acc;
    q5 := res.hihihi/y.hihihi; acc := q5*y; res := res - acc;
    q6 := res.hihihi/y.hihihi; acc := q6*y; res := res - acc;
    q7 := res.hihihi/y.hihihi; acc := q7*y; res := res - acc;
    q8 := res.hihihi/y.hihihi;
    fast_renorm(q0,q1,q2,q3,q4,q5,q6,q7,q8,
                res.hihihi,res.lohihi,res.hilohi,res.lolohi,
                res.hihilo,res.lohilo,res.hilolo,res.lololo);
    return res;
  end "/";

  function "/" ( x : octo_double; y : double_float ) return octo_double is

    ody : constant octo_double := create(y);
    res : constant octo_double := x/ody;

  begin
    return res;
  end "/";

  function "/" ( x : double_float; y : octo_double ) return octo_double is

    odx : constant octo_double := create(x);
    res : constant octo_double := odx/y;

  begin
    return res;
  end "/";

  function sqr ( x : octo_double ) return octo_double is
  begin
    return x*x;
  end sqr;

  function "**" ( x : octo_double; n : integer ) return octo_double is

    res,acc : octo_double;
    absn : natural;

  begin
    if n = 0 then
      res.hihihi := 1.0; res.lohihi := 0.0;
      res.hilohi := 0.0; res.lolohi := 0.0;
      res.hihilo := 0.0; res.lohilo := 0.0;
      res.hilolo := 0.0; res.lololo := 0.0;
    else
      if n > 0
       then absn := n;
       else absn := -n;
      end if;
      res.hihihi := x.hihihi; res.lohihi := x.lohihi;
      res.hilohi := x.hilohi; res.lolohi := x.lolohi;
      res.hihilo := x.hihilo; res.lohilo := x.lohilo;
      res.hilolo := x.hilolo; res.lololo := x.lololo;
      acc.hihihi := 1.0; acc.lohihi := 0.0;
      acc.hilohi := 0.0; acc.lolohi := 0.0;
      acc.hihilo := 0.0; acc.lohilo := 0.0;
      acc.hilolo := 0.0; acc.lololo := 0.0;
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
        acc.hihihi := res.hihihi; acc.lohihi := res.lohihi;
        acc.hilohi := res.hilohi; acc.lolohi := res.lolohi;
        acc.hihilo := res.hihilo; acc.lohilo := res.lohilo;
        acc.hilolo := res.hilolo; acc.lololo := res.lololo;
      end if;
      if n < 0 then
        res := 1.0/acc;          -- compute reciprocal
      else 
        res.hihihi := acc.hihihi; res.lohihi := acc.lohihi;
        res.hilohi := acc.hilohi; res.lolohi := acc.lolohi;
        res.hihilo := acc.hihilo; res.lohilo := acc.lohilo;
        res.hilolo := acc.hilolo; res.lololo := acc.lololo;
      end if;
    end if;
    return res;
  end "**";

  function ldexp ( x : octo_double; n : integer ) return octo_double is

    res : octo_double;
    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

  begin
    res.hihihi := C_ldexp(x.hihihi,n);
    res.lohihi := C_ldexp(x.lohihi,n);
    res.hilohi := C_ldexp(x.hilohi,n);
    res.lolohi := C_ldexp(x.lolohi,n);
    res.hihilo := C_ldexp(x.hihilo,n);
    res.lohilo := C_ldexp(x.lohilo,n);
    res.hilolo := C_ldexp(x.hilolo,n);
    res.lololo := C_ldexp(x.lololo,n);
    return res;
  end ldexp;

  function "**" ( x,y : octo_double ) return octo_double is
  begin
    return exp(y*log(x));
  end "**";

  function "**" ( x : octo_double; y : double_float ) return octo_double is

    od_y : constant octo_double := create(y);

  begin
    return x**od_y;
  end "**";

  function exp ( x : octo_double ) return octo_double is

  -- Strategy:  We first reduce the size of x by noting that
  --    exp(kr + m * log(2)) = 2^m * exp(r)^k
  -- where m and k are integers.  By choosing m appropriately
  -- we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
  -- evaluated using the familiar Taylor series.  Reducing the
  -- argument substantially speeds up the convergence.

    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

    res : octo_double;
    k : constant double_float := C_ldexp(1.0,20);
    inv_k : constant double_float := 1.0/k;
    e0 : constant double_float :=  2.71828182845904509E+00;
    e1 : constant double_float :=  1.44564689172925016E-16;
    e2 : constant double_float := -2.12771710803817676E-33;
    e3 : constant double_float :=  1.51563015984121914E-49;
    e4 : constant double_float := -9.33538137882084738E-66;
    e5 : constant double_float := -2.02185260335762112E-82;
    e6 : constant double_float :=  8.68351053192660636E-99;
    e7 : constant double_float := -7.14388345698342502E-115;
    exp1 : constant octo_double := Create(e0,e1,e2,e3,e4,e5,e6,e7);
    L0 : constant double_float :=  6.93147180559945286E-01;
    L1 : constant double_float :=  2.31904681384629956E-17;
    L2 : constant double_float :=  5.70770843841621207E-34;
    L3 : constant double_float := -3.58243221060181142E-50;
    L4 : constant double_float := -1.35216967579886296E-66;
    L5 : constant double_float :=  6.08063874024081391E-83;
    L6 : constant double_float :=  2.89550243323471469E-99;
    L7 : constant double_float :=  2.35138671214606540E-116;
    log2 : constant octo_double := Create(L0,L1,L2,L3,L4,L5,L6,L7);
    od_eps : constant double_float := 4.616489308892868e-128; -- 2^(-423)
    tol : constant double_float := inv_k*od_eps;
    m : constant double_float := double_float'floor(x.hihihi/L0 + 0.5);
    i_fac : array(0..14) of octo_double;
      -- inverse factorials for Taylor expansion
    p,s,t : octo_double;
    cnt : integer;

  begin
    if x.hihihi <= -709.0 then
      res := Create(0.0);
    elsif x.hihihi >= 709.0 then
      res := Create(-1.0);
    elsif is_zero(x) then
      res := Create(1.0);
    elsif is_one(x) then
      res := exp1;
    else
      i_fac(0) := Create(1.66666666666666657E-01, 9.25185853854297066E-18,
                         5.13581318503262866E-34, 2.85094902409834186E-50,
                         1.58259462429329970E-66, 8.78516495269210220E-83,
                         4.87674620280437299E-99, 2.70713795980339238E-115);
      i_fac(1) := Create(4.16666666666666644E-02, 2.31296463463574266E-18,
                         1.28395329625815716E-34, 7.12737256024585466E-51,
                         3.95648656073324926E-67, 2.19629123817302555E-83,
                         1.21918655070109325E-99,  6.76784489950843137E-116);
      i_fac(2) := Create(8.33333333333333322E-03, 1.15648231731787138E-19,
                         1.60494162032269652E-36, 2.22730392507682967E-53,
                         3.09100512557285111E-70, 4.28963132455669071E-87,
                         5.95305932959518212E-104, 8.26152941801411743E-121);
      i_fac(3) := Create(1.38888888888888894E-03, -5.30054395437357706E-20,
                        -1.73868675534958776E-36, -1.63335621172300840E-52,
                        -5.35774221765960817E-69, -5.03316742081318367E-85,
                        -1.65098178740773037E-101, -1.55096445613733005E-117);
      i_fac(4) := Create(1.98412698412698413E-04, 1.72095582934207053E-22,
                         1.49269123913941271E-40, 1.29470326746002471E-58,
                         1.12297607624339191E-76, 9.74026481209866379E-95,
                         8.44833301588918793E-113, 7.48646943949686349E-131);
      i_fac(5) := Create(2.48015873015873016E-05, 2.15119478667758816E-23,
                         1.86586404892426588E-41, 1.61837908432503088E-59,
                         1.40372009530423989E-77, 1.21753310151233297E-95,
                         1.05604162698614849E-113, 8.85811633556463122E-132);
      i_fac(6) := Create(2.75573192239858925E-06, -1.85839327404647208E-22,
                         8.49175460488199287E-39, -5.72661640789429621E-55,
                         2.61672391582886888E-71, -1.76464992319726308E-87,
                         8.06340311310246955E-104, -5.43774740551394713E-120);
      i_fac(7) := Create(2.75573192239858883E-07, 2.37677146222502973E-23,
                        -3.26318890334088294E-40, 1.61435111860404415E-56,
                        -9.99798078191450073E-74, -5.23082523455780306E-91,
                        -2.46767621898134459E-107, 8.08497474290649098E-124);
      i_fac(8) := Create(2.50521083854417202E-08, -1.44881407093591197E-24,
                         2.04267351467144546E-41, -8.49632672007163175E-58,
                        -9.08907343810409139E-75, -2.26065446324988261E-91,
                        -7.19805892199252488E-108, -6.40212493138514344E-125);
      i_fac(9) := Create(2.08767569878681002E-09, -1.20734505911325997E-25,
                         1.70222792889287100E-42, 1.41609532150396700E-58,
                         5.13820161376913309E-75, 1.44797661649508516E-91,
                         1.30256332445267049E-107, 3.72847700074999458E-124);
      i_fac(10) := Create(1.60590438368216133E-10, 1.25852945887520981E-26,
                         -5.31334602762985031E-43, 3.54021472597605528E-59,
                         -1.98567896059148884E-75, -1.02148490610754577E-91,
                          5.19442455358700245E-108, 2.03226501937076278E-124);
      i_fac(11) := Create(1.14707455977297245E-11, 2.06555127528307454E-28,
                          6.88907923246664603E-45, 5.72920002655109095E-61,
                         -3.65551458930952487E-78, -1.73752114063892946E-94,
                         -1.29464103159491163E-110, -4.67626589929544217E-127);
      i_fac(12) := Create(7.64716373181981641E-13, 7.03872877733453001E-30,
                         -7.82753927716258345E-48, 1.92138649443790242E-64,
                         -1.43027453759388186E-80, 7.76151305461206293E-97,
                          3.19012607380048399E-114, 1.78550560934448071E-130);
      i_fac(13) := Create(4.77947733238738525E-14, 4.39920548583408126E-31,
                         -4.89221204822661465E-49, 1.20086655902368901E-65,
                         -8.93921585996176165E-82, 4.85094565913253933E-98,
                          1.99382879612530249E-115, 1.11594100586527568E-131);
      i_fac(14) := Create(2.81145725434552060E-15, 1.65088427308614326E-31,
                         -2.87777179307447918E-50, 4.27110689256293549E-67,
                         -2.93287743014724397E-83, -1.23436334109628973E-99,
                         -2.14851021924804564E-118, 2.19396005239866827E-134);
      res := mul_pwr2(x - m*log2,inv_k);
      p := res*res;
      s := res + mul_pwr2(p,0.5);
      cnt := 0;
      loop
        p := p*res;
        t := i_fac(cnt)*p;
        s := s + t;
        exit when abs(t.hihihi) <= tol;
        cnt := cnt + 1;
        exit when (cnt >= 15);
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

  function log ( x : octo_double ) return octo_double is

  -- Strategy.  The Taylor series for log converges much more
  --   slowly than that of exp, due to the lack of the factorial
  --   term in the denominator.  Hence this routine instead tries
  --   to determine the root of the function f(x) = exp(x) - a
  --   using Newton iteration.  The iteration is given by
  --       x' = x - f(x)/f'(x)
  --          = x - (1 - a * exp(-x))
  --          = x + a * exp(-x) - 1.
  --   Five iterations are needed, since Newton's iteration
  --   approximately doubles the number of digits per iteration.

    res,acc : octo_double;

  begin
    if is_one(x) then
      res := Create(0.0);
    elsif x.hihihi <= 0.0 then
      put_line("td_log: argument is not positive");
      res := Create(-1.0);
    else
      res := Create(Standard_Mathematical_Functions.LN(x.hihihi));
      for i in 1..5 loop       -- r = r + x*exp(-r) - 1
        acc := x*exp(-res);    -- acc = x*exp(-res)      
        res := res + acc;      -- res = res + x*exp(-res)
        res := res - 1.0;
      end loop;
    end if;
    return res;
  end log;

  function log10 ( x : octo_double ) return octo_double is

    res : octo_double;
    log10_0 : constant double_float :=  2.30258509299404590E+00;
    log10_1 : constant double_float := -2.17075622338224935E-16;
    log10_2 : constant double_float := -9.98426245446577657E-33;
    log10_3 : constant double_float := -4.02335745445020638E-49;
    log10_4 : constant double_float :=  1.92889952896933719E-65;
    log10_5 : constant double_float := -5.21257011815125513E-82;
    log10_6 : constant double_float := -2.60373698986932938E-98;
    log10_7 : constant double_float :=  8.29741762082224918E-115;
    logten : constant octo_double
           := create(log10_0,log10_1,log10_2,log10_3,
                     log10_4,log10_5,log10_6,log10_7);

  begin
    res := log(x)/logten;
    return res;
  end log10;

-- ARITHMETICAL OPERATIONS AS PROCEDURES :

  procedure Add ( x : in out octo_double; y : in octo_double ) is
  begin
    x := x + y;
  end Add;

  procedure Sub ( x : in out octo_double; y : in octo_double ) is
  begin
    x := x - y;
  end Sub;

  procedure Min ( x : in out octo_double ) is
  begin
    x := -x;
  end Min;

  procedure Mul ( x : in out octo_double; y : in octo_double ) is
  begin
    x := x*y;
  end Mul;

  procedure Div ( x : in out octo_double; y : in octo_double ) is
  begin
    x := x/y;
  end Div;

-- DESTRUCTOR :

  procedure clear ( x : in out octo_double ) is
  begin
    x.hihihi := 0.0; x.lohihi := 0.0;
    x.hilohi := 0.0; x.lolohi := 0.0;
    x.hihilo := 0.0; x.lohilo := 0.0;
    x.hilolo := 0.0; x.lololo := 0.0;
  end clear;

end Octo_Double_Numbers;
