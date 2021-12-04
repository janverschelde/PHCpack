-- with text_io;                            use text_io;
-- with Standard_Mathematical_Functions;
with Double_Double_Basics;
with Double_Double_Numbers;
with Fast_Double_Renormalizations;       use Fast_Double_Renormalizations;

package body Hexa_Double_Numbers is

-- CONSTRUCTORS :

  function create ( i : integer ) return hexa_double is

    res : constant hexa_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( n : natural32 ) return hexa_double is

    res : constant hexa_double := create(double_float(n));

  begin
    return res;
  end create;

  function create ( n : natural64 ) return hexa_double is

    res : constant hexa_double := create(double_float(n));

  begin
    return res;
  end create;

  function create ( i : integer32 ) return hexa_double is

    res : constant hexa_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( i : integer64 ) return hexa_double is

    res : constant hexa_double := create(double_float(i));

  begin
    return res;
  end create;

  function create ( x : double_float ) return hexa_double is

    res : hexa_double;

  begin
    res.hihihihi := x;
    res.lohihihi := 0.0;
    res.hilohihi := 0.0;
    res.lolohihi := 0.0;
    res.hihilohi := 0.0;
    res.lohilohi := 0.0;
    res.hilolohi := 0.0;
    res.lololohi := 0.0;
    res.hihihilo := 0.0;
    res.lohihilo := 0.0;
    res.hilohilo := 0.0;
    res.lolohilo := 0.0;
    res.hihilolo := 0.0;
    res.lohilolo := 0.0;
    res.hilololo := 0.0;
    res.lolololo := 0.0;
    return res;
  end create;

  function create ( hihihihi,lohihihi,hilohihi,lolohihi : double_float;
                    hihilohi,lohilohi,hilolohi,lololohi : double_float;
                    hihihilo,lohihilo,hilohilo,lolohilo : double_float;
                    hihilolo,lohilolo,hilololo,lolololo : double_float )
                  return hexa_double is

    res : hexa_double;

  begin
    res.hihihihi := hihihihi;
    res.lohihihi := lohihihi;
    res.hilohihi := hilohihi;
    res.lolohihi := lolohihi;
    res.hihilohi := hihilohi;
    res.lohilohi := lohilohi;
    res.hilolohi := hilolohi;
    res.lololohi := lololohi;
    res.hihihilo := hihihilo;
    res.lohihilo := lohihilo;
    res.hilohilo := hilohilo;
    res.lolohilo := lolohilo;
    res.hihilolo := hihilolo;
    res.lohilolo := lohilolo;
    res.hilololo := hilololo;
    res.lolololo := lolololo;
    return res;
  end create;

  function "abs" ( x : hexa_double ) return hexa_double is

    res : hexa_double;

  begin
    if x.hihihihi < 0.0 then
      res.hihihihi := -x.hihihihi; res.lohihihi := -x.lohihihi;
      res.hilohihi := -x.hilohihi; res.lolohihi := -x.lolohihi;
      res.hihilohi := -x.hihilohi; res.lohilohi := -x.lohilohi;
      res.hilolohi := -x.hilolohi; res.lololohi := -x.lololohi;
      res.hihihilo := -x.hihihilo; res.lohihilo := -x.lohihilo;
      res.hilohilo := -x.hilohilo; res.lolohilo := -x.lolohilo;
      res.hihilolo := -x.hihilolo; res.lohilolo := -x.lohilolo;
      res.hilololo := -x.hilololo; res.lolololo := -x.lolololo;
    else
      res.hihihihi := x.hihihihi; res.lohihihi := x.lohihihi;
      res.hilohihi := x.hilohihi; res.lolohihi := x.lolohihi;
      res.hihilohi := x.hihilohi; res.lohilohi := x.lohilohi;
      res.hilolohi := x.hilolohi; res.lololohi := x.lololohi;
      res.hihihilo := x.hihihilo; res.lohihilo := x.lohihilo;
      res.hilohilo := x.hilohilo; res.lolohilo := x.lolohilo;
      res.hihilolo := x.hihilolo; res.lohilolo := x.lohilolo;
      res.hilololo := x.hilololo; res.lolololo := x.lolololo;
    end if;
    return res;
  end "abs";

  function AbsVal ( x : hexa_double ) return hexa_double is

    res : hexa_double;

  begin
    if x.hihihihi < 0.0 then
      res.hihihihi := -x.hihihihi; res.lohihihi := -x.lohihihi;
      res.hilohihi := -x.hilohihi; res.lolohihi := -x.lolohihi;
      res.hihilohi := -x.hihilohi; res.lohilohi := -x.lohilohi;
      res.hilolohi := -x.hilolohi; res.lololohi := -x.lololohi;
      res.hihihilo := -x.hihihilo; res.lohihilo := -x.lohihilo;
      res.hilohilo := -x.hilohilo; res.lolohilo := -x.lolohilo;
      res.hihilolo := -x.hihilolo; res.lohilolo := -x.lohilolo;
      res.hilololo := -x.hilololo; res.lolololo := -x.lolololo;
    else
      res.hihihihi := x.hihihihi; res.lohihihi := x.lohihihi;
      res.hilohihi := x.hilohihi; res.lolohihi := x.lolohihi;
      res.hihilohi := x.hihilohi; res.lohilohi := x.lohilohi;
      res.hilolohi := x.hilolohi; res.lololohi := x.lololohi;
      res.hihihilo := x.hihihilo; res.lohihilo := x.lohihilo;
      res.hilohilo := x.hilohilo; res.lolohilo := x.lolohilo;
      res.hihilolo := x.hihilolo; res.lohilolo := x.lohilolo;
      res.hilololo := x.hilololo; res.lolololo := x.lolololo;
    end if;
    return res;
  end AbsVal;

  function floor ( x : hexa_double ) return hexa_double is

    res : hexa_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16 : double_float;

  begin
    f0 := double_float'floor(x.hihihihi);
    f1 := 0.0; f2 := 0.0; f3 := 0.0; f4 := 0.0;
    f5 := 0.0; f6 := 0.0; f7 := 0.0; f8 := 0.0;
    f9 := 0.0; f10 := 0.0; f11 := 0.0; f12 := 0.0;
    f13 := 0.0; f14 := 0.0; f15 := 0.0; f16 := 0.0;
    if f0 = x.hihihihi then
      f1 := double_float'floor(x.lohihihi);
      if f1 = x.lohihihi then
        f2 := double_float'floor(x.hilohihi);
        if f2 = x.hilohihi then
          f3 := double_float'floor(x.lolohihi);
          if f3 = x.lolohihi then
            f4 := double_float'floor(x.hihilohi);
            if f4 = x.hihilohi then
              f5 := double_float'floor(x.lohilohi);
              if f5 = x.lohilohi then
                f6 := double_float'floor(x.hilolohi);
                if f6 = x.hilolohi then
                  f7 := double_float'floor(x.lololohi);
                  if f7 = x.lololohi then
                    f8 := double_float'floor(x.hihihilo);
                    if f8 = x.hihihilo then
                      f9 := double_float'floor(x.lohihilo);
                      if f9 = x.lohihilo then
                        f10 := double_float'floor(x.hilohilo);
                        if f10 = x.hilohilo then
                          f11 := double_float'floor(x.lolohilo);
                          if f11 = x.lolohilo then
                            f12 := double_float'floor(x.hihilolo);
                            if f12 = x.hihilolo then
                              f13 := double_float'floor(x.lohilolo);
                              if f13 = x.lohilolo then
                                f14 := double_float'floor(x.hilololo);
                                if f14 = x.hilololo then
                                  f15 := double_float'floor(x.lolololo);
                                end if;
                              end if;
                            end if;
                          end if;
                        end if;
                      end if;
                    end if;
                  end if;
                end if;
              end if;
            end if;
          end if;
        end if;
      end if;
    end if;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                res.hihihihi,res.lohihihi,res.hilohihi,res.lolohihi,
                res.hihilohi,res.lohilohi,res.hilolohi,res.lololohi,
                res.hihihilo,res.lohihilo,res.hilohilo,res.lolohilo,
                res.hihilolo,res.lohilolo,res.hilololo,res.lolololo);
    return res;
  end floor;

  function nint ( x : hexa_double ) return hexa_double is

    res : hexa_double;
    x0,x1,x2,x3,x4,x5,x6,x7,x8 : double_float;
    x9,x10,x11,x12,x13,x14,x15,x16 : double_float;

  begin
    x1 := 0.0; x2 := 0.0; x3 := 0.0; x4 := 0.0; x5 := 0.0;
    x6 := 0.0; x7 := 0.0; x8 := 0.0; x9 := 0.0; x10 := 0.0;
    x11 := 0.0; x12 := 0.0; x13 := 0.0; x14 := 0.0; x15 := 0.0;
    x16 := 0.0;
    x0 := Double_Double_Numbers.nint(x.hihihihi);
    if x0 = x.hihihihi then      -- first double is already an integer
      x1 := Double_Double_Numbers.nint(x.lohihihi);
      if x1 = x.lohihihi then    -- second double is already an integer
        x2 := Double_Double_Numbers.nint(x.hilohihi);
        if x2 = x.hilohihi then  -- third double is already an integer
          x3 := Double_Double_Numbers.nint(x.lolohihi);
          if x3 = x.lolohihi then  -- 4-th double is an integer
            x4 := Double_Double_Numbers.nint(x.hihilohi);
            if x4 = x.hihilohi then  -- 5-th double is an integer
              x5 := Double_Double_Numbers.nint(x.lohilohi);
              if x5 = x.lohilohi then  -- 6-th double is an integer
                x6 := Double_Double_Numbers.nint(x.hilolohi);
                if x6 = x.hilolohi then  -- 7-th double is an integer
                  x7 := Double_Double_Numbers.nint(x.lololohi);
                  if x7 = x.lololohi then  -- 8-th double is an integer
                    x8 := Double_Double_Numbers.nint(x.hihihilo);
                    if x8 = x.hihihilo then  -- 9-th double is an integer
                      x9 := Double_Double_Numbers.nint(x.lohihilo);
                      if x9 = x.lohihilo then  -- 10-th double is an integer
                        x10 := Double_Double_Numbers.nint(x.hilohilo);
                        if x10 = x.hilohilo then  -- 10-th double is integer
                          x11 := Double_Double_Numbers.nint(x.lolohilo);
                          if x11 = x.lolohilo then  -- 11-th double is integer
                            x12 := Double_Double_Numbers.nint(x.hihilolo);
                            if x12 = x.hihilolo then  -- 12-th double is int
                              x13 := Double_Double_Numbers.nint(x.lohilolo);
                              if x13 = x.lohilolo then  -- 13-th double is int
                                x14 := Double_Double_Numbers.nint(x.hilololo);
                                if x14 = x.hilololo then  -- 14-th is int
                                  x15 := Double_Double_Numbers.nint(x.lolololo);
                                else
                                  if abs(x14 - x.hilololo) = 0.5
                                    and x.lolololo < 0.0
                                   then x14 := x14 - 1.0;
                                  end if;
                                end if;
                              else
                                if abs(x13 - x.lohilolo) = 0.5
                                  and x.hilololo < 0.0
                                 then x13 := x13 - 1.0;
                                end if;
                              end if;
                            else
                              if abs(x12 - x.hihilolo) = 0.5
                                and x.lohilolo < 0.0
                               then x12 := x12 - 1.0;
                              end if;
                            end if;
                          else
                            if abs(x11 - x.lolohilo) = 0.5 and x.hihilolo < 0.0
                             then x11 := x11 - 1.0;
                            end if;
                          end if;
                        else
                          if abs(x10 - x.hilohilo) = 0.5 and x.lolohilo < 0.0
                           then x10 := x10 - 1.0;
                          end if;
                        end if;
                      else
                        if abs(x9 - x.lohihilo) = 0.5 and x.hilohilo < 0.0
                         then x9 := x9 - 1.0;
                        end if;
                      end if;
                    else
                      if abs(x8 - x.hihihilo) = 0.5 and x.lohihilo < 0.0
                       then x8 := x8 - 1.0;
                      end if;
                    end if;
                  else
                    if abs(x7 - x.lololohi) = 0.5 and x.hihihilo < 0.0
                     then x7 := x7 - 1.0;
                    end if;
                  end if;
                else
                  if abs(x6 - x.hilolohi) = 0.5 and x.lololohi < 0.0
                   then x6 := x6 - 1.0;
                  end if;
                end if;
              else
                if abs(x5 - x.lohilohi) = 0.5 and x.hilolohi < 0.0
                 then x5 := x5 - 1.0;
                end if;
              end if;
            else
              if abs(x4 - x.hihilohi) = 0.5 and x.lohilohi < 0.0
               then x4 := x4 - 1.0;
              end if;
            end if;
          else
            if abs(x3 - x.lolohihi) = 0.5 and x.hihilohi < 0.0
             then x3 := x3 - 1.0;
            end if;
          end if;
        else
          if abs(x2 - x.hilohihi) = 0.5 and x.lolohihi < 0.0
           then x2 := x2 - 1.0;
          end if;
        end if;
      else
        if abs(x1 - x.lohihihi) = 0.5 and x.hilohihi < 0.0
         then x1 := x1 - 1.0;
        end if;
      end if;
    else                     -- first double is not an integer
      if abs(x0 - x.hihihihi) = 0.5 and x.lohihihi < 0.0
       then x0 := x0 - 1.0;
      end if;
    end if;
    fast_renorm(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,
                res.hihihihi,res.lohihihi,res.hilohihi,res.lolohihi,
                res.hihilohi,res.lohilohi,res.hilolohi,res.lololohi,
                res.hihihilo,res.lohihilo,res.hilohilo,res.lolohilo,
                res.hihilolo,res.lohilolo,res.hilololo,res.lolololo);
    return res;
  end nint;

-- SELECTORS :

  function hihihihi_part ( x : hexa_double ) return double_float is
  begin
    return x.hihihihi;
  end hihihihi_part;

  function lohihihi_part ( x : hexa_double ) return double_float is
  begin
    return x.lohihihi;
  end lohihihi_part;

  function hilohihi_part ( x : hexa_double ) return double_float is
  begin
    return x.hilohihi;
  end hilohihi_part;

  function lolohihi_part ( x : hexa_double ) return double_float is
  begin
    return x.lolohihi;
  end lolohihi_part;

  function hihilohi_part ( x : hexa_double ) return double_float is
  begin
    return x.hihilohi;
  end hihilohi_part;

  function lohilohi_part ( x : hexa_double ) return double_float is
  begin
    return x.lohilohi;
  end lohilohi_part;

  function hilolohi_part ( x : hexa_double ) return double_float is
  begin
    return x.hilolohi;
  end hilolohi_part;

  function lololohi_part ( x : hexa_double ) return double_float is
  begin
    return x.lololohi;
  end lololohi_part;

  function hihihilo_part ( x : hexa_double ) return double_float is
  begin
    return x.hihihilo;
  end hihihilo_part;

  function lohihilo_part ( x : hexa_double ) return double_float is
  begin
    return x.lohihilo;
  end lohihilo_part;

  function hilohilo_part ( x : hexa_double ) return double_float is
  begin
    return x.hilohilo;
  end hilohilo_part;

  function lolohilo_part ( x : hexa_double ) return double_float is
  begin
    return x.lolohilo;
  end lolohilo_part;

  function hihilolo_part ( x : hexa_double ) return double_float is
  begin
    return x.hihilolo;
  end hihilolo_part;

  function lohilolo_part ( x : hexa_double ) return double_float is
  begin
    return x.lohilolo;
  end lohilolo_part;

  function hilololo_part ( x : hexa_double ) return double_float is
  begin
    return x.hilololo;
  end hilololo_part;

  function lolololo_part ( x : hexa_double ) return double_float is
  begin
    return x.lolololo;
  end lolololo_part;

-- TYPE CASTS :

  function to_int ( x : hexa_double ) return integer32 is
  begin
    return integer32(x.hihihihi);
  end to_int;

  function to_double ( x : hexa_double ) return double_float is
  begin
    return x.hihihihi;
  end to_double;

-- COMPARISON and COPYING :

  function is_zero ( x : hexa_double ) return boolean is
  begin
    return ((x.hihihihi = 0.0) and (x.lohihihi = 0.0) and
            (x.hilohihi = 0.0) and (x.lolohihi = 0.0) and
            (x.hihilohi = 0.0) and (x.lohilohi = 0.0) and
            (x.hilolohi = 0.0) and (x.lololohi = 0.0) and
            (x.hihihilo = 0.0) and (x.lohihilo = 0.0) and
            (x.hilohilo = 0.0) and (x.lolohilo = 0.0) and
            (x.hihilolo = 0.0) and (x.lohilolo = 0.0) and
            (x.hilololo = 0.0) and (x.lolololo = 0.0));
  end is_zero;

  function is_one ( x : hexa_double ) return boolean is
  begin
    return ((x.hihihihi = 1.0) and (x.lohihihi = 0.0) and
            (x.hilohihi = 0.0) and (x.lolohihi = 0.0) and
            (x.hihilohi = 0.0) and (x.lohilohi = 0.0) and
            (x.hilolohi = 0.0) and (x.lololohi = 0.0) and
            (x.hihihilo = 0.0) and (x.lohihilo = 0.0) and
            (x.hilohilo = 0.0) and (x.lolohilo = 0.0) and
            (x.hihilolo = 0.0) and (x.lohilolo = 0.0) and
            (x.hilololo = 0.0) and (x.lolololo = 0.0));
  end is_one;

  function is_positive ( x : hexa_double ) return boolean is
  begin
    return (x.hihihihi > 0.0);
  end is_positive;

  function is_negative ( x : hexa_double ) return boolean is
  begin
    return (x.hihihihi < 0.0);
  end is_negative;

  function equal ( x,y : hexa_double ) return boolean is
  begin
    return ((x.hihihihi = y.hihihihi) and (x.lohihihi = y.lohihihi) and
            (x.hilohihi = y.hilohihi) and (x.lolohihi = y.lolohihi) and
            (x.hihilohi = y.hihilohi) and (x.lohilohi = y.lohilohi) and
            (x.hilolohi = y.hilolohi) and (x.lololohi = y.lololohi) and
            (x.hihihilo = y.hihihilo) and (x.lohihilo = y.lohihilo) and
            (x.hilohilo = y.hilohilo) and (x.lolohilo = y.lolohilo) and
            (x.hihilolo = y.hihilolo) and (x.lohilolo = y.lohilolo) and
            (x.hilololo = y.hilololo) and (x.lolololo = y.lolololo));
  end equal;

  function equal ( x : hexa_double; y : double_float ) return boolean is
  begin
    return ((x.hihihihi = y) and (x.lohihihi = 0.0) and
            (x.hilohihi = 0.0) and (x.lolohihi = 0.0) and
            (x.hihilohi = 0.0) and (x.lohilohi = 0.0) and
            (x.hilolohi = 0.0) and (x.lololohi = 0.0) and
            (x.hihihilo = 0.0) and (x.lohihilo = 0.0) and
            (x.hilohilo = 0.0) and (x.lolohilo = 0.0) and
            (x.hihilolo = 0.0) and (x.lohilolo = 0.0) and
            (x.hilololo = 0.0) and (x.lolololo = 0.0));
  end equal;

  function "<" ( x,y : hexa_double ) return boolean is
  begin
    return ((x.hihihihi < y.hihihihi)
         or (x.hihihihi = y.hihihihi and x.lohihihi < y.lohihihi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi < y.hilohihi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi < y.lolohihi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi < y.hihilohi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi < y.lohilohi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi < y.hilolohi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi < y.lololohi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo < y.hihihilo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo < y.lohihilo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo < y.hilohilo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo = y.hilohilo and x.lolohilo < y.lolohilo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo = y.hilohilo and x.lolohilo = y.lolohilo and
             x.hihilolo < y.hihilolo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo = y.hilohilo and x.lolohilo = y.lolohilo and
             x.hihilolo = y.hihilolo and x.lohilolo < y.lohilolo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo = y.hilohilo and x.lolohilo = y.lolohilo and
             x.hihilolo = y.hihilolo and x.lohilolo = y.lohilolo and
             x.hilololo < y.hilololo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo = y.hilohilo and x.lolohilo = y.lolohilo and
             x.hihilolo = y.hihilolo and x.lohilolo = y.lohilolo and
             x.hilololo = y.hilololo and x.lolololo < y.lolololo));
  end "<";

  function "<" ( x : hexa_double; y : double_float ) return boolean is
  begin
    return ((x.hihihihi < y)
         or (x.hihihihi = y and x.lohihihi < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo = 0.0 and x.lolohilo < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo = 0.0 and x.lolohilo = 0.0 and
             x.hihilolo < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo = 0.0 and x.lolohilo = 0.0 and
             x.hihilolo = 0.0 and x.lohilolo < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo = 0.0 and x.lolohilo = 0.0 and
             x.hihilolo = 0.0 and x.lohilolo = 0.0 and
             x.hilololo < 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo = 0.0 and x.lolohilo = 0.0 and
             x.hihilolo = 0.0 and x.lohilolo = 0.0 and
             x.hilololo = 0.0 and x.lolololo < 0.0));
  end "<";

  function "<" ( x : double_float; y : hexa_double ) return boolean is
  begin
    return (y < x);
  end "<";

  function "<=" ( x,y : hexa_double ) return boolean is
  begin
    return (x < y) or equal(x,y);
  end "<=";

  function "<=" ( x : hexa_double; y : double_float ) return boolean is
  begin
    return (x < y) or equal(x,y);
  end "<=";

  function "<=" ( x : double_float; y : hexa_double ) return boolean is
  begin
    return (y >= x);
  end "<=";

  function ">" ( x,y : hexa_double ) return boolean is
  begin
    return ((x.hihihihi > y.hihihihi)
         or (x.hihihihi = y.hihihihi and x.lohihihi > y.lohihihi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi > y.hilohihi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi > y.lolohihi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi > y.hihilohi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi > y.lohilohi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi > y.hilolohi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi > y.lololohi)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo > y.hihihilo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo > y.lohihilo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo > y.hilohilo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo = y.hilohilo and x.lolohilo > y.lolohilo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo = y.hilohilo and x.lolohilo = y.lolohilo and
             x.hihilolo > y.hihilolo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo = y.hilohilo and x.lolohilo = y.lolohilo and
             x.hihilolo = y.hihilolo and x.lohilolo > y.lohilolo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo = y.hilohilo and x.lolohilo = y.lolohilo and
             x.hihilolo = y.hihilolo and x.lohilolo = y.lohilolo and
             x.hilololo > y.hilololo)
         or (x.hihihihi = y.hihihihi and x.lohihihi = y.lohihihi and
             x.hilohihi = y.hilohihi and x.lolohihi = y.lolohihi and
             x.hihilohi = y.hihilohi and x.lohilohi = y.lohilohi and
             x.hilolohi = y.hilolohi and x.lololohi = y.lololohi and
             x.hihihilo = y.hihihilo and x.lohihilo = y.lohihilo and
             x.hilohilo = y.hilohilo and x.lolohilo = y.lolohilo and
             x.hihilolo = y.hihilolo and x.lohilolo = y.lohilolo and
             x.hilololo = y.hilololo and x.lolololo > y.lolololo));
  end ">";

  function ">" ( x : hexa_double; y : double_float ) return boolean is
  begin
    return ((x.hihihihi > y)
         or (x.hihihihi = y and x.lohihihi > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo = 0.0 and x.lolohilo > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo = 0.0 and x.lolohilo = 0.0 and
             x.hihilolo > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo = 0.0 and x.lolohilo = 0.0 and
             x.hihilolo = 0.0 and x.lohilolo > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo = 0.0 and x.lolohilo = 0.0 and
             x.hihilolo = 0.0 and x.lohilolo = 0.0 and
             x.hilololo > 0.0)
         or (x.hihihihi = y and x.lohihihi = 0.0 and
             x.hilohihi = 0.0 and x.lolohihi = 0.0 and
             x.hihilohi = 0.0 and x.lohilohi = 0.0 and
             x.hilolohi = 0.0 and x.lololohi = 0.0 and
             x.hihihilo = 0.0 and x.lohihilo = 0.0 and
             x.hilohilo = 0.0 and x.lolohilo = 0.0 and
             x.hihilolo = 0.0 and x.lohilolo = 0.0 and
             x.hilololo = 0.0 and x.lolololo > 0.0));
  end ">";

  function ">" ( x : double_float; y : hexa_double ) return boolean is
  begin
    return (y < x);
  end ">";

  function ">=" ( x,y : hexa_double ) return boolean is
  begin
    return (x > y) or equal(x,y);
  end ">=";

  function ">=" ( x : hexa_double; y : double_float ) return boolean is
  begin
    return (x > y) or equal(x,y);
  end ">=";

  function ">=" ( x : double_float; y : hexa_double ) return boolean is
  begin
    return (y <= x);
  end ">=";

  procedure copy ( x : in hexa_double; y : in out hexa_double ) is
  begin
    y.hihihihi := x.hihihihi; y.lohihihi := x.lohihihi;
    y.hilohihi := x.hilohihi; y.lolohihi := x.lolohihi;
    y.hihilohi := x.hihilohi; y.lohilohi := x.lohilohi;
    y.hilolohi := x.hilolohi; y.lololohi := x.lololohi;
    y.hihihilo := x.hihihilo; y.lohihilo := x.lohihilo;
    y.hilohilo := x.hilohilo; y.lolohilo := x.lolohilo;
    y.hihilolo := x.hihilolo; y.lohilolo := x.lohilolo;
    y.hilololo := x.hilololo; y.lolololo := x.lolololo;
  end copy;

-- ARITHMETICAL OPERATIONS :

  function "+" ( x,y : hexa_double ) return hexa_double is

   -- ALGORITHM : baileyADouble_Double_Basics.fast<16,16,16> generated by CAMPARY.

    res : hexa_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8 : double_float;
    f9,f10,f11,f12,f13,f14,f15,f16,e : double_float;

  begin
    f16 := 0.0;
    Double_Double_Basics.two_sum(x.lolololo,y.lolololo,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.hilololo,y.hilololo,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.lohilolo,y.lohilolo,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.hihilolo,y.hihilolo,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.lolohilo,y.lolohilo,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.hilohilo,y.hilohilo,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.lohihilo,y.lohihilo,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.hihihilo,y.hihihilo,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.lololohi,y.lololohi,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.hilolohi,y.hilolohi,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.lohilohi,y.lohilohi,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.hihilohi,y.hihilohi,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.lolohihi,y.lolohihi,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.hilohihi,y.hilohihi,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.lohihihi,y.lohihihi,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(x.hihihihi,y.hihihihi,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                res.hihihihi,res.lohihihi,res.hilohihi,res.lolohihi,
                res.hihilohi,res.lohilohi,res.hilolohi,res.lololohi,
                res.hihihilo,res.lohihilo,res.hilohilo,res.lolohilo,
                res.hihilolo,res.lohilolo,res.hilololo,res.lolololo);
    return res;
  end "+";

  function "+" ( x : hexa_double; y : double_float ) return hexa_double is

    res : hexa_double;

  begin
    renorm_add1(x.hihihihi,x.lohihihi,x.hilohihi,x.lolohihi,
                x.hihilohi,x.lohilohi,x.hilolohi,x.lololohi,
                x.hihihilo,x.lohihilo,x.hilohilo,x.lolohilo,
                x.hihilolo,x.lohilolo,x.hilololo,x.lolololo,y,
                res.hihihihi,res.lohihihi,res.hilohihi,res.lolohihi,
                res.hihilohi,res.lohilohi,res.hilolohi,res.lololohi,
                res.hihihilo,res.lohihilo,res.hilohilo,res.lolohilo,
                res.hihilolo,res.lohilolo,res.hilololo,res.lolololo);
    return res;
  end "+";

  function "+" ( x : double_float; y : hexa_double ) return hexa_double is

    res : constant hexa_double := y + x;

  begin
    return res;
  end "+";

  function "+" ( x : hexa_double ) return hexa_double is

    res : hexa_double;

  begin
    res.hihihihi := x.hihihihi;
    res.lohihihi := x.lohihihi;
    res.hilohihi := x.hilohihi;
    res.lolohihi := x.lolohihi;
    res.hihilohi := x.hihilohi;
    res.lohilohi := x.lohilohi;
    res.hilolohi := x.hilolohi;
    res.lololohi := x.lololohi;
    res.hihihilo := x.hihihilo;
    res.lohihilo := x.lohihilo;
    res.hilohilo := x.hilohilo;
    res.lolohilo := x.lolohilo;
    res.hihilolo := x.hihilolo;
    res.lohilolo := x.lohilolo;
    res.hilololo := x.hilololo;
    res.lolololo := x.lolololo;
    return res;
  end "+";

  function "-" ( x : hexa_double ) return hexa_double is

    res : hexa_double;

  begin
    res.hihihihi := -x.hihihihi;
    res.lohihihi := -x.lohihihi;
    res.hilohihi := -x.hilohihi;
    res.lolohihi := -x.lolohihi;
    res.hihilohi := -x.hihilohi;
    res.lohilohi := -x.lohilohi;
    res.hilolohi := -x.hilolohi;
    res.lololohi := -x.lololohi;
    res.hihihilo := -x.hihihilo;
    res.lohihilo := -x.lohihilo;
    res.hilohilo := -x.hilohilo;
    res.lolohilo := -x.lolohilo;
    res.hihilolo := -x.hihilolo;
    res.lohilolo := -x.lohilolo;
    res.hilololo := -x.hilololo;
    res.lolololo := -x.lolololo;
    return res;
  end "-";

  function "-" ( x,y : hexa_double ) return hexa_double is

    mny : constant hexa_double := -y;
    res : constant hexa_double := x + mny;

  begin
    return res;
  end "-";

  function "-" ( x : hexa_double; y : double_float ) return hexa_double is

    mny : constant double_float := -y;
    res : constant hexa_double := x + mny;

  begin
    return res;
  end "-";

  function "-" ( x : double_float; y : hexa_double ) return hexa_double is

    mny : constant hexa_double := -y;
    res : constant hexa_double := x + mny;

  begin
    return res;
  end "-";

  function "*" ( x,y : hexa_double ) return hexa_double is

  -- ALGORITHM : baileyMul_fast<16,16,16> generated by CAMPARY.

    res : hexa_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8 : double_float;
    f9,f10,f11,f12,f13,f14,f15,f16,p,e : double_float;

  begin
    f16 :=  x.lohihihi*y.lolololo;
    f16 := f16 + x.hilohihi*y.hilololo;
    f16 := f16 + x.lolohihi*y.lohilolo;
    f16 := f16 + x.hihilohi*y.hihilolo;
    f16 := f16 + x.lohilohi*y.lolohilo;
    f16 := f16 + x.hilolohi*y.hilohilo;
    f16 := f16 + x.lololohi*y.lohihilo;
    f16 := f16 + x.hihihilo*y.hihihilo;
    f16 := f16 + x.lohihilo*y.lololohi;
    f16 := f16 + x.hilohilo*y.hilolohi;
    f16 := f16 + x.lolohilo*y.lohilohi;
    f16 := f16 + x.hihilolo*y.hihilohi;
    f16 := f16 + x.lohilolo*y.lolohihi;
    f16 := f16 + x.hilololo*y.hilohihi;
    f16 := f16 + x.lolololo*y.lohihihi;
    Double_Double_Basics.two_prod(x.hihihihi,y.lolololo,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.hilololo,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.lohilolo,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.hihilolo,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.lolohilo,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.hilohilo,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y.lohihilo,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lololohi,y.hihihilo,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihilo,y.lololohi,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihilo,y.hilolohi,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohilo,y.lohilohi,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohilo,y.hihilohi,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilolo,y.lolohihi,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilolo,y.hilohihi,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilololo,y.lohihihi,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolololo,y.hihihihi,p,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f15,p,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.hilololo,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.lohilolo,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.hihilolo,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.lolohilo,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.hilohilo,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.lohihilo,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y.hihihilo,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lololohi,y.lololohi,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihilo,y.hilolohi,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihilo,y.lohilohi,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohilo,y.hihilohi,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohilo,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilolo,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilolo,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilololo,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f14,p,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.lohilolo,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.hihilolo,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.lolohilo,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.hilohilo,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.lohihilo,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.hihihilo,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y.lololohi,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lololohi,y.hilolohi,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihilo,y.lohilohi,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihilo,y.hihilohi,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohilo,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohilo,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilolo,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilolo,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f13,p,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.hihilolo,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.lolohilo,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.hilohilo,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.lohihilo,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.hihihilo,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.lololohi,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y.hilolohi,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lololohi,y.lohilohi,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihilo,y.hihilohi,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihilo,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohilo,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohilo,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilolo,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f12,p,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.lolohilo,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.hilohilo,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.lohihilo,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.hihihilo,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.lololohi,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.hilolohi,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y.lohilohi,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lololohi,y.hihilohi,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihilo,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihilo,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohilo,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohilo,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f11,p,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.hilohilo,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.lohihilo,p,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f10,p,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.hihihilo,p,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f10,p,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.lololohi,p,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f10,p,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.hilolohi,p,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f10,p,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.lohilohi,p,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f10,p,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y.hihilohi,p,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f10,p,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lololohi,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f10,p,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihilo,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f10,p,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihilo,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f10,p,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohilo,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f10,p,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.lohihilo,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.hihihilo,p,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.lololohi,p,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.hilolohi,p,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.lohilohi,p,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.hihilohi,p,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lololohi,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihilo,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihilo,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f9,p,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.hihihilo,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.lololohi,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.hilolohi,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.lohilohi,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.hihilohi,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lololohi,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihilo,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f8,p,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.lololohi,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.hilolohi,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.lohilohi,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.hihilohi,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lololohi,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f7,p,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.hilolohi,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.lohilohi,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.hihilohi,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f6,p,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.lohilohi,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.hihilohi,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f5,p,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.hihilohi,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.lolohihi,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f4,p,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.lolohihi,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.hilohihi,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f3,p,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.hilohihi,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.lohihihi,p,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f2,p,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f2,p,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.lohihihi,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y.hihihihi,p,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_sum(f1,p,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y.hihihihi,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                res.hihihihi,res.lohihihi,res.hilohihi,res.lolohihi,
                res.hihilohi,res.lohilohi,res.hilolohi,res.lololohi,
                res.hihihilo,res.lohihilo,res.hilohilo,res.lolohilo,
                res.hihilolo,res.lohilolo,res.hilololo,res.lolololo);
    return res;
  end "*";

  function "*" ( x : hexa_double; y : double_float ) return hexa_double is

  -- ALGORITHM : baileyMul_fast<16,1,16> generated by CAMPARY.

    res : hexa_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8 : double_float;
    f9,f10,f11,f12,f13,f14,f15,f16,e : double_float;

  begin
    f16 := 0.0;
    Double_Double_Basics.two_prod(x.lolololo,y,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilololo,y,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilolo,y,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilolo,y,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohilo,y,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohilo,y,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihilo,y,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihilo,y,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lololohi,y,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilolohi,y,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohilohi,y,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihilohi,y,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lolohihi,y,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hilohihi,y,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.lohihihi,y,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;
    Double_Double_Basics.two_prod(x.hihihihi,y,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    Double_Double_Basics.two_sum(f10,e,f10,e);
    Double_Double_Basics.two_sum(f11,e,f11,e);
    Double_Double_Basics.two_sum(f12,e,f12,e);
    Double_Double_Basics.two_sum(f13,e,f13,e);
    Double_Double_Basics.two_sum(f14,e,f14,e);
    Double_Double_Basics.two_sum(f15,e,f15,e);
    f16 := f16 + e;

    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,
                res.hihihihi,res.lohihihi,res.hilohihi,res.lolohihi,
                res.hihilohi,res.lohilohi,res.hilolohi,res.lololohi,
                res.hihihilo,res.lohihilo,res.hilohilo,res.lolohilo,
                res.hihilolo,res.lohilolo,res.hilololo,res.lolololo);
    return res;
  end "*";

  function "*" ( x : double_float; y : hexa_double ) return hexa_double is

    res : constant hexa_double := y*x;

  begin
    return res;
  end "*";

  function Mul_pwr2 ( x : hexa_double; y : double_float )
                    return hexa_double is

    res : hexa_double;

  begin
    res.hihihihi := x.hihihihi*y; res.lohihihi := x.lohihihi*y;
    res.hilohihi := x.hilohihi*y; res.lolohihi := x.lolohihi*y;
    res.hihilohi := x.hihilohi*y; res.lohilohi := x.lohilohi*y;
    res.hilolohi := x.hilolohi*y; res.lololohi := x.lololohi*y;
    res.hihihilo := x.hihihilo*y; res.lohihilo := x.lohihilo*y;
    res.hilohilo := x.hilohilo*y; res.lolohilo := x.lolohilo*y;
    res.hihilolo := x.hihilolo*y; res.lohilolo := x.lohilolo*y;
    res.hilololo := x.hilololo*y; res.lolololo := x.lolololo*y;
    return res;
  end Mul_pwr2;

  procedure Mul_pwr2 ( x : in out hexa_double; y : in double_float ) is
  begin
    x.hihihihi := x.hihihihi*y; x.lohihihi := x.lohihihi*y;
    x.hilohihi := x.hilohihi*y; x.lolohihi := x.lolohihi*y;
    x.hihilohi := x.hihilohi*y; x.lohilohi := x.lohilohi*y;
    x.hilolohi := x.hilolohi*y; x.lololohi := x.lololohi*y;
    x.hihihilo := x.hihihilo*y; x.lohihilo := x.lohihilo*y;
    x.hilohilo := x.hilohilo*y; x.lolohilo := x.lolohilo*y;
    x.hihilolo := x.hihilolo*y; x.lohilolo := x.lohilolo*y;
    x.hilololo := x.hilololo*y; x.lolololo := x.lolololo*y;
  end Mul_pwr2;

  function "/" ( x,y : hexa_double ) return hexa_double is

    res,acc : hexa_double;
    q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16 : double_float;

  begin
    q0 := x.hihihihi/y.hihihihi;   acc := q0*y; res := x - acc;
    q1 := res.hihihihi/y.hihihihi; acc := q1*y; res := res - acc;
    q2 := res.hihihihi/y.hihihihi; acc := q2*y; res := res - acc;
    q3 := res.hihihihi/y.hihihihi; acc := q3*y; res := res - acc;
    q4 := res.hihihihi/y.hihihihi; acc := q4*y; res := res - acc;
    q5 := res.hihihihi/y.hihihihi; acc := q5*y; res := res - acc;
    q6 := res.hihihihi/y.hihihihi; acc := q6*y; res := res - acc;
    q7 := res.hihihihi/y.hihihihi; acc := q7*y; res := res - acc;
    q8 := res.hihihihi/y.hihihihi; acc := q8*y; res := res - acc;
    q9 := res.hihihihi/y.hihihihi; acc := q9*y; res := res - acc;
    q10 := res.hihihihi/y.hihihihi; acc := q10*y; res := res - acc;
    q11 := res.hihihihi/y.hihihihi; acc := q11*y; res := res - acc;
    q12 := res.hihihihi/y.hihihihi; acc := q12*y; res := res - acc;
    q13 := res.hihihihi/y.hihihihi; acc := q13*y; res := res - acc;
    q14 := res.hihihihi/y.hihihihi; acc := q14*y; res := res - acc;
    q15 := res.hihihihi/y.hihihihi; acc := q15*y; res := res - acc;
    q16 := res.hihihihi/y.hihihihi;
    fast_renorm(q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,
                res.hihihihi,res.lohihihi,res.hilohihi,res.lolohihi,
                res.hihilohi,res.lohilohi,res.hilolohi,res.lololohi,
                res.hihihilo,res.lohihilo,res.hilohilo,res.lolohilo,
                res.hihilolo,res.lohilolo,res.hilololo,res.lolololo);
    return res;
  end "/";

  function "/" ( x : hexa_double; y : double_float ) return hexa_double is

    hdy : constant hexa_double := create(y);
    res : constant hexa_double := x/hdy;

  begin
    return res;
  end "/";

  function "/" ( x : double_float; y : hexa_double ) return hexa_double is

    hdx : constant hexa_double := create(x);
    res : constant hexa_double := hdx/y;

  begin
    return res;
  end "/";

  function sqr ( x : hexa_double ) return hexa_double is
  begin
    return x*x;
  end sqr;

  function "**" ( x : hexa_double; n : integer ) return hexa_double is

    res,acc : hexa_double;
    absn : natural;

  begin
    if n = 0 then
      res.hihihihi := 1.0; res.lohihihi := 0.0;
      res.hilohihi := 0.0; res.lolohihi := 0.0;
      res.hihilohi := 0.0; res.lohilohi := 0.0;
      res.hilolohi := 0.0; res.lololohi := 0.0;
      res.hihihilo := 0.0; res.lohihilo := 0.0;
      res.hilohilo := 0.0; res.lolohilo := 0.0;
      res.hihilolo := 0.0; res.lohilolo := 0.0;
      res.hilololo := 0.0; res.lolololo := 0.0;
    else
      if n > 0
       then absn := n;
       else absn := -n;
      end if;
      res.hihihihi := x.hihihihi; res.lohihihi := x.lohihihi;
      res.hilohihi := x.hilohihi; res.lolohihi := x.lolohihi;
      res.hihilohi := x.hihilohi; res.lohilohi := x.lohilohi;
      res.hilolohi := x.hilolohi; res.lololohi := x.lololohi;
      res.hihihilo := x.hihihilo; res.lohihilo := x.lohihilo;
      res.hilohilo := x.hilohilo; res.lolohilo := x.lolohilo;
      res.hihilolo := x.hihilolo; res.lohilolo := x.lohilolo;
      res.hilololo := x.hilololo; res.lolololo := x.lolololo;
      acc.hihihihi := 1.0; acc.lohihihi := 0.0;
      acc.hilohihi := 0.0; acc.lolohihi := 0.0;
      acc.hihilohi := 0.0; acc.lohilohi := 0.0;
      acc.hilolohi := 0.0; acc.lololohi := 0.0;
      acc.hihihilo := 0.0; acc.lohihilo := 0.0;
      acc.hilohilo := 0.0; acc.lolohilo := 0.0;
      acc.hihilolo := 0.0; acc.lohilolo := 0.0;
      acc.hilololo := 0.0; acc.lolololo := 0.0;
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
        acc.hihihihi := res.hihihihi; acc.lohihihi := res.lohihihi;
        acc.hilohihi := res.hilohihi; acc.lolohihi := res.lolohihi;
        acc.hihilohi := res.hihilohi; acc.lohilohi := res.lohilohi;
        acc.hilolohi := res.hilolohi; acc.lololohi := res.lololohi;
        acc.hihihilo := res.hihihilo; acc.lohihilo := res.lohihilo;
        acc.hilohilo := res.hilohilo; acc.lolohilo := res.lolohilo;
        acc.hihilolo := res.hihilolo; acc.lohilolo := res.lohilolo;
        acc.hilololo := res.hilololo; acc.lolololo := res.lolololo;
      end if;
      if n < 0 then
        res := 1.0/acc;          -- compute reciprocal
      else 
        res.hihihihi := acc.hihihihi; res.lohihihi := acc.lohihihi;
        res.hilohihi := acc.hilohihi; res.lolohihi := acc.lolohihi;
        res.hihilohi := acc.hihilohi; res.lohilohi := acc.lohilohi;
        res.hilolohi := acc.hilolohi; res.lololohi := acc.lololohi;
        res.hihihilo := acc.hihihilo; res.lohihilo := acc.lohihilo;
        res.hilohilo := acc.hilohilo; res.lolohilo := acc.lolohilo;
        res.hihilolo := acc.hihilolo; res.lohilolo := acc.lohilolo;
        res.hilololo := acc.hilololo; res.lolololo := acc.lolololo;
      end if;
    end if;
    return res;
  end "**";

  function ldexp ( x : hexa_double; n : integer ) return hexa_double is

    res : hexa_double;
    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

  begin
    res.hihihihi := C_ldexp(x.hihihihi,n);
    res.lohihihi := C_ldexp(x.lohihihi,n);
    res.hilohihi := C_ldexp(x.hilohihi,n);
    res.lolohihi := C_ldexp(x.lolohihi,n);
    res.hihilohi := C_ldexp(x.hihilohi,n);
    res.lohilohi := C_ldexp(x.lohilohi,n);
    res.hilolohi := C_ldexp(x.hilolohi,n);
    res.lololohi := C_ldexp(x.lololohi,n);
    res.hihihilo := C_ldexp(x.hihihilo,n);
    res.lohihilo := C_ldexp(x.lohihilo,n);
    res.hilohilo := C_ldexp(x.hilohilo,n);
    res.lolohilo := C_ldexp(x.lolohilo,n);
    res.hihilolo := C_ldexp(x.hihilolo,n);
    res.lohilolo := C_ldexp(x.lohilolo,n);
    res.hilololo := C_ldexp(x.hilololo,n);
    res.lolololo := C_ldexp(x.lolololo,n);
    return res;
  end ldexp;

  function "**" ( x,y : hexa_double ) return hexa_double is
  begin
    return exp(y*log(x));
  end "**";

  function "**" ( x : hexa_double; y : double_float ) return hexa_double is

    hd_y : constant hexa_double := create(y);

  begin
    return x**hd_y;
  end "**";

  function exp ( x : hexa_double ) return hexa_double is
  begin
    return x;
  end exp;

  function log ( x : hexa_double ) return hexa_double is
  begin
    return x;
  end log;

  function log10 ( x : hexa_double ) return hexa_double is
  begin
    return x;
  end log10;

-- ARITHMETICAL OPERATIONS AS PROCEDURES :

  procedure Add ( x : in out hexa_double; y : in hexa_double ) is
  begin
    x := x + y;
  end Add;

  procedure Sub ( x : in out hexa_double; y : in hexa_double ) is
  begin
    x := x - y;
  end Sub;

  procedure Min ( x : in out hexa_double ) is
  begin
    x := -x;
  end Min;

  procedure Mul ( x : in out hexa_double; y : in hexa_double ) is
  begin
    x := x * y;
  end Mul;

  procedure Div ( x : in out hexa_double; y : in hexa_double ) is
  begin
    x := x / y;
  end Div;

-- DESTRUCTOR :

  procedure clear ( x : in out hexa_double ) is
  begin
    x.hihihihi := 0.0; x.lohihihi := 0.0;
    x.hilohihi := 0.0; x.lolohihi := 0.0;
    x.hihilohi := 0.0; x.lohilohi := 0.0;
    x.hilolohi := 0.0; x.lololohi := 0.0;
    x.hihihilo := 0.0; x.lohihilo := 0.0;
    x.hilohilo := 0.0; x.lolohilo := 0.0;
    x.hihilolo := 0.0; x.lohilolo := 0.0;
    x.hilololo := 0.0; x.lolololo := 0.0;
  end clear;

end Hexa_Double_Numbers;
