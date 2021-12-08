with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with Double_Double_Basics;
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

  function to_double_double ( x : hexa_double ) return double_double is

    res : constant double_double := create(x.hihihihi,x.lohihihi);

  begin
    return res;
  end to_double_double;

  function to_triple_double ( x : hexa_double ) return triple_double is

    res : constant triple_double
        := create(x.hihihihi,x.lohihihi,x.hilohihi);

  begin
    return res;
  end to_triple_double;

  function to_quad_double ( x : hexa_double ) return quad_double is

    res : constant quad_double
        := create(x.hihihihi,x.lohihihi,x.hilohihi,x.lolohihi);

  begin
    return res;
  end to_quad_double;

  function to_penta_double ( x : hexa_double ) return penta_double is

    res : constant penta_double
        := create(x.hihihihi,x.lohihihi,x.hilohihi,x.lolohihi,
                  x.hihilohi);

  begin
    return res;
  end to_penta_double;

  function to_octo_double ( x : hexa_double ) return octo_double is

    res : constant octo_double
        := create(x.hihihihi,x.lohihihi,x.hilohihi,x.lolohihi,
                  x.hihilohi,x.lohilohi,x.hilolohi,x.lololohi);

  begin
    return res;
  end to_octo_double;

  function to_deca_double ( x : hexa_double ) return deca_double is

    res : constant deca_double
        := create(x.hihihihi,x.lohihihi,x.hilohihi,x.lolohihi,
                  x.hihilohi,x.lohilohi,x.hilolohi,x.lololohi,
                  x.hihihilo,x.lohihilo);

  begin
    return res;
  end to_deca_double;

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

  -- ALGORITHM : baileyADouble_Double_Basics.fast<16,16,16>
  --   generated by CAMPARY.

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

  -- Strategy:  We first reduce the size of x by noting that
  --    exp(kr + m * log(2)) = 2^m * exp(r)^k
  -- where m and k are integers.  By choosing m appropriately
  -- we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
  -- evaluated using the familiar Taylor series.  Reducing the
  -- argument substantially speeds up the convergence.

    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

    res : hexa_double;
    k : constant double_float := C_ldexp(1.0,137);
    inv_k : constant double_float := 1.0/k;
    e0  : constant double_float :=  2.71828182845904509E+00;
    e1  : constant double_float :=  1.44564689172925016E-16;
    e2  : constant double_float := -2.12771710803817676E-33;
    e3  : constant double_float :=  1.51563015984121914E-49;
    e4  : constant double_float := -9.33538137882084738E-66;
    e5  : constant double_float := -2.02185260335762112E-82;
    e6  : constant double_float :=  8.68351053192660636E-99;
    e7  : constant double_float := -7.14388345698345838E-115;
    e8  : constant double_float := -4.29981073684482334E-131;
    e9  : constant double_float := -5.09120833096307459E-149;
    e10 : constant double_float := -2.64514695832580870E-166;
    e11 : constant double_float := -1.03856280207021200E-182;
    e12 : constant double_float := -6.60048905371551111E-199;
    e13 : constant double_float := -3.54999918598759126E-215;
    e14 : constant double_float := -2.95381612024769761E-232;
    e15 : constant double_float :=  1.59387336409917663E-248;
    exp1 : constant hexa_double
         := Create(e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15);
    L0  : constant double_float :=  6.93147180559945286E-01;
    L1  : constant double_float :=  2.31904681384629956E-17;
    L2  : constant double_float :=  5.70770843841621207E-34;
    L3  : constant double_float := -3.58243221060181142E-50;
    L4  : constant double_float := -1.35216967579886296E-66;
    L5  : constant double_float :=  6.08063874024081391E-83;
    L6  : constant double_float :=  2.89550243323471469E-99;
    L7  : constant double_float :=  2.35138671214564106E-116;
    L8  : constant double_float :=  4.45977441701428101E-133;
    L9  : constant double_float := -3.06993326323252717E-149;
    L10 : constant double_float := -2.01514744619668319E-165;
    L11 : constant double_float :=  1.61853434886374103E-182;
    L12 : constant double_float := -1.30949780474544622E-198;
    L13 : constant double_float :=  6.66518827858982421E-215;
    L14 : constant double_float := -2.93171211597727003E-231;
    L15 : constant double_float :=  7.85979956399546348E-248;
    log2 : constant hexa_double
         := Create(L0,L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15);
    hd_eps : constant double_float := 5.327993384780537E-256; -- 2.0^(-53*16)
    tol : constant double_float := inv_k*hd_eps;
    m : constant double_float := double_float'floor(x.hihihihi/L0 + 0.5);
    i_fac : array(0..14) of hexa_double;
      -- inverse factorials for Taylor expansion
    p,s,t : hexa_double;
    cnt : integer;

  begin
    if x.hihihihi <= -709.0 then
      res := Create(0.0);
    elsif x.hihihihi >= 709.0 then
      res := Create(-1.0);
    elsif is_zero(x) then
      res := Create(1.0);
    elsif is_one(x) then
      res := exp1;
    else
      i_fac(0) := Create(1.66666666666666657E-01,  9.25185853854297066E-18,
                         5.13581318503262866E-34,  2.85094902409834186E-50,
                         1.58259462429329970E-66,  8.78516495269210220E-83,
                         4.87674620280437299E-99,  2.70713795980339238E-115,
                         1.22261156727740753E-131, 9.74717870535265007E-148,
                        -5.82587312349853683E-164,-3.77554461049033400E-180,
                         1.83914509634085812E-196, 8.66044561926341746E-213,
                        -1.68140325688177841E-229,-5.75245515446792921E-246);
      i_fac(1) := Create(4.16666666666666644E-02,  2.31296463463574266E-18,
                         1.28395329625815716E-34,  7.12737256024585466E-51,
                         3.95648656073324926E-67,  2.19629123817302555E-83,
                         1.21918655070109325E-99,  6.76784489950843137E-116,
                        -1.03153360308043797E-132,-5.89409547576532608E-149,
                        -4.08934285046332034E-165,-8.83640772366872816E-182,
                        -5.68848573224037734E-198, 2.71953447459104285E-214,
                        -3.24691718649123437E-232, 2.54478647527963024E-248);
      i_fac(2) := Create(8.33333333333333322E-03,  1.15648231731787138E-19,
                         1.60494162032269652E-36,  2.22730392507682967E-53,
                         3.09100512557285111E-70,  4.28963132455669071E-87,
                         5.95305932959518212E-104, 8.26152941801187828E-121,
                         6.99698728255496442E-137, 1.88696466546395415E-153,
                        -2.58015324267443182E-170, 7.21304190564117080E-187,
                        -3.61752526611418417E-203,-2.27802497696402197E-219,
                         1.37624330020480649E-235,-6.28150291922469969E-252);
      i_fac(3) := Create(1.38888888888888894E-03, -5.30054395437357706E-20,
                        -1.73868675534958776E-36, -1.63335621172300840E-52,
                        -5.35774221765960817E-69, -5.03316742081318367E-85,
                        -1.65098178740773037E-101,-1.55096445613733005E-117,
                        -3.43844534360146001E-134,-6.61258976828610821E-151,
                        -9.68942927329279048E-168, 1.92447249796470790E-184,
                        -3.23357920876208907E-201,-9.66906015862940204E-218,
                         6.62172062935478925E-234,-2.69099039861369768E-250);
      i_fac(4) := Create(1.98412698412698413E-04,  1.72095582934207053E-22,
                         1.49269123913941271E-40,  1.29470326746002471E-58,
                         1.12297607624339191E-76,  9.74026481209866379E-95,
                         8.44833301588918793E-113, 7.48649096649661558E-131,
                         1.19813886009252660E-147,-3.19798992768170296E-164,
                         1.48764734138575271E-181,-2.76634008526703165E-198,
                         1.03846458794745048E-214, 7.26161589010453353E-231,
                         1.03370541212648439E-248, 0.00000000000000000E+00);
      i_fac(5) := Create(2.48015873015873016E-05,  2.15119478667758816E-23,
                         1.86586404892426588E-41,  1.61837908432503088E-59,
                         1.40372009530423989E-77,  1.21753310151233297E-95,
                         1.05604162698614849E-113, 8.85811370812076948E-132,
                         1.42738989191293614E-148, 8.43786203392313070E-165,
                        -2.75041093672864639E-181, 4.59243836924367153E-198,
                        -2.04641686548929214E-214,-1.00869684999106651E-231,
                         1.17660941215217309E-248,-2.77566378900000154E-265);
      i_fac(6) := Create(2.75573192239858925E-06, -1.85839327404647208E-22,
                         8.49175460488199287E-39, -5.72661640789429621E-55,
                         2.61672391582886888E-71, -1.76464992319726308E-87,
                         8.06340311310246955E-104,-5.43774740551394713E-120,
                        -1.97632978408579155E-136, 8.47573041938171586E-153,
                         1.12351901997189503E-169, 5.66711386774374456E-186,
                        -4.36357064022533437E-203, 2.23343831877544446E-219,
                        -1.18464196118406929E-235,-1.86181533903132709E-252);
      i_fac(7) := Create(2.75573192239858883E-07, 2.37677146222502973E-23,
                        -3.26318890334088294E-40, 1.61435111860404415E-56,
                        -9.99798078191450073E-74, -5.23082523455780306E-91,
                        -2.46767621898134459E-107, 8.08497474290649098E-124,
                         2.08472959271617542E-140,-8.13216235851919747E-157,
                        -2.10079199559597373E-173,-6.76162337253620171E-190,
                        -1.78736415467946995E-206,-7.38806966587211950E-223,
                        -5.23239643850046933E-239,-3.75104040824711363E-255);
      i_fac(8) := Create(2.50521083854417202E-08, -1.44881407093591197E-24,
                         2.04267351467144546E-41, -8.49632672007163175E-58,
                        -9.08907343810409139E-75, -2.26065446324988261E-91,
                        -7.19805892199252488E-108,-6.40212493138514344E-125,
                        -4.80964612854587844E-141, 2.76952341964975665E-157,
                         1.58592999339168693E-173,-1.43721185168648149E-190,
                         1.15038993246231282E-206, 3.58005012020971148E-223,
                         6.23805819878531708E-240,-4.36828280297612803E-256);
      i_fac(9) := Create(2.08767569878681002E-09, -1.20734505911325997E-25,
                         1.70222792889287100E-42,  1.41609532150396700E-58,
                         5.13820161376913309E-75,  1.44797661649508516E-91,
                         1.30256332445267049E-107, 3.72847700074999458E-124,
                        -7.08461312493450429E-141,-3.21433824881373091E-157,
                        -1.04576131833441555E-174, 1.10694447546812161E-190,
                        -5.86922957336079147E-207,-2.86550988437943358E-223,
                         1.67917720852385565E-239, 1.87411792142416342E-256);
      i_fac(10) := Create(1.60590438368216133E-10,  1.25852945887520981E-26,
                         -5.31334602762985031E-43,  3.54021472597605528E-59,
                         -1.98567896059148884E-75, -1.02148490610754577E-91,
                          5.19442455358700245E-108, 2.03226501937076278E-124,
                          2.04111228153276298E-140,-3.09982828459628974E-157,
                          3.32392835861557285E-174,-1.33335163123651162E-190,
                          1.66995016156050641E-207, 4.58726118892580877E-224,
                         -1.47581171014520755E-240, 8.51983503178078183E-258);
      i_fac(11) := Create(1.14707455977297245E-11,  2.06555127528307454E-28,
                          6.88907923246664603E-45,  5.72920002655109095E-61,
                         -3.65551458930952487E-78, -1.73752114063892946E-94,
                         -1.29464103159491163E-110,-4.67626589929544217E-127,
                          2.08558412597439078E-143, 1.54340416983905211E-159,
                         -1.03282539198955174E-175, 4.33390566027236715E-192,
                         -1.38117575424250920E-208, 2.77742358833030132E-225,
                         -4.26440138393592775E-242,-3.06998801663533946E-258);
      i_fac(12) := Create(7.64716373181981641E-13,  7.03872877733453001E-30,
                         -7.82753927716258345E-48,  1.92138649443790242E-64,
                         -1.43027453759388186E-80,  7.76151305461206293E-97,
                          3.19012607380048399E-114, 1.78550560934447991E-130,
                          4.52147920813402651E-147, 2.59219293685287978E-163,
                         -7.70983451628364634E-180, 4.35654131810895619E-197,
                         -1.90120616037127669E-213, 1.12129042759357876E-229,
                         -8.81146651499022840E-246, 1.49430480158690228E-262);
      i_fac(13) := Create(4.77947733238738525E-14,  4.39920548583408126E-31,
                         -4.89221204822661465E-49,  1.20086655902368901E-65,
                         -8.93921585996176165E-82,  4.85094565913253933E-98,
                          1.99382879612530249E-115, 1.11594100586529995E-131,
                          1.66496936695943139E-148, 1.13949171616546479E-164,
                         -2.13988072119987192E-181,-1.92763218164551551E-198,
                          1.59664607050423981E-214,-1.87610722392907552E-231,
                          1.09829039805601509E-247, 5.27793172589625328E-264);
      i_fac(14) := Create(2.81145725434552060E-15,  1.65088427308614326E-31,
                         -2.87777179307447918E-50,  4.27110689256293549E-67,
                         -2.93287743014724397E-83, -1.23436334109628973E-99,
                         -2.14851021924804564E-118, 2.19396005239905784E-134,
                          1.10961028311563215E-150,-4.17632175886513988E-167,
                          1.13113853924118729E-183, 9.81506423798373456E-200,
                          5.70764923661908468E-216,-2.33208429934520284E-233,
                         -1.91194635886479302E-250,-1.04355593700700599E-266);
      res := mul_pwr2(x - m*log2,inv_k);
      p := res*res;
      s := res + mul_pwr2(p,0.5);
      cnt := 0;
      loop
        p := p*res;
        t := i_fac(cnt)*p;
        s := s + t;
        exit when abs(t.hihihihi) <= tol;
        cnt := cnt + 1;
        exit when (cnt >= 15);
      end loop;
      for i in 1..137 loop -- 137 times s = mul_pwr2(s,2.0) + sqr(s);
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

  function log ( x : hexa_double ) return hexa_double is

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

    res,acc : hexa_double;

  begin
    if is_one(x) then
      res := Create(0.0);
    elsif x.hihihihi <= 0.0 then
      put_line("hd_log: argument is not positive");
      res := Create(-1.0);
    else
      res := Create(Standard_Mathematical_Functions.LN(x.hihihihi));
      for i in 1..6 loop       -- r = r + x*exp(-r) - 1
        acc := x*exp(-res);    -- acc = x*exp(-res)      
        res := res + acc;      -- res = res + x*exp(-res)
        res := res - 1.0;
      end loop;
    end if;
    return res;
  end log;

  function log10 ( x : hexa_double ) return hexa_double is

    res : hexa_double;
    log10_0  : constant double_float :=  2.30258509299404590E+00;
    log10_1  : constant double_float := -2.17075622338224935E-16;
    log10_2  : constant double_float := -9.98426245446577657E-33;
    log10_3  : constant double_float := -4.02335745445020638E-49;
    log10_4  : constant double_float :=  1.92889952896933719E-65;
    log10_5  : constant double_float := -5.21257011815125513E-82;
    log10_6  : constant double_float := -2.60373698986932938E-98;
    log10_7  : constant double_float :=  8.29741762082190114E-115;
    log10_8  : constant double_float := -4.10600981629892264E-131;
    log10_9  : constant double_float := -5.26148805167540611E-148;
    log10_10 : constant double_float := -3.18682909306959149E-164;
    log10_11 : constant double_float := -1.19106223015530956E-180;
    log10_12 : constant double_float :=  9.10673440282089099E-197;
    log10_13 : constant double_float :=  3.22698486627368573E-213;
    log10_14 : constant double_float := -1.75738641610789914E-229;
    log10_15 : constant double_float :=  2.87909858432172129E-246;
    logten : constant hexa_double
           := create(log10_0,log10_1,log10_2,log10_3,
                     log10_4,log10_5,log10_6,log10_7,
                     log10_8,log10_9,log10_10,log10_11,
                     log10_12,log10_13,log10_14,log10_15);
  begin
    res := log(x)/logten;
    return res;
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
