package body Three_Way_Minima is

  function Minimum ( a,b,c : double_float ) return double_float is
  begin
    if a < b then    -- the minimum is either a or c
      if a < c
       then return a;
       else return c;
      end if;
    else             -- the minimum is either b or c
      if b < c
       then return b;
       else return c;
      end if;
    end if;
  end Minimum;

  function Minimum ( a,b,c : double_double ) return double_double is
  begin
    if a < b then    -- the minimum is either a or c
      if a < c
       then return a;
       else return c;
      end if;
    else             -- the minimum is either b or c
      if b < c
       then return b;
       else return c;
      end if;
    end if;
  end Minimum;

  function Minimum ( a,b,c : quad_double ) return quad_double is
  begin
    if a < b then    -- the minimum is either a or c
      if a < c
       then return a;
       else return c;
      end if;
    else             -- the minimum is either b or c
      if b < c
       then return b;
       else return c;
      end if;
    end if;
  end Minimum;

  procedure Minimum ( a,b,c : in double_float; abcmin : out double_float;
                      anb,bnb,cnb : in out natural32 ) is
  begin
    if a < b then    -- the minimum is either a or c
      if a < c
       then abcmin := a; anb := anb + 1;
       else abcmin := c; cnb := cnb + 1;
      end if;
    else             -- the minimum is either b or c
      if b < c
       then abcmin := b; bnb := bnb + 1;
       else abcmin := c; cnb := cnb + 1;
      end if;
    end if;
  end Minimum;

  procedure Minimum ( a,b,c : in double_double; abcmin : out double_double;
                      anb,bnb,cnb : in out natural32 ) is
  begin
    if a < b then    -- the minimum is either a or c
      if a < c
       then abcmin := a; anb := anb + 1;
       else abcmin := c; cnb := cnb + 1;
      end if;
    else             -- the minimum is either b or c
      if b < c
       then abcmin := b; bnb := bnb + 1;
       else abcmin := c; cnb := cnb + 1;
      end if;
    end if;
  end Minimum;

  procedure Minimum ( a,b,c : in quad_double; abcmin : out quad_double;
                      anb,bnb,cnb : in out natural32 ) is
  begin
    if a < b then    -- the minimum is either a or c
      if a < c
       then abcmin := a; anb := anb + 1;
       else abcmin := c; cnb := cnb + 1;
      end if;
    else             -- the minimum is either b or c
      if b < c
       then abcmin := b; bnb := bnb + 1;
       else abcmin := c; cnb := cnb + 1;
      end if;
    end if;
  end Minimum;

  procedure Bounded_Update
              ( t,step : in out double_float;
                endt,minstep : in double_float ) is

    t_backup : constant double_float := t;
    endt_bound : constant double_float := endt - minstep;

  begin
    t := t + step;
    if t >= endt_bound then
      t := endt;
      step := endt - t_backup;
    end if;
  end Bounded_Update;

  procedure Bounded_Update
              ( t,step : in out double_double;
                endt,minstep : in double_float ) is

    t_backup : constant double_double := t;
    endt_bound : constant double_float := endt - minstep;
    dd_endt : constant double_double := create(endt);

  begin
    t := t + step;
    if t >= endt_bound then
      t := dd_endt;
      step := dd_endt - t_backup;
    end if;
  end Bounded_Update;

  procedure Bounded_Update
              ( t,step : in out quad_double;
                endt,minstep : in double_float ) is

    t_backup : constant quad_double := t;
    endt_bound : constant double_float := endt - minstep;
    qd_endt : constant quad_double := create(endt);

  begin
    t := t + step;
    if t >= endt_bound then
      t := qd_endt;
      step := qd_endt - t_backup;
    end if;
  end Bounded_Update;

end Three_Way_Minima;
