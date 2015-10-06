with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Mathematical_Functions;
with Characters_and_Numbers;

package body Double_Double_Numbers_io is

  procedure scan_characters
              ( file : in file_type; s : in out string;
                cnt : out natural ) is

    c : character;

  begin
    cnt := 0;
    while not end_of_file(file) loop
      get(file,c);
      cnt := cnt + 1;
      s(cnt) := c;
      exit when end_of_line(file);
      exit when ((c = ASCII.CR) or (c = ASCII.LF));
      exit when (cnt >= s'last);
    end loop;
  end scan_characters;

  procedure scan_for_integer 
              ( s : in string; p : in out integer; i : out integer32 ) is

    sign : integer := 0;
    done : boolean := false;

  begin
    i := 0;
    while s(p) = ' ' loop      -- skip spaces
      p := p + 1;
      exit when (p > s'last);
    end loop;
    if p <= s'last then
      loop
        case s(p) is
          when '-' | '+' =>
            if s(p) = '-'
             then sign := -1;
             else sign := +1;
            end if;
          when '0' => i := 10*i;
          when '1' => i := 10*i + 1;
          when '2' => i := 10*i + 2;
          when '3' => i := 10*i + 3;
          when '4' => i := 10*i + 4;
          when '5' => i := 10*i + 5;
          when '6' => i := 10*i + 6;
          when '7' => i := 10*i + 7;
          when '8' => i := 10*i + 8;
          when '9' => i := 10*i + 9;
          when others => done := true; -- still adjust for sign
        end case;
        exit when done;
        p := p + 1;
        exit when (p > s'last);
      end loop;
      if sign = -1
       then i := -i;
      end if;
    end if;
  end scan_for_integer;

  procedure read ( s : in string; ps : in out integer;
                   x : out double_double; fail : out boolean ) is

    ch : character;
    done : boolean := false;
    d : natural32;
    cd : double_float;
    nd : natural32 := 0;
    point : integer32 := -1;
    sign : integer32 := 0;
    e : integer32 := 0;
    acc : double_double;

  begin
    x := Create(0.0,0.0); fail := false;
    while s(ps) = ' ' loop  -- skip leading spaces
      ps := ps + 1;
    end loop;
    while not done and ps <= s'last loop
      ch := s(ps);
      if (ch >= '0' and ch <= '9') then
        d := Characters_and_Numbers.Convert(ch);
        cd := double_float(d);
        x := 10.0*x + cd;
        nd := nd + 1;
      else
        case ch is
          when '.' =>
            if point >= 0 
             then fail := true; return;
             else point := integer32(nd);
            end if;
          when '-' | '+' =>
            if (sign /= 0 and nd > 0)
             then fail := true; return;
            end if;
            if ch = '-'
             then sign := -1;
             else sign := +1;
            end if;
          when 'E' | 'e' =>
            ps := ps + 1;
            scan_for_integer(s,ps,e);
            done := true;
          when others => done := true;
        end case;
      end if;
      ps := ps + 1;
    end loop;
    if point >= 0
     then e := e - (integer32(nd) - point);
    end if;
    if e /= 0 then
      acc := Create(10.0);
      acc := acc**integer(e);
      Mul(x,acc);
    end if;
    if sign = -1
     then x := -x;
    end if;
  end read;

  procedure read ( s : in string;
                   x : out double_double; fail : out boolean ) is

    ps : integer := s'first;

  begin
    read(s,ps,x,fail);   
  end read;

  procedure to_digits ( x : in double_double; precision : in natural32;
                        s : out string; expn : out integer32 ) is

    Dp1 : constant natural := natural(precision + 1);
    acc,tmp : double_double;
    d,e : integer32;    -- e is exponent
    ind : integer;
    log10x : double_float;
    v0 : constant natural := character'pos('0');

  begin
    s(s'first) := ' '; expn := 0;
    if hi_part(x) = 0.0 then          -- x equals zero
      for k in 1..integer(precision) loop
        s(k) := '0';
      end loop;
    else
      acc := abs(x);
      log10x := Standard_Mathematical_Functions.LOG10(hi_part(acc));
      e := integer32(double_float'floor(log10x));
      if (e < -300) then
        tmp := Create(10.0)**integer(300);
        Mul(acc,tmp);
        tmp := Create(10.0)**(integer(e)+300);
        Div(acc,tmp);
      elsif (e > 300) then
        acc := ldexp(acc,-53);
        tmp := Create(10.0)**integer(e);
        Div(acc,tmp);
        acc := ldexp(acc,+53);
      else
        tmp := Create(10.0)**integer(e);
        Div(acc,tmp);
      end if;
      if acc > 10.0 then -- fix exponent if we are off by one
        Div(acc,10.0); e := e+1;
      elsif acc < 1.0 then
        Mul(acc,10.0); e := e-1;
      end if;
      if (acc >= 10.0 or acc < 1.0)
       then put_line("to_digits: cannot compute exponent"); return;
      end if;
      for k in s'first..s'first+Dp1-1 loop -- extract the digits
        d := integer32(hi_part(acc));
        Sub(acc,double_float(d));         
        Mul(acc,10.0);
        s(k) := character'val(integer(d) + v0);
      end loop;
      -- changed s'first into s'first+1 in analogy with quad_double_numbers_io
      for k in reverse s'first+1..s'first+Dp1-1 loop -- fix out of range digits
        if s(k) < '0' then
          d := character'pos(s(k-1)) - 1; s(k-1) := character'val(d);
          d := character'pos(s(k)) + 10;  s(k) := character'val(d);
        elsif s(k) > '9' then
          d := character'pos(s(k-1)) + 1; s(k-1) := character'val(d);
          d := character'pos(s(k)) - 10;  s(k) := character'val(d);
        end if;
      end loop;
      if s(s'first) <= '0'
       then put_line("to_digits: nonpositive leading digit"); return;
      end if;
      if s(Dp1) >= '5' then       -- round, handle carry
        d := character'pos(s(Dp1-1)) + 1; s(Dp1-1) := character'val(d);
        ind := Dp1-1;
        while (ind > 0 and s(ind) > '9') loop
          while (s(ind) > '9') loop
            d := character'pos(s(ind)) - 10; s(ind) := character'val(d);
            ind := ind - 1;
            d := character'pos(s(ind)) + 1;  s(ind) := character'val(d);
            exit when (ind <= 1);
          end loop;
          exit when (ind <= 1);
        end loop;
      end if;
      if s(s'first) > '9' then   -- if first digit is 10, shift everything
        e := e + 1;
        for i in reverse s'first+2..integer(precision) loop
          s(i) := s(i-1);
        end loop;
        s(s'first) := '1';
        s(s'first+1) := '0';
      end if;
      expn := e;
    end if;
  end to_digits;

  procedure to_string ( x : in double_double; precision,width : in natural32;
                        fixed,showpos,uppercase : in boolean;
                        fill : in character;
                        s : out string; endstring : out natural ) is

  -- Note : nan and inf are ignored...

    acc : double_double;
    cnt : natural32 := 0;      -- counts #characters written to string s
    d,e : integer32 := 0;
    ind : natural32;
    off : integer32;
    ps : integer := s'first; -- current position in string s

  begin
    if is_negative(x) then
      s(ps) := '-'; ps := ps + 1; cnt := cnt + 1;
    elsif showpos then
      s(ps) := '+'; ps := ps + 1; cnt := cnt + 1;
    end if;
    if is_zero(x) then
      s(ps) := '0'; ps := ps + 1; cnt := cnt + 1;
      if precision > 0 then
        s(ps) := '.'; ps := ps + 1; cnt := cnt + 1;
        for i in 1..precision loop
          s(ps) := '0'; ps := ps + 1; cnt := cnt + 1;
        end loop;
      end if;
    else         -- nonzero case
      if not fixed then
        off := 1;
      else
        acc := floor(log10(abs(x)));
        off := integer32(to_int(acc) + 1);
      end if;
      d := integer32(precision) + off;
      if fixed and (d <= 0) then
        s(ps) := '0'; ps := ps + 1; cnt := cnt + 1;
        if precision > 0 then
          s(ps) := '.'; ps := ps + 1; cnt := cnt + 1;
          for i in 1..precision loop
            s(ps) := '0'; ps := ps + 1; cnt := cnt + 1;
          end loop;
        end if;
      else
        declare
          t : string(1..integer(d+1));
        begin
          ind := natural32(t'first);
          to_digits(x,natural32(d),t,e);
          if fixed then
            if off > 0 then
              for i in 0..integer(off-1) loop
                s(ps) := t(t'first+i); ps := ps + 1; cnt := cnt + 1;
                ind := ind + 1;
              end loop;
              if precision > 0 then
                s(ps) := '.'; ps := ps + 1; cnt := cnt + 1;
                for j in 1..precision loop
                  s(ps) := t(integer(ind)); ps := ps + 1; cnt := cnt + 1;
                  ind := ind + 1;
                end loop;
              end if;
            else
              s(ps) := '0'; ps := ps + 1; cnt := cnt + 1;
              s(ps) := '.'; ps := ps + 1; cnt := cnt + 1;
              if off < 0 then
                for i in 0..(-off) loop
                  s(ps) := '0'; ps := ps + 1; cnt := cnt + 1;
                end loop;
                for j in 1..d loop
                  s(ps) := t(integer(ind)); ps := ps + 1; cnt := cnt + 1;
                  ind := ind + 1;
                end loop;
              end if;
            end if;
          else
            s(ps) := t(t'first); ps := ps + 1; cnt := cnt + 1;
            if precision > 0 then
              s(ps) := '.'; ps := ps + 1; cnt := cnt + 1;
            end if;
            for i in 1..integer(precision) loop
              s(ps) := t(t'first+i); ps := ps + 1; cnt := cnt + 1;
            end loop;
          end if;
        end;
      end if;
    end if;
    if not fixed then  -- fill in the exponent part
      if uppercase
       then s(ps) := 'E'; ps := ps + 1; cnt := cnt + 1;
       else s(ps) := 'e'; ps := ps + 1; cnt := cnt + 1;
      end if;
      if e < 0
       then s(ps) := '-'; ps := ps + 1; cnt := cnt + 1;
       else s(ps) := '+'; ps := ps + 1; cnt := cnt + 1;
      end if;
      if e < 0
       then e := -e;
      end if;
      if (e >= 100) then
        ind := natural32(e/100);
        s(ps) := Characters_and_Numbers.Convert_Decimal(ind);
        ps := ps + 1; cnt := cnt + 1;
        e := e - integer32(100*ind);
      end if;
      ind := natural32(e/10);
      s(ps) := Characters_and_Numbers.Convert_Decimal(ind);
      ps := ps + 1; cnt := cnt + 1;
      e := e - integer32(10*ind);
      s(ps) := Characters_and_Numbers.Convert_Decimal(natural32(e));
      ps := ps + 1; cnt := cnt + 1;
    end if;
    if cnt >= width then
      endstring := ps-1;
    else                    -- fill in the blanks
      d := integer32(width) - integer32(cnt);
      for i in 0..integer(cnt-1) loop
        s(s'first+i+integer(d)) := s(s'first+i);
      end loop;
      for i in 0..integer(d-1) loop
        s(s'first+i) := fill;
      end loop;
      endstring := integer(width);
    end if;
  end to_string;

  procedure write ( x : in double_double; precision : in natural32 ) is

    s : string(1..integer(precision)+10);
    ends : integer;

  begin
    to_string(x,precision,0,false,false,true,' ',s,ends);
    put(s(1..ends));
  end write;

  procedure get ( x : in out double_double ) is

    temp : string(1..80);
    cnt : natural;
    fail : boolean;

  begin
    get_line(temp,cnt);
    read(temp(1..cnt),x,fail);
  end get;

  procedure get ( file : in file_type; x : in out double_double ) is

    s : string(1..256);
    ind : natural := 0;
    fail : boolean;

  begin
    scan_characters(file,s,ind);
    read(s(1..ind),x,fail);
  end get;

  procedure get ( x,y : in out double_double ) is

    temp : string(1..160);
    cnt : natural;
    ps : integer := 1;
    fail : boolean;

  begin
    get_line(temp,cnt);
    read(temp(1..cnt),ps,x,fail);
    read(temp(ps..cnt),ps,y,fail);
  end get;

  procedure get ( file : in file_type; x,y : in out double_double ) is

    s : string(1..512);
    ind : natural := 0;
    ps : integer := 1;
    fail : boolean;

  begin
    scan_characters(file,s,ind);
    read(s(1..ind),ps,x,fail);
    read(s(ps..ind),ps,y,fail);
  end get;

  procedure put ( x : in double_double ) is
  begin
    put(standard_output,x,32);
  end put;

  procedure put ( file : in file_type; x : in double_double ) is
  begin
    put(file,x,32);
  end put;

  procedure put ( x : in double_double; precision : in natural32 ) is
  begin
    put(standard_output,x,precision);
  end put;

  procedure put ( file : in file_type; x : in double_double;
                  precision : in natural32 ) is

    s : string(1..integer(precision)+10);
    ends : integer;

  begin
    to_string(x,precision,0,false,false,true,' ',s,ends);
    put(file,s(1..ends));
  end put;

end Double_Double_Numbers_io;
