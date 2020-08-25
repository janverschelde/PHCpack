with Characters_and_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers_io;

package body Octo_Double_Numbers_io is

  procedure read ( s : in string; ps : in out integer;
                   x : out octo_double; fail : out boolean ) is

    ch : character;
    done : boolean := false;
    d : natural32;
    cd : double_float;
    nd : natural32 := 0;
    point : integer32 := -1;
    sign : integer32 := 0;
    e : integer32 := 0;
    acc : octo_double;

  begin
    x := Create(0.0); fail := false;
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
            Double_Double_Numbers_io.scan_for_integer(s,ps,e);
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
      x := x*acc; -- Mul(x,acc);
    end if;
    if sign = -1
     then x := -x;
    end if;
  end read;

  procedure read ( s : in string; x : out octo_double; fail : out boolean ) is

    ps : integer := s'first;

  begin
    read(s,ps,x,fail);
  end read;

  procedure to_digits ( x : in octo_double; precision : in natural32;
                        s : out string; expn : out integer32 ) is
  begin
    null;
  end to_digits;

  procedure to_string ( x : in octo_double; precision,width : in natural32;
                        fixed,showpos,uppercase : in boolean;
                        fill : in character;
                        s : out string; endstring : out natural ) is
  begin
    null;
  end to_string;

  procedure write ( x : in octo_double; precision : in natural32 ) is
  begin
    null;
  end write;

  procedure get ( x : in out octo_double ) is

    temp : string(1..320);
    cnt : natural;
    fail : boolean;

  begin
    get_line(temp,cnt);
    read(temp(1..cnt),x,fail);
  end get;

  procedure get ( file : in file_type; x : in out octo_double ) is

    s : string(1..1024);
    ind : natural := 0;
    fail : boolean;

  begin
    Double_Double_Numbers_io.scan_characters(file,s,ind);
    read(s(1..ind),x,fail);
  end get;

  procedure get ( x,y : in out octo_double ) is

    temp : string(1..640);
    cnt : natural;
    ps : integer := 1;
    fail : boolean;

  begin
    get_line(temp,cnt);
    read(temp(1..cnt),ps,x,fail);
    read(temp(ps..cnt),ps,y,fail);
  end get;

  procedure get ( file : in file_type; x,y : in out octo_double ) is

    s : string(1..2048);
    ind : natural := 0;
    ps : integer := 1;
    fail : boolean;

  begin
    Double_Double_Numbers_io.scan_characters(file,s,ind);
    read(s(1..ind),ps,x,fail);
    read(s(ps..ind),ps,y,fail);
  end get;

  procedure put ( x : in octo_double ) is
  begin
    put(standard_output,x,128);
  end put;

  procedure put ( x : in octo_double; precision : in natural32 ) is
  begin
    put(standard_output,x,precision);
  end put;

  procedure put ( file : in file_type; x : in octo_double ) is
  begin
    put(file,x,128);
  end put;

  procedure put ( file : in file_type; x : in octo_double;
                  precision : in natural32 ) is

    s : string(1..integer(precision)+10);
    ends : integer;

  begin
    to_string(x,precision,0,false,false,true,' ',s,ends);
    put(file,s(1..ends));
  end put;

end Octo_Double_Numbers_io;
