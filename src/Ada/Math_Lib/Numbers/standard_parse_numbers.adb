with Characters_and_Numbers;             use Characters_and_Numbers;

package body Standard_Parse_Numbers is

  function Is_to_Skip ( ch : character ) return boolean is
  begin
    return ((ch = ' ') or (ch = ASCII.CR) or (ch = ASCII.LF));
  end Is_to_Skip;

  procedure Skip_Spaces_and_CR
              ( file : in file_type; ch : in out character ) is
  begin
    while Is_to_Skip(ch) loop
      exit when end_of_file(file);
      get(file,ch);
    end loop;
  end Skip_Spaces_and_CR;

  procedure Skip_Spaces_and_CR
              ( s : in string; p : in out integer ) is
  begin
    if ((p >= s'first) and (p <= s'last)) then
      while Is_to_Skip(s(p)) loop
        p := p + 1;
        exit when (p > s'last);
      end loop;
    end if;
  end Skip_Spaces_and_CR;

  procedure Parse ( file : in file_type;
                    char : in out character; i : out integer32;
                    ni : out natural32; sign : out character ) is

    res : integer32 := 0;
    min : boolean := false;
    cnt,temp : natural32 := 0;

  begin
    i := 0; ni := 0; sign := '+';
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    if (char = '+') or (char = '-') then
      min := (char = '-');
      sign := char;
      if end_of_file(file) then return; end if;
      get(file,char); 
    end if;
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    loop
      temp := Convert(char);
      if temp < 10 then
        if cnt < 9 then
          res := res*10 + integer32(temp);
          cnt := cnt + 1;
        else
          null;
        end if;
        if end_of_file(file) then return; end if;
        get(file,char);
      else exit;
      end if;
    end loop;
    if min
     then i := -res;
     else i := res;
    end if;
    ni := cnt;
  end Parse;

  procedure Parse_also_Brackets
               ( file : in file_type;
                 char : in out character; i : out integer32;
                 ni : out natural32; sign : out character ) is

    bracket : boolean := false;

  begin
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    if char = '(' then
      bracket := true;
      get(file,char);
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    end if;
    Standard_Parse_Numbers.Parse(file,char,i,ni,sign);
    if bracket then
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
      if char = ')' then
        get(file,char);
      end if;
    end if;
  end Parse_also_Brackets;

  procedure Parse ( s : in string; p : in out integer;
                    i : out integer32; ni : out natural32;
                    sign : out character ) is

    res : integer32 := 0;
    min : boolean := false;
    cnt,temp : natural32 := 0;

  begin
    i := 0; ni := 0; sign := '+';
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    if (p > s'last) then return; end if;
    if (s(p) = '+') or (s(p) = '-') then
      min := (s(p) = '-');
      sign := s(p);
      p := p + 1;
    end if;
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    loop
      temp := Convert(s(p));
      if temp < 10 then
        if cnt < 9 then
          res := res*10 + integer32(temp);
          cnt := cnt + 1;
        else
          null;
        end if;
        p := p + 1;
      else
        exit;
      end if;
      exit when (p > s'last);
    end loop;
    if min
     then i := -res;
     else i := res;
    end if;
    ni := cnt;
  end Parse;

  procedure Parse_also_Brackets
               ( s : in string; p : in out integer;
                 i : out integer32; ni : out natural32;
                 sign : out character ) is

    bracket : boolean := false;

  begin
    i := 0; ni := 0; sign := '+';
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    if (p > s'last) then return; end if;
    if (s(p) = '(') then
      bracket := true;
      p := p + 1;
      if (p > s'last) then return; end if;
      Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
      if (p > s'last) then return; end if;
    end if;
    Standard_Parse_Numbers.Parse(s,p,i,ni,sign);
    if bracket then
      Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
      if s(p) = ')' then
        p := p + 1;
      end if;
    end if;
  end Parse_also_Brackets;

  procedure Parse ( file : in file_type;
                    char : in out character; i1,i2 : out integer32;
                    ni1,ni2 : out natural32; sign : out character ) is

    res1,res2 : integer32 := 0;
    min : boolean := false;
    k1,k2,temp : natural32 := 0;

  begin
    i1 := 0; i2 := 0; ni1 := 0; ni2 := 0; sign := '+';
    Skip_Spaces_and_CR(file,char);
    if (char = '+') or (char = '-') then
      min := (char = '-');
      sign := char;
      if end_of_file(file) then return; end if;
      get(file,char); 
    end if;
    Skip_Spaces_and_CR(file,char);
    loop
      temp := Convert(char);
      if temp < 10 then
        if k1 < 9 then
          res1 := res1*10 + integer32(temp);
          k1 := k1 + 1;
        elsif k2 < 9 then
          res2 := res2*10 + integer32(temp);
          k2 := k2 + 1;
        else
          null;  -- skip the rest of the numbers
        end if;      
        exit when end_of_file(file);
        get(file,char);
      else
        exit;
      end if;
    end loop;
    if min
     then i1 := -res1; i2 := -res2;
     else i1 := res1;  i2 := res2;
    end if;
    ni1 := k1;
    ni2 := k2;
  end Parse;

  procedure Parse_also_Brackets
               ( file : in file_type;
                 char : in out character; i1,i2 : out integer32;
                 ni1,ni2 : out natural32; sign : out character ) is

    bracket : boolean := false;

  begin
    Skip_Spaces_and_CR(file,char);
    if (char = '(') then
      bracket := true;
      get(file,char);
      Skip_Spaces_and_CR(file,char);
    end if;
    Standard_Parse_Numbers.Parse(file,char,i1,i2,ni1,ni2,sign);
    if bracket then
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
      if char = ')' then
        get(file,char);
      end if;
    end if;
  end Parse_also_Brackets;

  procedure Parse ( s : in string; p : in out integer;
                    i1,i2 : out integer32; ni1,ni2 : out natural32;
                    sign : out character ) is

    res1,res2 : integer32 := 0;
    min : boolean := false;
    k1,k2,temp : natural32 := 0;

  begin
    i1 := 0; i2 := 0; ni1 := 0; ni2 := 0; sign := '+';
    Skip_Spaces_and_CR(s,p);
    if (p > s'last) then return; end if;
    if (s(p) = '+') or (s(p) = '-') then
      min := (s(p) = '-');
      sign := s(p);
      p := p + 1;
    end if;
    Skip_Spaces_and_CR(s,p);
    if (p > s'last) then return; end if;
    loop
      temp := Convert(s(p));
      if temp < 10 then
        if k1 < 9 then
          res1 := res1*10 + integer32(temp);
          k1 := k1 + 1;
        elsif k2 < 9 then
          res2 := res2*10 + integer32(temp);
          k2 := k2 + 1;
        else
          null;  -- skip the rest of the numbers
        end if;
        p := p + 1;
      else
        exit;
      end if;
      exit when (p > s'last);
    end loop;
    if min
     then i1 := -res1; i2 := -res2;
     else i1 := res1;  i2 := res2;
    end if;
    ni1 := k1;
    ni2 := k2;
  end Parse;

  procedure Parse_also_Brackets
              ( s : in string; p : in out integer;
                i1,i2 : out integer32; ni1,ni2 : out natural32;
                sign : out character ) is

    bracket : boolean := false;

  begin
    i1 := 0; i2 := 0; ni1 := 0; ni2 := 0; sign := '+';
    Skip_Spaces_and_CR(s,p);
    if (p > s'last) then return; end if;
    if (s(p) = '(') then
      bracket := true;
      p := p + 1;
      if (p > s'last) then return; end if;
      Skip_Spaces_and_CR(s,p);
      if (p > s'last) then return; end if;
    end if;
    Standard_Parse_Numbers.Parse(s,p,i1,i2,ni1,ni2,sign);
    if bracket then
      Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
      if s(p) = ')' then
        p := p + 1;
      end if;
    end if;
  end Parse_also_Brackets;

  procedure Parse ( file : in file_type;
                    char : in out character; f : out double_float ) is

    f_int1,f_int2,f_quot1,f_quot2,expo,expo2 : integer32 := 0;
    f_int,f_quot : double_float := 0.0;
    k1,k2,nq1,nq2,np1,np2,temp : natural32 := 0;
    sign : character;
    min : boolean;

  begin
    Parse(file,char,f_int1,f_int2,np1,np2,sign);
    f_int := double_float(f_int1) * 10.0**integer(np2) + double_float(f_int2);
    min := (sign = '-');
    case char is
      when '.' =>
        get(file,char);       -- skip the point
        temp := Convert(char);
        if temp < 10 then
          Parse(file,char,f_quot1,f_quot2,nq1,nq2,sign);
          f_quot := double_float(f_quot1) * 10.0**integer(nq2) 
                  + double_float(f_quot2);
        end if;
        if (char = 'E') or (char = 'e') then
          get(file,char); -- skip the 'E'/'e'
          Parse(file,char,expo,expo2,k1,k2,sign);
        end if;
      when 'E' => 
        get(file,char); -- skip the 'E'/'e'
        Parse(file,char,expo,expo2,k1,k2,sign);
      when 'e' => 
        get(file,char); -- skip the 'E'/'e'
        Parse(file,char,expo,expo2,k1,k2,sign);
      when others => null;
    end case; 
    if min then
      if (f_int = 0.0) and (f_quot = 0.0) and (nq1 = 0) and (np1 = 0)
       then f := -1.0;   --  "-x" = -1*x 
       else f := ( f_int 
          - f_quot*10.0**(-integer(nq1)-integer(nq2)) )*10.0**integer(expo);
      end if;
    else
      f := ( f_int
           + f_quot*10.0**(-integer(nq1)-integer(nq2)) )*10.0**integer(expo);
    end if;
  end Parse;

  procedure Parse ( s : in string; p : in out integer;
                    f : out double_float ) is

    f_int1,f_int2,f_quot1,f_quot2,expo,expo2 : integer32 := 0;
    f_int,f_quot : double_float := 0.0;
    k1,k2,nq1,nq2,np1,np2,temp : natural32 := 0;
    sign : character;
    min : boolean;

  begin
    Parse(s,p,f_int1,f_int2,np1,np2,sign);
    f_int := double_float(f_int1) * 10.0**integer(np2) + double_float(f_int2);
    f := 0.0; min := (sign = '-');
    if p > s'last then return; end if;
    case s(p) is
      when '.' =>
        p := p + 1;                              -- skip the point
        temp := Convert(s(p));
        if temp < 10 then
          Parse(s,p,f_quot1,f_quot2,nq1,nq2,sign);
          f_quot := double_float(f_quot1) * 10.0**integer(nq2)
                  + double_float(f_quot2);
        end if;
        if (s(p) = 'E') or (s(p) = 'e') then
          p := p + 1;                            -- skip the 'E'/'e'
          Parse(s,p,expo,expo2,k1,k2,sign);
        end if;
      when 'E' | 'e' => 
        p := p + 1;                              -- skip the 'E'/'e'
        Parse(s,p,expo,expo2,k1,k2,sign);
      when others => null;
    end case; 
    if min then
      if (f_int = 0.0) and (f_quot = 0.0) and (nq1 = 0) and (np1 = 0)
       then f := -1.0;   --  "-x" = -1*x 
       else f := ( f_int
          - f_quot*10.0**(-integer(nq1)-integer(nq2)) )*10.0**integer(expo);
      end if;
    else
      f := ( f_int
           + f_quot*10.0**(-integer(nq1)-integer(nq2)) )*10.0**integer(expo);
    end if;
  end Parse;

  procedure Parse ( file : in file_type;
                    char : in out character; c : out Complex_Number ) is
 
    f1,f2 : double_float;

  begin
    Parse(file,char,f1);
    if char = '/' then
      get(file,char);           -- skip the '/'
      Parse(file,char,f2);
      c := Create(f1/f2);
    else 
      c := Create(f1);
    end if;
  exception
    when constraint_error => raise INFINITE_NUMBER;
  end Parse;

  procedure Parse ( s : in string; p : in out integer;
                    c : out Complex_Number ) is
 
    f1,f2 : double_float;

  begin
    Parse(s,p,f1);
    if s(p) = '/' then
      p := p + 1;               -- skip the '/'
      Parse(s,p,f2);
      c := Create(f1/f2);
    else
      c := Create(f1);
    end if;
  exception
    when constraint_error => raise INFINITE_NUMBER;
  end Parse;

end Standard_Parse_Numbers;
