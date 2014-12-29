with integer_io;
with String_Parsing;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Characters_and_Numbers;
with Standard_Complex_Vector_Strings;

package body Standard_Solution_Strings is

-- PART I: compute the lengths and write the strings

  function Length_Number ( f : double_float ) return natural32 is

    s : string(1..22) := (1..22 => ' ');

  begin
    put(s,f);
    if (s(1) = ' ' and (s(2) = ' ' or s(2) = '-'))
     then return 21;
     else return 22;
    end if;
  end Length_Number;

  function Length_Number ( c : Complex_Number ) return natural32 is
  begin
    return Length_Number(REAL_PART(c)) + 2 + Length_Number(IMAG_PART(c));
  end Length_Number;

  function Write_Number ( c : Complex_Number ) return string is

    re : string(1..integer(Length_Number(REAL_PART(c))));
    im : string(1..integer(Length_Number(IMAG_PART(c))));

  begin
    put(re,REAL_PART(c));
    put(im,IMAG_PART(c));
    declare
      result : constant string := re & "  " & im;
    begin
      return result;
    end;
  end Write_Number;

  function Number_of_Digits ( i : natural32 ) return natural32 is
  begin
    if i < 10
     then return 1;
     else return 1 + Number_of_Digits(i/10);
    end if;
  end Number_of_Digits;

  function Length_Intro ( s : Solution ) return natural32 is

    st : constant natural32 := 4 + Length_Number(s.t) + 1;
    sm : constant natural32 := 4 + Number_of_Digits(natural32(s.m)) + 1;
    bt : constant natural32 := 21;

  begin
    return st + sm + bt;
  end Length_Intro;

  function Write_Intro ( t : Complex_Number; m : integer32 ) return string is

    st : constant string := "t : " & Write_Number(t) & ASCII.LF;
    sm : constant string := "m : " & Characters_and_Numbers.Convert(m)
                                   & ASCII.LF;
    bt : constant string := "the solution for t :" & ASCII.LF;

  begin
    return st & sm & bt; 
  end Write_Intro;

  function Write_Intro ( s : Solution ) return string is
  begin
    return Write_Intro(s.t,s.m);
  end Write_Intro;

  function Length_Symbol ( i : natural32 ) return natural32 is
  begin
    if Symbol_Table.Number < i then
      return 1 + Number_of_Digits(i);
    else
      declare
        sb : constant Symbol_Table.Symbol := Symbol_Table.Get(i);
      begin
        for k in sb'range loop
          if sb(k) = ' '
           then return natural32(k-1);
          end if;
        end loop;
        return natural32(sb'last);
      end;
    end if;
  end Length_Symbol;

  function Write_Symbol ( i : natural32 ) return string is
  begin
    if Symbol_Table.Number < i then
      declare
        res : constant string := "x" & Characters_and_Numbers.nConvert(i);
      begin
        return res;
      end;
    else
      declare
        sb : constant Symbol_Table.Symbol := Symbol_Table.Get(i);
      begin
        for k in sb'range loop
          if sb(k) = ' '
           then return sb(sb'first..k-1);
          end if;
        end loop;
        return sb;
      end;
    end if;
  end Write_Symbol;

  function Length_Components
              ( k : natural32; v : Vector; accu : natural32 )
              return natural32 is

  -- DESCRIPTION :
  --   If k equals zero, then the accumulator is returned,
  --   otherwise, the length of the k-th component is added to the accu.

  begin
    if k = 0 then
      return accu;
    else
      declare
        sb : constant natural32 := 1 + Length_Symbol(k) + 3;
        nb : constant natural32 := Length_Number(v(integer32(k))) + 1;
        new_accu : constant natural32 := sb + nb;
      begin
        return Length_Components(k-1,v,new_accu + accu);
      end;
    end if;
  end Length_Components;

  function Write_Components
              ( k : natural32; v : Vector; accu : string )
              return string is

  -- DESCRIPTION :
  --   If k equals zero, then the accumulator is returned,
  --   otherwise, the k-th component is added to the accu.

  begin
    if k = 0 then
      return accu;
    else
      declare
        sb : constant string := " " & Write_Symbol(k) & " : ";
        nb : constant string := Write_Number(v(integer32(k))) & ASCII.LF;
        new_accu : constant string := sb & nb;
      begin
        return Write_Components(k-1,v,new_accu & accu);
      end;
    end if;
  end Write_Components;

  function Add_Linefeed ( s : string ) return string is

  -- DESCRIPTION :
  --   If s(s'last) /= ASCII.LF, then it will be added to the string s,
  --   otherwise the string s is returned.

  begin
    if s(s'last) = ASCII.LF
     then return s;
     else return s & ASCII.LF;
    end if;
  end Add_Linefeed;

  function Write_Components
             ( k,n : integer32; xv,accu : string ) return string is

  -- DESCRIPTION :
  --   If k equals n, or is larger than n, then the accu is returned,
  --   otherwise, the k-th symbol is added to the string, along with
  --   the number on the current line of xv. 

  begin
    if k > n then  -- be careful if n equals 1 !!!
      return accu;
    else
      declare
        sb : constant string := " " & Write_Symbol(natural32(k)) & " : ";
        pos : constant integer
            := Standard_Complex_Vector_Strings.Next_Linefeed(xv);
        nb : constant string := Add_Linefeed(xv(xv'first..pos));
        new_accu : constant string := sb & nb;
      begin
        if k = n then
          return accu & new_accu;
        else -- skip the linefeed in xv
          return Write_Components(k+1,n,xv(pos+1..xv'last),accu & new_accu);
        end if;
      end;
    end if;
  end Write_Components;

  function Length_Vector ( v : Vector ) return natural32 is
  begin
    return Length_Components(natural32(v'last),v,0);
  end Length_Vector;

  function Length_Vector ( s : Solution ) return natural32 is
  begin
    return Length_Components(natural32(s.n),s.v,0);
  end Length_Vector;

  function Write_Vector ( v : Vector ) return string is
  begin
    return Write_Components(natural32(v'last),v,"");
  end Write_Vector;

  function Write_Vector ( s : Solution ) return string is
  begin
    return Write_Components(natural32(s.n),s.v,"");
  end Write_Vector;

  function Length_Diagnostics return natural32 is
  begin
    return 59;
  end Length_Diagnostics;

  function Write_Diagnostics ( err,rco,res : double_float ) return string is

    s_err,s_rco,s_res : string(1..10);

  begin
    put(s_err,err,3);
    put(s_rco,rco,3);
    put(s_res,res,3);
    declare
      result : constant string
             := "== err : " & s_err & " = rco : " & s_rco 
              & " = res : " & s_res & " =";
    begin
      return result;
    end;
  end Write_Diagnostics;

  function Write_Diagnostics ( s : Solution ) return string is
  begin
    return Write_Diagnostics(s.err,s.rco,s.res);
  end Write_Diagnostics;

  function Length ( s : Solution ) return natural32 is
  begin
    return Length_Intro(s) + Length_Vector(s) + Length_Diagnostics;
  end Length;

  function Write ( s : Solution ) return string is

    tm : constant string := Write_Intro(s);
    sv : constant string := Write_Vector(s.v);
    dg : constant string := Write_Diagnostics(s.err,s.rco,s.res);

  begin
    return tm & sv & dg;
  end Write;

  function Write ( t : Complex_Number; n,m : integer32; xv : string; 
                   err,rco,res : double_float ) return string is

    result : constant string
           := Write_Intro(t,m)
            & Write_Components(1,n,xv,"")
            & Write_Diagnostics(err,rco,res);

  begin
    return result;
  end Write;

-- PART II: parse strings into solutions

  procedure Parse_Intro
               ( s : in string; k : in out integer; 
                 t : out Complex_Number; m : out integer32;
                 fail : out boolean ) is

    ind : integer := String_Parsing.Scan(s(k..s'last),":");
    last : integer;

  begin
    t := Create(0.0); m := 0; fail := false;
    if ind > 0 then
      declare
        re,im : double_float := 0.0;
      begin
        get(s(ind+1..s'last),re,last);
        get(s(last+1..s'last),im,last);
        t := Create(re,im);
      exception
        when others => fail := true; return;
      end;
      ind := String_Parsing.Scan(s(last+1..s'last),":");
      if ind > 0 then
        declare
        begin
          integer_io.get(s(ind+1..s'last),integer(m),last);
          k := last;
        exception
          when others => fail := true; return;
        end;
      end if;
    end if;
  exception
    when others => fail := true; return;
  end Parse_Intro;

  procedure Parse_Symbol ( s : in String; k : in out integer;
                           sb : out Symbol; fail : out boolean ) is

    ind : integer := sb'first-1;

  begin
    fail := false;
    for i in sb'range loop
      sb(i) := ' ';
    end loop;
    while s(k) = ' ' or s(k) = ASCII.LF or s(k) = ASCII.CR loop
      k := k+1;
      exit when k > s'last;
    end loop;
    while s(k) /= ' ' loop
      ind := ind + 1;
      sb(ind) := s(k);
      k := k + 1;
      exit when k > s'last;
    end loop;
    if k <= s'last then
      while s(k) /= ':' loop
        k := k + 1;
      end loop;
    end if;
  exception
    when others => fail := true; return;
  end Parse_Symbol;

  procedure Parse_Vector
               ( s : in string; k : in out integer; n : in natural32;
                 v : out Standard_Complex_Vectors.Vector;
                 fail : out boolean ) is

    re,im : double_float := 0.0;
    ind : integer := k;
    sb : Symbol;
    sbk : integer32;

  begin
    fail := false;
    for i in 1..integer32(n) loop
      v(i) := Create(0.0);
    end loop;
    for i in 1..integer32(n) loop
      Parse_Symbol(s,ind,sb,fail);
      exit when fail;
      get(s(ind+1..s'last),re,ind);
      get(s(ind+1..s'last),im,ind);
      if Symbol_Table.Number < n then
        v(i) := Create(re,im);
        if Symbol_Table.Empty
         then Symbol_Table.Init(n);
        end if;
        Symbol_Table.Add(sb);
      else
        sbk := integer32(Symbol_Table.get(sb));
        if sbk > 0 and sbk <= integer32(n) then
          v(sbk) := Create(re,im);
        else
          v(i) := Create(re,im); 
        end if;
      end if;
      while s(ind) /= ASCII.LF and ind < s'last loop
        ind := ind + 1;
      end loop;
    end loop;
    k := ind;
  exception
    when others => fail := true; return;
  end Parse_Vector;

  procedure Parse_Diagnostics
               ( s : in string; k : in out integer;
                 err,rco,res : out double_float; fail : out boolean ) is

    ind : integer := String_Parsing.Scan(s(k..s'last),"err :");
    last : integer;

  begin
    err := 0.0; rco := 0.0; res := 0.0; fail := false;
    if ind > 0 then
      declare
      begin
        get(s(ind+1..s'last),err,last);
      exception
        when others => fail := true; return;
      end;
      ind := String_Parsing.Scan(s(last..s'last),"rco :");
      if ind < 0 then
        rco := 0.0; res := 0.0;       
      else
        declare
        begin
          get(s(ind+1..s'last),rco,last);
        exception
          when others => fail := true; return;
        end;
        ind := String_Parsing.Scan(s(last..s'last),"res :");
        if ind < 0 then
          res := 0.0;
        else
          declare
          begin
            get(s(ind+1..s'last),res,last);
          exception
            when others => fail := true; return;
          end;
        end if;
      end if;
    end if;
    k := last;
  exception
    when others => fail := true; return;
  end Parse_Diagnostics;

  procedure Parse ( s : in string; k : in out integer; n : in natural32;
                    sol : out Solution; fail : out boolean ) is

    ind : integer := s'first;
  
  begin
    Parse_Intro(s,ind,sol.t,sol.m,fail);
    ind := String_Parsing.Scan(s(ind+1..s'last),":");
    if not fail then
      ind := ind + 1;  -- scan marks the position of the ":"
      Parse_Vector(s,ind,n,sol.v,fail);
      if not fail then
        Parse_Diagnostics(s,ind,sol.err,sol.rco,sol.res,fail);
      end if;
    end if;
  exception 
    when others => fail := true; return;
  end Parse;

end Standard_Solution_Strings;
