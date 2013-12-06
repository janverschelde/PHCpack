with integer_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with String_Parsing;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Parse_Numbers;
with Standard_Solution_Strings;

package body Multprec_Solution_Strings is

  procedure get ( s : in string; f : out Floating_Number;
                  last : out integer ) is

  -- DESCRIPTION :
  --   This wrapper for parsing multiprecision numbers achieves
  --   the same interface as in the body of Standard_Solution_Strings. 

   -- df : double_float := 0.0;
    p : integer;

  begin
   -- get(s,df,p);
   -- last := p;
   -- f := Create(df);
    p := s'first;
    Multprec_Parse_Numbers.Parse(s,p,f);
    last := p;
  end get;

-- PART I: compute the lengths and write the strings

  function Length_Number ( f : Floating_Number ) return natural32 is
  begin
    return Character_Size(f);
  end Length_Number;

  function Length_Number ( c : Complex_Number ) return natural32 is

    re : Floating_Number := REAL_PART(c);
    im : Floating_Number := IMAG_PART(c);
    res : constant natural32 := Length_Number(re) + 2 + Length_Number(im);

  begin
    Clear(re); Clear(im);
    return res;
  end Length_Number;

  function Write_Number ( c : Complex_Number ) return string is

  -- NOTE : the initialization of the strings with blank space is
  --   required for numbers that occupy less than the predicted amount.

    re_c : Floating_Number := REAL_PART(c);
    im_c : Floating_Number := IMAG_PART(c);
    len_re_c : constant integer := integer(Length_Number(re_c));
    len_im_c : constant integer := integer(Length_Number(im_c));
    re : string(1..len_re_c) := (1..len_re_c => ' ');
    im : string(1..len_im_c) := (1..len_im_c => ' ');

  begin
    put(re,re_c); Clear(re_c);
    put(im,im_c); Clear(im_c);
    declare
      result : constant string := re & "  " & im;
    begin
      return result;
    end;
  end Write_Number;

  function Number_of_Digits ( i : natural32 ) return natural32 is
  begin
    return Standard_Solution_Strings.Number_of_Digits(i);
  end Number_of_Digits;

  function Length_Intro ( s : Solution ) return natural32 is

    st : constant natural32 := 4 + Length_Number(s.t) + 1;
    sm : constant natural32 := 4 + Number_of_Digits(natural32(s.m)) + 1;
    bt : constant natural32 := 21;

  begin
    return st + sm + bt;
  end Length_Intro;

  function Write_Intro ( s : Solution ) return string is

    st : constant string := "t : " & Write_Number(s.t) & ASCII.LF;
    sm : constant string := "m : " & Convert(s.m) & ASCII.LF;
    bt : constant string := "the solution for t :" & ASCII.LF;

  begin
    return st & sm & bt; 
  end Write_Intro;

  function Length_Symbol ( i : natural32 ) return natural32 is
  begin
    return Standard_Solution_Strings.Length_Symbol(i);
  end Length_Symbol;

  function Write_Symbol ( i : natural32 ) return string is
  begin
    return Standard_Solution_Strings.Write_Symbol(i);
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

  function Write_Diagnostics ( err,rco,res : Floating_Number ) return string is

    s_err,s_rco,s_res : string(1..10);
    df_err : constant double_float := Round(err);
    df_rco : constant double_float := Round(rco);
    df_res : constant double_float := Round(res);

  begin
    put(s_err,df_err,3);
    put(s_rco,df_rco,3);
    put(s_res,df_res,3);
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

-- PART II: parse strings into solutions

  procedure Parse_Intro
               ( s : in string; k : in out integer;
                 t : out Complex_Number; m : out integer32;
                 fail : out boolean ) is

    ind : integer := String_Parsing.Scan(s(k..s'last),":");
    last : integer;

  begin
    t := Multprec_Complex_Numbers.Create(integer(0));
    m := 0; fail := false;
    if ind > 0 then
      declare
        re,im : Floating_Number;
      begin
        get(s(ind+1..s'last),re,last);
        get(s(last+1..s'last),im,last);
        t := Create(re,im);
        Clear(re); Clear(im);
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
  begin
    Standard_Solution_Strings.Parse_Symbol(s,k,sb,fail);
  end Parse_Symbol;

  procedure Parse_Vector
               ( s : in string; k : in out integer; n : in natural32;
                 v : out Multprec_Complex_Vectors.Vector;
                 fail : out boolean ) is

    re,im : Floating_Number;
    ind : integer := k;
    sb : Symbol;
    sbk : integer32;

  begin
    fail := false;
    for i in 1..integer32(n) loop
      v(i) := Multprec_Complex_Numbers.Create(integer(0));
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
                 err,rco,res : out Floating_Number; fail : out boolean ) is

    ind : integer := String_Parsing.Scan(s(k..s'last),"err :");
    last : integer;

  begin
    err := Create(0.0);
    rco := Create(0.0); 
    res := Create(0.0); fail := false;
    if ind > 0 then
      declare
      begin
        get(s(ind+1..s'last),err,last);
      exception
        when others => fail := true; return;
      end;
      ind := String_Parsing.Scan(s(last..s'last),"rco :");
      if ind < 0 then
        rco := Create(0.0); res := Create(0.0);       
      else
        declare
        begin
          get(s(ind+1..s'last),rco,last);
        exception
          when others => fail := true; return;
        end;
        ind := String_Parsing.Scan(s(last..s'last),"res :");
        if ind < 0 then
          res := Create(0.0);
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
    k := ind;
  exception 
    when others => fail := true; return;
  end Parse;

end Multprec_Solution_Strings;
