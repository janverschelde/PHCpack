--with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Parse_Numbers;

package body Multprec_Parse_Numbers is

  procedure Parse ( file : in file_type; char : in out character;
                    i : out Integer_Number;
                    n : out natural32; sign : out character ) is

    res : Integer_Number := Create(integer32(0));
    min : boolean := false;
    cnt : natural32 := 0;
    temp : integer32 := 0;

  begin
    n := 0; sign := '+';
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    if (char = '+') or (char = '-') then
      min := (char = '-');
      sign := char;
      if end_of_file(file) then return; end if;
      get(file,char); 
    end if;
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    loop
      temp := integer32(Convert(char));
      if temp < 10 then
        Mul(res,10);
        Add(res,temp);
        cnt := cnt + 1;
        exit when end_of_file(file);
        get(file,char);
      else
        exit;
      end if;
    end loop;
    if min and not Equal(res,0)
     then Multprec_Integer_Numbers.Min(res);
    end if;
    i := res;
    n := cnt;
  end Parse;

  procedure Parse_also_Brackets
                  ( file : in file_type; char : in out character;
                    i : out Integer_Number;
                    n : out natural32; sign : out character ) is

    bracket : boolean := false;

  begin
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    if char = '(' then
      bracket := true;
      get(file,char);
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    end if;
    Multprec_Parse_Numbers.Parse(file,char,i,n,sign);
    if bracket then
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
      if char = ')' then
        get(file,char);
      end if;
    end if;
  end Parse_also_Brackets;

  procedure Parse ( s : in string; p : in out integer;
                    i : out Integer_Number;
                    n : out natural32; sign : out character ) is

    res : Integer_Number := Create(integer32(0));
    min : boolean := false;
    cnt : natural32 := 0;
    temp : integer32 := 0;

  begin
    n := 0; sign := '+';
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    if (p > s'last) then return; end if;
    if (s(p) = '+') or (s(p) = '-') then
      min := (s(p) = '-');
      sign := s(p);
      p := p + 1;
    end if;
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    if (p > s'last) then return; end if;
    loop
      temp := integer32(Convert(s(p)));
      if temp < 10 then
        Mul(res,10);
        Add(res,temp);
        cnt := cnt + 1;
        p := p + 1;
      else
        exit;
      end if;
      exit when (p > s'last);
    end loop;
    if min and not Equal(res,0)
     then Multprec_Integer_Numbers.Min(res);
    end if;
    i := res;
    n := cnt;
  end Parse;

  procedure Parse_also_Brackets
                  ( s : in string; p : in out integer;
                    i : out Integer_Number;
                    n : out natural32; sign : out character ) is

    bracket : boolean := false;

  begin
    n := 0; sign := '+';
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    if (p > s'last) then return; end if;
    if (s(p) = '(') then
      bracket := true;
      p := p + 1;
      if (p > s'last) then return; end if;
      Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
      if (p > s'last) then return; end if;
    end if;
    Multprec_Parse_Numbers.Parse(s,p,i,n,sign);
    if bracket then
      Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
      if s(p) = ')' then
        p := p + 1;
      end if;
    end if;
  end Parse_also_Brackets;

  procedure Parse ( file : in file_type; char : in out character;
                    f : out Floating_Number ) is

    int,quot,expo : Integer_Number;
    f_int,f_quot,acc : Floating_Number := Create(0.0);
    nq,np,ne,temp : natural32 := 0;
    sign : character;
    min : boolean;

  begin
    Parse(file,char,int,np,sign);
    f_int := Create(int);
    min := (sign = '-');
    case char is
      when '.'  =>
        get(file,char);                       -- skip the point
        temp := Convert(char);
        if temp < 10 then
          Parse(file,char,quot,nq,sign);
          f_quot := Create(quot);
        end if;
        if (char = 'E') or (char = 'e') then
          get(file,char);                     -- skip the 'E' or 'e'
          Parse(file,char,expo,ne,sign);
        end if;
      when 'E' | 'e' => 
        get(file,char);                       -- skip the 'E' or 'e'
        Parse(file,char,expo,ne,sign);
      when others => null;
    end case; 
    if min then
      if Equal(f_int,0.0) and Equal(f_quot,0.0) 
        and (nq = 0) and (np = 0) then
        f := Create(integer32(-1));   --  "-x" = -1*x 
      else --f := ( f_int - f_quot*10.0**(-nq) )*10.0**expo ;
        acc := f_quot*10.0**(-integer(nq));
        Multprec_Floating_Numbers.Set_Size
          (f_int,Size_Fraction(f_int)+Decimal_to_Size(nq)+1);
        Sub(f_int,acc); Clear(acc);
        acc := 10.0**expo;
        Mul(f_int,acc); Clear(acc);
        f := f_int;
      end if;
    else --f := ( f_int + f_quot*10.0**(-nq) )*10.0**expo ;
      acc := f_quot*10.0**(-integer(nq));
      Multprec_Floating_Numbers.Set_Size
        (f_int,Size_Fraction(f_int)+Decimal_to_Size(nq)+1);
      Add(f_int,acc); Clear(acc);
      acc := 10.0**expo;
      Mul(f_int,acc); Clear(acc);
      f := f_int;
    end if;
  end Parse;

  procedure Parse ( s : in string; p : in out integer;
                    f : out Floating_Number ) is

    int,quot,expo : Integer_Number;
    f_int,f_quot,acc : Floating_Number := Create(0.0);
    nq,np,ne,temp : natural32 := 0;
    sign : character;
    min : boolean;

  begin
   -- put("starting to parse integer ...");
    Parse(s,p,int,np,sign);
   -- put("done parsing integer, int = "); put(int); new_line;
   -- put("p = "); put(p,1); new_line;
    f_int := Create(int);
    min := (sign = '-');
    case s(p) is
      when '.' => 
        p := p + 1;                          -- skip the point
        temp := Convert(s(p));
        if temp < 10 then
          Parse(s,p,quot,nq,sign);
          f_quot := Create(quot);
        end if;
        if (s(p) = 'E') or (s(p) = 'e') then
          p := p + 1;                        -- skip the 'E'/'e'
          Parse(s,p,expo,ne,sign);
        end if;
      when 'E' | 'e' => 
        p := p + 1;                          -- skip the 'E'/'e'
        Parse(s,p,expo,ne,sign);
      when others => null;
    end case; 
   -- put_line("past the case ...");
    if min then
     -- put_line("sign is negative");
      if Equal(f_int,0.0) and Equal(f_quot,0.0)
         and (nq = 0) and (np = 0) then
        f := Create(integer32(-1));   --  "-x" = -1*x 
      else --f := ( f_int - f_quot*10.0**(-nq) )*10.0**expo ;
        acc := f_quot*10.0**(-integer(nq));
        Multprec_Floating_Numbers.Set_Size
          (f_int,Size_Fraction(f_int)+Decimal_to_Size(nq)+1);
        Sub(f_int,acc); Clear(acc);
        acc := 10.0**expo;
        Mul(f_int,acc); Clear(acc);
        f := f_int;
      end if;
    else --f := ( f_int + f_quot*10.0**(-nq) )*10.0**expo ;
     -- put_line("sign is positive");
     -- put("nq = "); put(nq,1); new_line;
      if nq = 0
       then acc := f_quot;
       else acc := f_quot*10.0**(-integer(nq));
      end if;
     -- put("Size_Fraction(f_int) = ");
     -- put(Size_Fraction(f_int),1); new_line;
     -- put("Decimal_to_Size(nq) = ");
     -- put(Decimal_to_Size(nq),1); new_line;
     -- put("Size_Fraction(f_int)+Decimal_to_Size(nq)+1 = ");
     -- put(Size_Fraction(f_int)+Decimal_to_Size(nq)+1,1); new_line;
     -- put("calling set_size...");
      Multprec_Floating_Numbers.Set_Size
        (f_int,Size_Fraction(f_int)+Decimal_to_Size(nq)+1);
     -- put_line("after set_size...");
      Add(f_int,acc); Clear(acc);
      acc := 10.0**expo;
      Mul(f_int,acc); Clear(acc);
      f := f_int;
    end if;
--  exception 
--    when others => put_line("exception occurred when parsing float ...");
--                   raise;
  end Parse;

  procedure Parse ( file : in file_type; size : in natural32;
                    char : in out character; c : out Complex_Number ) is

    f1,f2 : Floating_Number;

  begin
    Parse(file,char,f1);
    Set_Size(f1,size);
    if char = '/' then
      get(file,char);           -- skip the '/'
      Parse(file,char,f2);
      Set_Size(f2,size);
      c := Create(f1/f2);
    else
      Set_Size(f1,size);
      c := Create(f1);
    end if;
  exception
    when constraint_error => raise INFINITE_NUMBER;
  end Parse;

  procedure Parse ( s : in string; size : in natural32; p : in out integer;
                    c : out Complex_Number ) is

    f1,f2 : Floating_Number;

  begin
   -- put("in parse with p = "); put(p,1); new_line;
    Parse(s,p,f1);
    Set_Size(f1,size);
   -- put("s'last = "); put(s'last,1); new_line;
   -- put("p = "); put(p,1); new_line;
    if s(p) = '/' then
      p := p + 1;               -- skip the '/'
      Set_Size(f2,size);
      Parse(s,p,f2);
      c := Create(f1/f2);
    else 
      Set_Size(f1,size);
      c := Create(f1);
    end if;
  exception
    when constraint_error => raise INFINITE_NUMBER;
  end Parse;

end Multprec_Parse_Numbers;
