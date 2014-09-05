with Characters_and_Numbers;             use Characters_and_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;

package body Multprec_Integer_Numbers_io is

  procedure get ( lc : in out character; i : in out Integer_Number ) is
  begin
    get(Standard_Input,lc,i);
  end get;

  procedure get ( file : in file_type;
                  lc : in out character;  i : in out Integer_Number ) is

    n : Natural_Number;
    plus : boolean;

  begin
    if lc = ' '
     then Skip_Spaces(file,lc);
    end if;
    if lc = '+' then
      plus := true;
      get(file,lc);
    elsif lc = '-' then
      plus := false;
      get(file,lc);
    else
      plus := true;
    end if;
    get(file,lc,n);
    Clear(i);
    i := Create(n);
    if not plus
     then Min(i);
    end if;
  end get;

  procedure get ( i : in out Integer_Number ) is
  begin
    get(Standard_Input,i);
  end get;

  procedure get ( file : in file_type; i : in out Integer_Number ) is

    c : character := ' ';

  begin
    Skip_Spaces(file,c);
    get(file,c,i);
  end get;

  procedure get ( lc : in out character; i : in out Integer_Number;
                  nonnegative : out boolean ) is
  begin
    get(Standard_Input,lc,i,nonnegative);
  end get;

  procedure get ( file : in file_type;
                  lc : in out character;  i : in out Integer_Number;
                  nonnegative : out boolean ) is

    n : Natural_Number;
    plus : boolean;

  begin
    if lc = ' '
     then Skip_Spaces(file,lc);
    end if;
    if lc = '+' then
      plus := true;
      get(file,lc);
    elsif lc = '-' then
      plus := false;
      get(file,lc);
    else
      plus := true;
    end if;
    get(file,lc,n);
    Clear(i);
    i := Create(n);
    if not plus
     then Set_Min(i);
    end if;
    nonnegative := not plus;
  end get;

  procedure get ( s : in string; i : in out Integer_Number ) is

    minus : boolean;
    n : Natural_Number;
    pos : integer := s'first;

  begin
    while s(pos) = ' ' loop
      pos := pos + 1;
      exit when pos > s'last;
    end loop;
    if pos <= s'last then
      if s(pos) = '-' then
        minus := true;
        get(s(pos+1..s'last),n);
      else
        minus := false;
        if s(pos) = '+'
         then get(s(pos+1..s'last),n);
         else get(s(pos..s'last),n);
        end if;
      end if;
      i := Create(n);
      if minus
       then Min(i);
      end if;
    end if;
  end get;

  procedure put ( i : in Integer_Number ) is
  begin
    put(Standard_Output,i);
  end put;

  procedure put ( file : in file_type; i : in Integer_Number ) is
  begin
    if Empty(i) then
      put(file,"0");
    else
      if Negative(i)
       then put(file,"-");
      end if;
      put(file,Unsigned(i));
    end if;
  end put;

  procedure put ( s : out string; i : in Integer_Number ) is
  begin
    if Empty(i) then
      s(s'first) := '0';
    else
      if Negative(i) then
        s(s'first) := '-';
        put(s(s'first+1..s'last),Unsigned(i));
      else
        put(s,Unsigned(i));
      end if;
    end if;
  end put;

  function One_if_Negative ( i : Integer_Number ) return natural32 is

  -- DESCRIPTION :
  --   Returns one if i < 0, returns 0 otherwise.

  begin
    if i < 0
     then return 1;
     else return 0;
    end if;
  end One_if_Negative;

  function Convert_to_String ( i : Integer_Number ) return string is

    dp : constant natural32 := Multprec_Integer_Numbers.Decimal_Places(i);
    dp1 : constant natural32 := dp + One_if_Negative(i);
    res : string(1..integer(dp1));

  begin
    if dp = 0 then
      return "0";
    else
      put(res,i);
      return res;
    end if;
  end Convert_to_String;

  procedure put ( i : in Integer_Number; dp : in natural32 ) is
  begin
    put(Standard_Output,i,dp);
  end put;
  
  procedure put ( file : in file_type;
                  i : in Integer_Number; dp : in natural32 ) is
  begin
    if Empty(i) then
      for k in 1..dp-1 loop
        put(file," ");
      end loop;
      put(file,"0");
    else
      if Negative(i) then
        for k in 1..dp-Decimal_Places(i)-1 loop
          put(file," ");
        end loop;
        put(file,"-");
        put(file,Unsigned(i));
      else
        put(file,Unsigned(i),dp);
      end if;
     -- put(file,Unsigned(i));
    end if;  
  end put;

end Multprec_Integer_Numbers_io;
