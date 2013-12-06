with Characters_and_Numbers;             use Characters_and_Numbers;
with Multprec_Natural64_Numbers;         use Multprec_Natural64_Numbers;
with Multprec_Natural64_Numbers_io;      use Multprec_Natural64_Numbers_io;

package body Multprec_Integer64_Numbers_io is

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

  begin
    if s(s'first) = '-' then
      minus := true;
      get(s(s'first+1..s'last),n);
    else
      minus := false;
      if s(s'first) = '+'
       then get(s(s'first+1..s'last),n);
       else get(s,n);
      end if;
    end if;
    i := Create(n);
    if minus
     then Min(i);
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

end Multprec_Integer64_Numbers_io;
