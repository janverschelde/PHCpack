with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;

package body Standard_Complex_Numbers_io is

  procedure get ( c : in out Complex_Number ) is
  begin
    get(Standard_Input,c);
  end get;

  procedure get ( file : in file_type; c : in out Complex_Number ) is

    x,y : double_float := 0.0;

  begin
    get(file,x); get(file,y);
    c := Create(x,y);
  end get;

  procedure get ( s : in string; c : in out Complex_Number;
                  last : out integer ) is

    x,y : double_float := 0.0;

  begin
    get(s,x,last);
    while s(last) /= ' ' loop  -- make sure to start at dividing spaces
      last := last + 1;
      exit when last >= s'last;
    end loop;
    get(s(last..s'last),y,last);
    c := Create(x,y);
  end get;

  procedure put ( c : in Complex_Number ) is
  begin
    put(Standard_Output,c);
  end put;

  procedure put ( file : in file_type; c : in Complex_Number ) is
  begin
    put(file,REAL_PART(c));
    put(file,"  ");
    put(file,IMAG_PART(c));
  end put;

  procedure put ( s : out string; c : in Complex_Number ) is

    sre,sim : string(1..22) := (1..22 => ' '); -- careful with E-100

  begin
    put(sre,REAL_PART(c));
    put(sim,IMAG_PART(c));
    s := sre & "  " & sim;
  end put;

  procedure put ( c : in Complex_Number; fore,aft,exp : in natural32 ) is
  begin
    put(Standard_Output,c,fore,aft,exp);
  end put;

  procedure put ( file : in file_type;
                  c : in Complex_Number; fore,aft,exp : in natural32 ) is
  begin
    put(file,REAL_PART(c),fore,aft,exp);
    put(file,"  ");
    put(file,IMAG_PART(c),fore,aft,exp);
  end put;

  function Trim_Trailing_Spaces ( s : string ) return string is

  -- DESCRIPTION :
  --   Returns the string s with trailing spaces removed.

    last : integer := s'last;

  begin
    while s(last) = ' ' loop
      last := last - 1;
      exit when last <= s'first;
    end loop;
    return s(s'first..last);
  end Trim_Trailing_Spaces;

  function Trim_Leading_Spaces ( s : string ) return string is

  -- DESCRIPTION :
  --   Returns the string s with leading spaces removed.

    first : integer := s'first;

  begin
    while s(first) = ' ' loop
      first := first + 1;
      exit when first >= s'last;
    end loop;
    if s(first) = '-'
     then return s(first..s'last);
     else return (' ' & s(first..s'last)); -- insert space for +
    end if;
  end Trim_Leading_Spaces;

  function Trim_Spaces ( s : string ) return string is

  -- DESCRIPTION :
  --   Trims leading and trailing spaces of s.

  begin
    return Trim_Leading_Spaces(Trim_Trailing_Spaces(s));
  end Trim_Spaces;

  procedure put ( s : out string;
                  c : in Complex_Number; aft,exp : in natural32 ) is

    sre,sim : string(1..22) := (1..22 => ' ');

  begin
    put(sre,REAL_PART(c),aft,exp);
    put(sim,IMAG_PART(c),aft,exp);
    s := Trim_Spaces(sre) & "  " & Trim_Spaces(sim);
  end put;

  procedure put ( c : in Complex_Number; dp : in natural32 ) is
  begin
    put(c,dp,dp,dp);
  end put;

  procedure put ( file : in file_type;
                  c : in Complex_Number; dp : in natural32 ) is
  begin
    put(file,c,dp,dp,dp);
  end put;

  procedure put ( s : out string;
                  c : in Complex_Number; dp : in natural32 ) is

    sre,sim : string(1..22) := (1..22 => ' ');

  begin
    s := (s'range => ' ');
    put(sre,REAL_PART(c),dp);
    put(sim,IMAG_PART(c),dp);
    declare
      t : constant string
        := Trim_Spaces(sre) & "  " & Trim_Spaces(sim);
    begin
      s(t'range) := t;
    end;
  end put;

end Standard_Complex_Numbers_io;
