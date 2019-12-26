with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;

package body Multprec_Complex_Numbers_io is

  procedure get ( x : in out Complex_Number ) is
  begin
    get(Standard_Input,x);
  end get;

  procedure get ( file : in file_type; x : in out Complex_Number ) is

    xre,xim : Floating_Number;

  begin
    get(file,xre);
    get(file,xim);
    x := Create(xre,xim);
    Clear(xre);
    Clear(xim);
  end get;

  procedure get ( s : in string; c : in out Complex_Number;
                  last : out integer ) is

    x,y : Floating_Number;

  begin
    get(s,x,last);
    while s(last) /= ' ' loop  -- make sure to start at dividing spaces
      last := last + 1;
      exit when last >= s'last;
    end loop;
    get(s(last..s'last),y,last);
    c := Create(x,y);
    Clear(x); Clear(y);
  end get;

  procedure put ( x : in Complex_Number ) is
  begin
    put(Standard_Output,x);
  end put;

  procedure put ( file : in file_type; x : in Complex_Number ) is

    xre,xim : Floating_Number;

  begin
    xre := REAL_PART(x);
    xim := IMAG_PART(x);
    put(file,xre);
    put(file,"  ");
    put(file,xim);
    Clear(xre);
    Clear(xim);
  end put;

  procedure put ( c : in Complex_Number; fore,aft,exp : in natural32 ) is
  begin
    put(Standard_Output,c,fore,aft,exp);
  end put;

  procedure put ( file : in file_type;
                  c : in Complex_Number; fore,aft,exp : in natural32 ) is

    cre,cim : Floating_Number;

  begin
    cre := REAL_PART(c);
    cim := IMAG_PART(c);
    put(file,cre,fore,aft,exp);
    put(file,"  ");
    put(file,cim,fore,aft,exp);
    Clear(cre);
    Clear(cim);
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

  function Character_Size ( c : Complex_Number ) return natural32 is

    res : natural32;
    x : Floating_Number := REAL_PART(c);
    y : Floating_Number := IMAG_PART(c);

  begin
    res := Character_Size(x) + 2 + Character_Size(y);
    Clear(x);
    Clear(y);
    return res;
  end Character_Size;

  procedure put ( s : out string; c : in Complex_Number ) is

    x : Floating_Number := REAL_PART(c);
    y : Floating_Number := IMAG_PART(c);
    sx : string(1..integer(Character_Size(x)));
    sy : string(1..integer(Character_Size(y)));
  
  begin
    put(sx,x); Clear(x);
    put(sy,y); Clear(y);
    s := sx & "  " & sy;
  end put;

end Multprec_Complex_Numbers_io;
