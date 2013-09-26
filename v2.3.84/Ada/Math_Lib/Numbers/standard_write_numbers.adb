with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;

package body Standard_Write_Numbers is

  function Is_Imag ( c : Complex_Number ) return boolean is
  begin
    return ( REAL_PART(c) = 0.0 );
  end Is_Imag;

  function Is_Real ( c : Complex_Number ) return boolean is
  begin
    return ( IMAG_PART(c) = 0.0 );
  end Is_Real;

  function Is_Integer ( f : double_float ) return boolean is
  begin
    return ( (f - double_float(integer(f))) = 0.0 );
  exception
    when constraint_error => return false;
  end Is_Integer;

  function Length ( n : integer32 ) return natural32 is
  begin
    for i in 1..10 loop
      if n < integer32(10.0**i) then
        return natural32(i);
      end if;
    end loop;
    return 11;
  end Length;
 
  procedure Write_Number ( file : in file_type; i : in integer32;
                           cnt : out natural32 ) is
  begin
    put(file,i,1);
    cnt := Length(i);
  end Write_Number;

  function Length ( n : double_float ) return natural32 is
  begin
    if Is_Integer(n)
     then return Length(integer32(n));
     else return 21;
    end if;
  end Length;

  procedure Write_Number ( file : in file_type; f : in double_float;
                           cnt : out natural32 ) is
  begin
    if Is_Integer(f)
     then Write_Number(file,integer32(f),cnt);
     else put(file,f); cnt := 21;
    end if;
  end Write_Number;

  function Length ( n : Complex_Number ) return natural32 is

    res : natural32;

  begin
    if Is_Real(n) then
      res := Length(n);
    elsif Is_Imag(n) then
      res := Length(n) + 2;
    else
      res := 4 + Length(REAL_PART(n));
      if IMAG_PART(n) = 1.0 then
        res := res + 1;
      elsif IMAG_PART(n) = -1.0 then
        res := res + 3;
      else
        res := res + Length(IMAG_PART(n)) + 2;
      end if;
    end if;
    return res;
  end Length;

  procedure Write_Number ( file : in file_type; c : in Complex_Number;
                           cnt : out natural32 ) is

    nc : natural32;

  begin
    if Is_Real(c) then
      Write_Number(file,REAL_PART(c),cnt);
    elsif Is_Imag(c) then
      Write_Number(file,IMAG_PART(c),cnt);
      put(file,"*i"); cnt := cnt + 2;
    else
      put(file,"("); cnt := 1;
      Write_Number(file,REAL_PART(c),nc);
      if IMAG_PART(c) > 0.0
       then put(file," +");
       else put(file," -");
      end if;
      cnt := cnt + nc + 2;
      if IMAG_PART(c) = 1.0 then
        put(file,"i"); cnt := cnt + 1;
      elsif IMAG_PART(c) = -1.0 then
        put(file,"i"); cnt := cnt + 1;
      else
        Write_Number(file,abs(IMAG_PART(c)),nc);
        put(file,"*i"); cnt := cnt + nc + 2;
      end if;
      put(file,")"); cnt := cnt + 1;
    end if;
  end Write_Number;

  procedure Write_Coefficient ( file : in file_type; c : in Complex_Number;
                                cnt : out natural32 ) is
  begin
    if Equal(c,Create(-1.0)) then
      put(file,'-'); cnt := 1;
    elsif Equal(c,Create(0.0,1.0)) then
      put(file,"i*"); cnt := 2;
    elsif Equal(c,Create(0.0,-1.0)) then
      put(file,"-i*"); cnt := 3;
    elsif not Equal(c,Create(1.0)) then
      Write_Number(file,c,cnt); put(file,'*'); cnt := cnt + 1;
    else -- the case where c = +1, nothing is written to file
      cnt := 0;
    end if;
  end Write_Coefficient;

  procedure Write_Plus ( file : in file_type; c : in Complex_Number;
                         cnt : out natural32 ) is
  begin
    if (Is_Real(c) and then REAL_PART(c) > 0.0)
      or else (Is_Imag(c) and then IMAG_PART(c) > 0.0)
      or else (not Is_Real(c) and then not Is_Imag(c))
     then cnt := 1; put(file,'+');
     else cnt := 0;
    end if;
  end Write_Plus;

end Standard_Write_Numbers;
