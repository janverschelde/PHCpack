with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers_io;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;

package body Strings_and_Numbers is

  function Is_Imag ( c : Standard_Complex_Numbers.Complex_Number )
                   return boolean is
  begin
    return ( Standard_Complex_Numbers.REAL_PART(c) = 0.0 );
  end Is_Imag;

  function Is_Imag ( c : Multprec_Complex_Numbers.Complex_Number )
                   return boolean is

    re : constant Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    zero : constant Floating_Number := Create(0.0);

  begin
    return Equal(re,zero);
  end Is_Imag;

  function Is_Real ( c : Standard_Complex_Numbers.Complex_Number )
                   return boolean is
  begin
    return ( Standard_Complex_Numbers.IMAG_PART(c) = 0.0 );
  end is_real;

  function Is_Real ( c : Multprec_Complex_Numbers.Complex_Number )
                   return boolean is

    im : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);
    res : constant boolean := Equal(im,0.0);

  begin
    Clear(im);
    return res;
  end is_real;

  function Is_Integer ( f : double_float ) return boolean is
  begin
    return ( (f - double_float(integer(f))) = 0.0 );
  exception
    when constraint_error => return false;
  end is_integer;

  function Is_Unit ( c : Standard_Complex_Numbers.Complex_Number )
                   return boolean is

    use Standard_Complex_Numbers;

  begin
    if not Is_Real(c) then
      return false;
    elsif not Is_Integer(REAL_PART(c)) then
      return false;
    else
      declare
        i : constant integer := integer(REAL_PART(c));
      begin
        if (i = 1) or (i = -1)
         then return true;
         else return false;
        end if;
      end;
    end if;     
  end Is_Unit;

  function Is_Unit ( c : Multprec_Complex_Numbers.Complex_Number )
                   return boolean is
  begin
    if not Is_Real(c) then
      return false;
    else
      declare
        re : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
        res : constant boolean := Equal(re,1.0);
      begin
        Clear(re);
        return res;
      end;
    end if; 
  end Is_Unit;

  function Sign ( c : Standard_Complex_Numbers.Complex_Number )
                return integer is
  begin
    return integer(Standard_Complex_Numbers.REAL_PART(c));
  end Sign;

  function Sign ( c : Multprec_Complex_Numbers.Complex_Number )
                return integer is
   
    res : integer;
    re : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);

  begin
    if re > 0.0 then
      res := +1;
    elsif Equal(re,0.0) then
      res := 0;
    else
      res := -1;
    end if;
    Clear(re);
    return res;
  end Sign;

  function Truncate ( f : double_float ) return integer32 is

    tol : constant double_float := 1.0E-8;
    i : integer32 := integer32(f);

  begin
    if i >= 0 then
      if double_float(i) > f + tol
       then i := i-1;
      end if;
 --    else if double_float(i) < f - tol --          then i := i+1;
 --         end if;
    end if;
    return i;
  end Truncate;

  function Trim_Zeros ( s : string ) return string is
  begin
    for i in reverse s'range loop
      if s(i) /= '0'
       then return s(s'first..i);
      end if;
    end loop;
    return s(s'first..s'first);
  end Trim_Zeros;

  --function MyConvert ( f : double_float ) return string is

  -- NOTE: an earlier version omitted the abs() in the r1
  --   and crashed for negative numbers.
  --   The other problem with this conversion is the inaccuracy:
  --   -2.3 gets converted to the string "-2.2999999999999998".
  --   Therefore, the Convert below is used.

  --  a : constant integer := Truncate(f);
  --  r1 : constant double_float := abs(f - double_float(a));
  --  re8a : constant double_float := r1*1.0E+8;
  --  b1 : constant integer := Truncate(re8a);
  --  sb1 : constant string := Convert(b1); 
  --  s8b1 : string(1..8);
  -- -- dif1 : constant natural := 8 - sb1'last;
  --  dif1 : constant integer := 8 - sb1'last;
  --  r2 : constant double_float := re8a - double_float(b1);
  --  re8b : constant double_float := r2*1.0E+8;
  --  b2 : constant integer := Truncate(re8b);

  --begin
  --  if dif1 <= 0 then
  --    s8b1 := sb1;
  --  else
  --    s8b1(1..dif1) := (1..dif1 => '0');
  --    for i in sb1'range loop
  --      s8b1(dif1+i) := sb1(i);
  --    end loop;
  --  end if;
  --  if b2 <= 0 then
  --    return Convert(a) & "." & Trim_Zeros(s8b1);
  --  else
  --    declare
  --      sb2 : constant string := Convert(b2);
  --      s8b2 : string(1..8);
  --      dif2 : constant natural := 8 - sb2'last;
  --    begin
  --      if dif2 = 0 then
  --        return Convert(a) & "." & s8b1 & Trim_Zeros(sb2);
  --      else
  --        s8b2(1..dif2) := (1..dif2 => '0');
  --        for i in sb2'range loop
  --          s8b2(dif2+i) := sb2(i);
  --        end loop;
  --        return Convert(a) & "." & s8b1 & Trim_Zeros(s8b2);
  --      end if;
  --    end;
  --  end if;
  --end MyConvert;

  function Convert ( f : double_float ) return string is

    e : constant natural := 25;
    s : string(1..e) := (1..e => ' ');
    k : natural;

  begin
    Standard_Floating_Numbers_io.put(s,f);
    k := s'first;
    while (s(k) = ' ') loop      -- to trim leading spaces
      k := k+1;
      exit when (k > s'last);
    end loop;
    declare
      res : constant string(1..e-k+1) := s(k..e);
    begin
      return res;
    end;
  end Convert;

  function Trim_Backspaces ( s : string ) return string is

  -- DESCRIPTION :
  --   In order to patch the wrong "Character_Size" for certain numbers,
  --   like "3.0", this function returns a string with trailing spaces
  --   removed.  This is a patch to the function "Character_Size"
  --   of Multprec_Floating_Numbers.
 
    slast : integer := s'last;

  begin
    while s(slast) = ' ' and slast >= s'first loop
      slast := slast - 1;
    end loop;
    return s(s'first..slast);
  end Trim_Backspaces;

  function Convert ( f : Floating_Number ) return string is

    nc : constant natural32 := Character_Size(f);
    sf : String(1..integer(nc)) := (1..integer(nc) => ' ');

  begin
    put(sf,f);
    return Trim_Backspaces(sf);
  end Convert;

  function Convert
             ( c : Standard_Complex_Numbers.Complex_Number ) return string is

    use Standard_Complex_Numbers;

    re : constant double_float := REAL_PART(c);
    im : constant double_float := IMAG_PART(c);

  begin
    if im >= 0.0
     then return "(" & Convert(re) & " + " & Convert(im) & "*i)";
     else return "(" & Convert(re) & Convert(im) & "*i)";
    end if;
  end Convert;

  function Convert
             ( c : Multprec_Complex_Numbers.Complex_Number ) return string is

    use Multprec_Complex_Numbers;

    re : Floating_Number := REAL_PART(c);
    im : Floating_Number := IMAG_PART(c);

  begin
    if im > 0.0 or Equal(im,0.0) then
      declare
        res : constant string 
            := "(" & Convert(re) & " + " & Convert(im) & "*i)";
      begin
        Clear(re); Clear(im);
        return res;
      end;
    else 
      declare
        res : constant string := "(" & Convert(re) & Convert(im) & "*i)";
      begin
        Clear(re); Clear(im);
        return res;
      end;
    end if;
  end Convert;

  function Signed_Constant
             ( c : Standard_Complex_Numbers.Complex_Number ) return string is

    use Standard_Complex_Numbers;

    re : double_float;

  begin
    if Is_Real(c) then
      re := REAL_PART(c);
      if Is_Integer(re) then
        declare
          i : constant integer32 := integer32(re);
        begin
          if i >= 0
           then return " + " & Convert(i);
           else return " - " & Convert(-i);
          end if;
        end;
      else
        if re >= 0.0
         then return " + " & Convert(re);
         else return " - " & Convert(abs(re));
        end if;
      end if;
    else
      return "+" & Convert(c);
    end if;
  end Signed_Constant;

  function Signed_Constant
             ( c : Multprec_Complex_Numbers.Complex_Number ) return string is

    re : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    im : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);

  begin
    if not Equal(im,0.0) then -- c is not a real number
      declare
        res : constant string := "+" & Convert(c);
      begin
        Clear(re); Clear(im);
        return res;
      end;
    else                      -- c is a real number
      if re > 0.0 then
        declare
          res : constant string := "+" & Convert(re);
        begin
          Clear(re); Clear(im);
          return res;
        end;
      else
        Min(re);
        declare
          res : constant string := "-" & Convert(re);
        begin
          Clear(re); Clear(im);
          return res;
        end;
      end if;
    end if;
  end Signed_Constant;

  function Unsigned_Constant
             ( c : Standard_Complex_Numbers.Complex_Number ) return string is

    use Standard_Complex_Numbers;

    re : double_float;

  begin
    if Is_Real(c) then
      re := REAL_PART(c);
      if Is_Integer(re) then
        declare
          i : constant integer32 := integer32(re);
        begin
          if i >= 0
           then return Convert(i);
           else return " - " & Convert(-i);
          end if;
        end;
      else
        if re >= 0.0
         then return Convert(re);
         else return " - " & Convert(abs(re));
        end if;
      end if;
    else
      return Convert(c);
    end if;
  end Unsigned_Constant;

  function Unsigned_Constant
             ( c : Multprec_Complex_Numbers.Complex_Number ) return string is

    re : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    im : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);

  begin
    if not Equal(im,0.0) then   -- c is not a real number
      declare
        res : constant string := Convert(c);
      begin
        Clear(re); Clear(im);
        return res;
      end;
    else                        -- c is a real number
      if re > 0.0 then
        declare
          res : constant string := Convert(re);
        begin
          Clear(re); Clear(im);
          return res;
        end;
      else
        Min(re);
        declare
          res : constant string := " - " & Convert(re);
        begin
          Clear(re); Clear(im);
          return res;
        end;
      end if;
    end if;
  end Unsigned_Constant;

  function Signed_Coefficient
             ( c : Standard_Complex_Numbers.Complex_Number ) return string is

    use Standard_Complex_Numbers;

    re : double_float;

  begin
    if Is_Real(c) then
      re := REAL_PART(c);
      if Is_Integer(re) then
        declare
          i : constant integer32 := integer32(re);
        begin
          if i >= 0 then
            if i = 1
             then return " + ";
             else return " + " & Convert(i);
            end if;
          else
            if i = -1
             then return " - ";
             else return " - " & Convert(-i);
            end if;
          end if;
        end;
      else
        if re >= 0.0
         then return " + " & Convert(re);
         else return " - " & Convert(abs(re));
        end if;
      end if;
    else
      return " + " & Convert(c);
    end if;
  end Signed_Coefficient;

  function Signed_Coefficient
             ( c : Multprec_Complex_Numbers.Complex_Number ) return string is

    re : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    im : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);

  begin
    if not Equal(im,0.0) then  -- c is not a real number
      declare
        res : constant string := " + " & Convert(c);
      begin
        Clear(re); Clear(im);
        return res;
      end;
    else                          -- c is a real number
      if re > 0.0 then
        declare
          res : constant string := " + " & Convert(re);
        begin
          Clear(re); Clear(im);
          return res;
        end;
      else
        Min(re);
        declare
          res : constant string := " - " & Convert(re);
        begin
          Clear(re); Clear(im);
          return res;
        end;
      end if;
    end if;
  end Signed_Coefficient;

  function Unsigned_Coefficient
             ( c : Standard_Complex_Numbers.Complex_Number ) return string is

    use Standard_Complex_Numbers;

    re : double_float;

  begin
    if Is_Real(c) then
      re := REAL_PART(c);
      if Is_Integer(re) then
        declare
          i : constant integer32 := integer32(re);
        begin
          if i >= 0 then
            if i = 1
             then return "";
             else return Convert(i);
            end if;
          elsif i = -1 then
            return " - ";
          else
            return " - " & Convert(-i);
          end if;
        end;
      elsif re >= 0.0 then
        return Convert(re);
      else
        return " - " & Convert(abs(re));
      end if;
    else
      return " + " & Convert(c);
    end if;
  end Unsigned_Coefficient;

  function Unsigned_Coefficient
             ( c : Multprec_Complex_Numbers.Complex_Number ) return string is

    re : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    im : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);

  begin
    if not Equal(im,0.0) then
      declare
        res : constant string := " + " & Convert(c);
      begin
        Clear(re); Clear(im);
        return res;
      end;
    else
      if re > 0.0 then
        declare
          res : constant string := Convert(re);
        begin
          Clear(re); Clear(im);
          return res;
        end;
      else
        Min(re);
        declare
          res : constant string := " - " & Convert(re);
        begin
          Clear(re); Clear(im);
          return res;
        end;
      end if;
    end if;
  end Unsigned_Coefficient;

end Strings_and_Numbers;
