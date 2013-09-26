with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Characters_and_Numbers;            use Characters_and_Numbers;

package body Multprec_Natural64_Numbers_io is

-- IMPORTANT NOTICE :
--   The choice of base is assumed to be decimal.
--   The maximum number of digits is one thousand times expo, see "maxl".

-- CONSTANTS :

  expo : constant natural32 := Multprec_Natural64_Numbers.Exponent;
  maxl : constant natural32 := 1000;

-- DATA STRUCTURE :

  type Array_of_Strings is
    array ( natural32 range <> ) of String(1..natural(expo));

-- BASIC PRIMITIVES FOR INPUT/OUTPUT :

  function Convert ( s : Array_of_Strings ) return Array_of_Naturals is

    res : Array_of_Naturals(s'range) := (s'range => 0);

  begin
    for i in reverse s'range loop
      res(res'last-i) := Convert(s(i));
    end loop;
    return res;
  end Convert;

  function Size ( len : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Given the number of characters read, the size of the natural number
  --   will be determined.

    res : constant natural32 := len/expo;

  begin
    if res = 0 then
      return res;
    else
      if expo*res = len
       then return res-1;
       else return res;
      end if;
    end if;
  end Size;

  function Create ( l : in natural32; s : in String ) return Array_of_Strings is

  -- DESCRIPTION :
  --   Partitions the string in blocks, according to the base.

    res : Array_of_Strings(0..l);
    ind : natural32 := l;
    cnt : natural32 := 0;

  begin
    for i in res'range loop
      res(i) := (res(i)'range => ' ');
    end loop;
    for i in reverse s'range loop
      cnt := cnt + 1;
      if cnt <= expo then
        res(ind)(natural(expo-cnt+1)) := s(i);
      else
        ind := ind-1;
        cnt := 1;
        res(ind)(natural(expo-cnt+1)) := s(i);
      end if;
    end loop;
    return res;
  end Create;

  procedure Write_Block ( file : in file_type; n : in natural64 ) is

  -- DESCRIPTION :
  --   This procedure writes the leading zeros.

    nbz : natural32 := 0;
    acc : natural64 := 0;

  begin
    if n = 0 then
      for i in 1..expo loop
        put(file,"0");
      end loop;
    else
      acc := 10;
      for i in 1..(expo-1) loop
        if n < acc
         then nbz := expo-i;
         else acc := acc*10;
        end if;
        exit when (nbz /= 0);
      end loop;
      for i in 1..nbz loop
        put(file,"0");
      end loop;
      put(file,n,1);
    end if;
  end Write_Block;

  procedure Write_Zero_Block ( file : in file_type ) is
  begin
    for i in 1..expo loop
      put(file,"0");
    end loop;
  end Write_Zero_Block;

-- INPUT ROUTINES :

  procedure get ( file : in file_type;
                  lc : in out character; n : in out Natural_Number ) is

    s : String(1..natural(maxl));
    cvn : natural32 := Convert(lc);
    cnt : natural := 0;

  begin
    while cvn < 10 loop
      cnt := cnt+1;
      s(cnt) := lc;
      exit when End_of_Line(file) or (cnt = s'last);
     -- Skip_Spaces(file,lc);
      Skip_Underscores(file,lc);
      cvn := Convert(lc);
    end loop;
    declare
      sz : constant natural32 := Size(natural32(cnt));
      sn : constant Array_of_Strings(0..sz) := Create(sz,s(1..cnt));
      an : constant Array_of_Naturals(0..sz) := Convert(sn);
    begin
      Clear(n);
      n := Create(an);
    end;
  end get;

  procedure get ( lc : in out character; n : in out Natural_Number ) is
  begin
    get(Standard_Input,lc,n);
  end get;

  procedure get ( n : in out Natural_Number ) is
  begin
    get(Standard_Input,n);
  end get;

  procedure get ( file : in file_type; n : in out Natural_Number ) is

    c : character := ' ';

  begin
    Skip_Spaces(file,c);
    get(file,c,n);
  end get;

  function Coefficient_Size ( s : string ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of coefficients the string s may hold.

    res : constant natural32 := natural32(s'last)/expo;

  begin
    return res;
  end Coefficient_Size;

  procedure get ( s : in string; n : in out Natural_Number ) is

    nb : constant natural32 := Coefficient_Size(s);
    coeff : Array_of_Naturals(0..nb);
    cffind : natural32 := 0;
    buffer : string(1..natural(expo));
    bufind : natural := natural(expo+1);

  begin
    coeff(nb) := 0;
    for i in reverse s'range loop
      if bufind > 1 then
        bufind := bufind - 1;
        buffer(bufind) := s(i);
      else
        coeff(cffind) := Convert(buffer);
        cffind := cffind + 1;
        bufind := natural(expo);
        buffer(bufind) := s(i);
      end if;
    end loop;
    coeff(cffind) := Convert(buffer(bufind..natural(expo)));
    n := Create(coeff);
  end get;

-- OUTPUT ROUTINES :

  procedure put ( n : in Natural_Number ) is
  begin
    put(Standard_Output,n);
  end put;

  procedure put ( file : in file_type; n : in Natural_Number ) is

  -- NOTE : the blocks can be separated by underscores.
  --   In principal, other symbols could be used, however, only underscores
  --   are skipped when processing a natural number.
  
    first : boolean := true;    -- first nonzero, leading block still to write
    coeff : natural64;

  begin
    if Empty(n) then
      put(file,"0");
    else
      for i in reverse 0..Size(n) loop
        coeff := Coefficient(n,i);
        if coeff /= 0 then
          if first
           then put(file,coeff,1); first := false;
           else Write_Block(file,coeff);
          end if;
        elsif not first then
          Write_Zero_Block(file);
        -- else skip leading zeros
        end if;
       -- if (not first and (i>0))   -- leading block written and not at end
       --  then put(file,"_");       -- so, write a separator symbol
       -- end if;
      end loop;
      if first 
       then put(file,"0");         -- there was no nonzero block, so n=0.
      end if;
    end if;
  end put;

  procedure put ( s : out string; n : Natural_Number ) is

    ind : integer := s'first-1;
    first : boolean := true;
    coeff : natural64;

  begin
    if Empty(n) then
      s(s'first) := '0';
    else
      for i in reverse 0..Size(n) loop
        coeff := Coefficient(n,i);
        if coeff /= 0 then
          declare
            sn : constant string := nConvert(coeff);
          begin
            if first then
              first := false;
            else
              for i in 1..natural(expo)-sn'last loop
                ind := ind + 1;
                s(ind) := '0';
              end loop;
            end if; 
            for i in sn'range loop
              ind := ind + 1;
              s(ind) := sn(i);
            end loop;
          end;
        elsif not first then
          for i in 1..expo loop
            ind := ind + 1;
            s(ind) := '0';
          end loop;
        end if;
      end loop;
    end if;
  end put;

  procedure put ( n : in Array_of_Naturals ) is
  begin
    put(Standard_Output,n);
  end put;

  procedure put ( file : in file_type; n : in Array_of_Naturals ) is
  begin
    for i in reverse n'range loop
      if n(i) = 0
       then Write_Zero_Block(file);
       else Write_Block(file,n(i));
      end if;
     -- if i > 0
     --  then put(file,"_");
     -- end if;
    end loop;
  end put;

  procedure put ( s : out string; n : in Array_of_Naturals ) is

    ind : integer := s'first-1;

  begin
    for i in reverse n'range loop
      if n(i) = 0 then
        for i in 1..expo loop
          ind := ind + 1;
          s(ind) := '0';
        end loop;
      else
        declare
          sn : constant string := Characters_and_Numbers.nConvert(n(i));
        begin
          for i in 1..natural(expo)-sn'last loop
            ind := ind + 1;
            s(ind) := '0';
          end loop;
          for i in sn'range loop
            ind := ind + 1;
            s(ind) := sn(i);   
          end loop;
        end;
      end if;
    end loop;
  end put;

  procedure put ( n : in Natural_Number; dp : in natural32 ) is
  begin
    put(Standard_Output,n,dp);
  end put;

  procedure put ( file : in file_type;
                  n : in Natural_Number; dp : in natural32 ) is
  begin
    for i in 1..(dp-Decimal_Places(n)) loop
      put(file," ");
    end loop;
    put(file,n);
  end put;

end Multprec_Natural64_Numbers_io;
