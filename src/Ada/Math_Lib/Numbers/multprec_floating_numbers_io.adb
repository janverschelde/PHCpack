with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Natural_Coefficients;      use Multprec_Natural_Coefficients;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Multprec_Parse_Numbers;

package body Multprec_Floating_Numbers_io is

-- NOTE :
--   No exceptions are raised when the input format is incorrect.

-- AUXILIARIES FOR OUTPUT :

  function Head ( i : Integer_Number ) return integer32 is

  -- DESCRIPTION :
  --   Returns the leading decimal number of the i, can be negative.

    res : integer32;
    wrk : Integer_Number;

  begin
    if Multprec_Integer_Numbers.Positive(i) then
      Copy(i,wrk);
      while wrk > 9 loop
        Div(wrk,10);
      end loop;
      res := Create(wrk);
      Clear(wrk);
    elsif Negative(i) then
      Copy(i,wrk);
      while wrk < -9 loop
        Div(wrk,10);
      end loop;
      res := Create(wrk);
      Clear(wrk);
    else
      res := 0;
    end if;
    return res;
  end Head;

  function Tail ( i : Integer_Number ) return Natural_Number is

  -- DESCRIPTION :
  --   Returns the decimals after the leading decimal of i.
  --   The number on return has the same size as i.

    res : Natural_Number;
    res_rep : Array_of_Naturals(0..Size(i));
    acc,wrk,prod : Integer_Number;
    cnt : natural32 := 0;
    r : integer32;

  begin
    if Multprec_Integer_Numbers.Positive(i) then
      Copy(i,wrk);
      while wrk > 9 loop
        cnt := cnt+1;
        Div(wrk,10,r);
        if r /= 0 then
          prod := Create(r);
          for i in 1..cnt-1 loop
            Mul(prod,10);
          end loop;
          Add(acc,prod);
          Clear(prod);
        end if;
      end loop;
      Clear(wrk);
    elsif Negative(i) then
      Copy(i,wrk);
      while wrk < -9 loop
        cnt := cnt+1;
        Div(wrk,10,r);
        if r /= 0 then
          prod := Create(r);
          for i in 1..cnt-1 loop
            Mul(prod,10);
          end loop;
          Add(acc,prod);
          Clear(prod);
        end if;
      end loop;
      Clear(wrk);
    else
      acc := Create(integer32(0));
    end if;
    for i in res_rep'range loop
      res_rep(i) := Coefficient(acc,i);
    end loop;
    Clear(acc);
    res := Create(res_rep);
    return res;
  end Tail;

  procedure Write_Number ( file : in file_type; n : in natural32;
                           cnt : in out natural32 ) is
  begin
    if n < 10 then
      put(file,n,1);
      cnt := cnt - 1;
    else
      Write_Number(file,n/10,cnt);
      if cnt > 0 then
        put(file,n mod 10,1);
        cnt := cnt - 1;
      end if;
    end if;
  end Write_Number;

  procedure Write_Block ( file : in file_type; n : in natural32;
                          cnt : in out natural32 ) is

  -- DESCRIPTION :
  --   This procedure writes the leading zeros, not exceeding cnt.

    expo : constant natural32 := Multprec_Natural_Coefficients.Exponent;
    nbz,acc : natural32 := 0;

  begin
    if n = 0 then
      for i in 1..expo loop
        put(file,"0");
        cnt := cnt - 1;
        exit when (cnt = 0);
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
        cnt := cnt - 1;
        exit when (cnt = 0);
      end loop;
      if cnt > 0
       then Write_Number(file,n,cnt);
      end if;
    end if;
  end Write_Block;

  procedure Write_Zero_Block ( file : in file_type; cnt : in out natural32 ) is

  -- DESCRIPTION :
  --   Writes as many zeros as there are in one block, not exceeding cnt.

    expo : constant natural32 := Multprec_Natural_Coefficients.Exponent;

  begin
    for i in 1..expo loop
      put(file,"0");
      cnt := cnt - 1;
      exit when (cnt = 0);
    end loop;
  end Write_Zero_Block;

  procedure put ( file : in file_type;
                  n : in Natural_Number; dp : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the natural number n using dp decimal places.
  --   If n is too long, then the leading dp decimal places will be written,
  --   otherwise zeros will be added.

    deciplan : constant natural32 := Decimal_Places(n);
    first : boolean;
    coeff,cnt : natural32;

  begin
    if deciplan <= dp then
      put(file,n);
      for i in 1..(dp-deciplan) loop
        put(file,"0");
      end loop;
    else
      first := true;
      cnt := dp;
      for i in reverse 0..Size(n) loop
        coeff := Coefficient(n,i);
        if coeff /= 0 then
          if first then
            Write_Number(file,coeff,cnt);
            first := false;
          else
            Write_Block(file,coeff,cnt);
          end if;
        elsif not first then
          Write_Zero_Block(file,cnt);
     -- else skip leading zeros in representation
        end if;
        exit when (cnt = 0);
      end loop;
    end if;
  end put;

-- INPUT ROUTINES :

  procedure get ( f : in out Floating_Number ) is
  begin
    get(Standard_Input,f);
  end get;

  procedure get ( f : in out Floating_Number; c : in out character ) is
  begin
    get(Standard_Input,f,c);
  end get;

  procedure get ( file : in file_type; f : in out Floating_Number;
                  c : in out character ) is

    fraction,decifrac,exponent : Integer_Number;
    shifexpo : integer32 := 0;
    deciplac : natural32;
    nonneg : boolean;

  begin
    get(file,c,fraction,nonneg);
    if c = '.' then
      get(file,c);
      while c = '0' loop
        Mul(fraction,10);
        shifexpo := shifexpo - 1;
        exit when End_of_Line(file);
        get(file,c);
      end loop;
      if Convert(c) < 10 then
        get(file,c,decifrac);
        deciplac := Decimal_Places(decifrac);
        shifexpo := shifexpo - integer32(deciplac);
        for i in 1..deciplac loop
          Mul(fraction,10);
        end loop;
        if Negative(fraction)
         then Sub(fraction,decifrac);
         else Add(fraction,decifrac);
        end if;
      end if;
    end if;
    if c = 'E' or c = 'e' then
      c := ' ';
      get(file,c,exponent);
      if Equal(fraction,0) then
        Clear(exponent);
        exponent := Create(integer32(0));
      elsif shifexpo /= 0 then
        Add(exponent,shifexpo);
      end if;
    elsif Equal(fraction,0) then
      exponent := Create(integer32(0));
    else
      exponent := Create(shifexpo);
    end if;
    f := Create(fraction,exponent);
  end get;

  procedure get ( file : in file_type; f : in out Floating_Number ) is

    c : character := ' ';

  begin
    get(file,f,c);
  end get;

  procedure get ( s : in string; f : in out Floating_Number;
                  last : out integer ) is

    pos : integer := s'first;

  begin
    Multprec_Parse_Numbers.parse(s,pos,f);
    last := pos;
  end get;

-- OUTPUT ROUTINES :

  procedure put ( f : in Floating_Number ) is
  begin
    put(Standard_Output,f);
  end put;

  procedure put ( file : in file_type; f : in Floating_Number ) is

    frac,expo : Integer_Number;
    tafr : Natural_Number;
    decifrac,decitafr : natural32;
    hd : integer32;

  begin
    Copy(Fraction(f),frac);  
    if Equal(frac,0) then
      put(file,"0");
    else
      hd := Head(frac);
      put(file,hd,1); put(file,".");
      tafr := Tail(Fraction(f));
      decifrac := Decimal_Places(frac);
      decitafr := Decimal_Places(tafr);
      if not Equal(tafr,0) then
        for i in 1..decifrac-decitafr-1 loop
          put(file,"0");
        end loop;
        put(file,tafr);
      else
        put(file,"0");
      end if;
      expo := Exponent(f) + integer32(decifrac - 1);
      if not Equal(expo,0) then
        put(file,"E");
        if expo > 0
         then put(file,"+");
        end if;
        put(file,expo);
      end if;
      Clear(tafr);
      Clear(expo);
    end if;
    Clear(frac);
  end put;

  procedure put ( f : in Floating_Number; fore,aft,exp : in natural32 ) is
  begin
    put(Standard_Output,f,fore,aft,exp);
  end put;

  procedure put ( file : in file_type;
                  f : in Floating_Number; fore,aft,exp : in natural32 ) is

    frac,expo : Integer_Number;
    tafr : Natural_Number;
    decifrac,decitafr,deciexpo,cnt : natural32;
    hd : integer32;

  begin
    Copy(Fraction(f),frac);
    if Equal(frac,0) then
      for i in 1..(fore-1) loop
        put(file," ");
      end loop;
      put(file,"0.");
      for i in 1..aft loop
        put(file,"0");
      end loop;
      put(file,"E+");
      for i in 1..(exp-1) loop
        put(file,"0");
      end loop;
    else
      hd := Head(frac);
      if hd > 0 then
        for i in 1..(fore-1) loop
          put(file," ");
        end loop;
      else
        for i in 1..(fore-2) loop
          put(file," ");
        end loop;
      end if;
      put(file,hd,1); put(file,".");
      tafr := Tail(Fraction(f));
      decifrac := Decimal_Places(frac);
      decitafr := Decimal_Places(tafr);
      if Equal(tafr,0) then
        for i in 1..aft loop
          put(file,"0");
        end loop;
      else
        cnt := aft;
        for i in 1..decifrac-decitafr-1 loop
          put(file,"0");
          cnt := cnt-1;
          exit when (cnt = 0);
        end loop;
        if cnt /= 0
         then put(file,tafr,cnt);
        end if;
      end if;
      expo := Exponent(f) + integer32(decifrac - 1);
      deciexpo := Decimal_Places(expo);
      put(file,"E");
      if Equal(expo,0) then
        put(file,"+");
        for i in 1..(exp-1) loop
          put(file,"0");
        end loop;
      else
        if expo > 0
         then put(file,"+");
         else put(file,"-");
        end if;
        for i in 1..(exp-deciexpo-1) loop
          put(file,"0");
        end loop;
        put(file,Unsigned(expo));
      end if;
      Clear(tafr);
      Clear(expo);
    end if;
    Clear(frac);
  end put;

  procedure put ( f : in Floating_Number; dp : in natural32 ) is
  begin
    put(f,dp,dp,dp);
  end put;

  procedure put ( file : in file_type;
                  f : in Floating_Number; dp : in natural32 ) is
  begin
    put(file,f,dp,dp,dp);
  end put;

  function Character_Size ( f : Floating_Number ) return natural32 is

    dpf : natural32 := Decimal_Places(Unsigned(Fraction(f)));
    ef : Integer_Number := Exponent(f) + integer32(dpf-1);
    dpe : constant natural32 := Decimal_Places(Unsigned(ef));
    res : natural32;

  begin
    if dpf = 1
     then dpf := dpf + 1;
    end if;
    if Equal(fraction(f),0) then
      return 1;
    elsif Negative(Fraction(f)) then
      res := dpf + 2;
    else
      res := dpf + 1;
    end if;
    if not Equal(ef,0)
     then res := res + dpe + 2;
    end if;
    Clear(ef);
    return res;
  end Character_Size;

  procedure put ( s : out string; f : in Floating_Number ) is

    dp : constant natural32 := Decimal_Places(Unsigned(Fraction(f)));
    ef : constant Integer_Number := Exponent(f) + integer32(dp-1);
    hd : constant integer32 := Head(Fraction(f));
    tl : constant Natural_Number := Tail(Fraction(f));
    ind,dpf : natural;

  begin
    if dp = 1
     then dpf := 2;
     else dpf := natural(dp);
    end if;
    if Equal(fraction(f),0) then
      s(s'first) := '0';
    else
      if hd < 0 then
        s(s'first) := '-';
        s(s'first+1) := Convert_Decimal(natural32(-hd));
        s(s'first+2) := '.';
        if Equal(tl,0) then
          s(s'first+3) := '0';
          ind := s'first+4;
        else
          ind := s'first+3;
          for i in 1..dp-Decimal_Places(tl)-1 loop
            s(ind) := '0';
            ind := ind + 1;
          end loop;
          put(s(ind..s'first+dpf+1),tl);
          ind := s'first+dpf+2;
        end if;
      else
        s(s'first) := Convert_Decimal(natural32(hd));
        s(s'first+1) := '.';
        if Equal(tl,0) then
          s(s'first+2) := '0';
          ind := s'first+3;
        else
          ind := s'first+2;
          for i in 1..dp-Decimal_Places(tl)-1 loop
            s(ind) := '0';
            ind := ind + 1;
          end loop;
          put(s(ind..s'first+dpf),tl);
          ind := s'first+dpf+1;
        end if;
      end if;
      if not Equal(ef,0) then
        s(ind) := 'E';
        if Negative(ef)
         then put(s(ind+1..s'last),ef);
         else s(ind+1) := '+';
              put(s(ind+2..s'last),ef);
        end if;
      end if;
    end if;
  end put;

end Multprec_Floating_Numbers_io;
