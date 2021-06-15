with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Multprec_Natural64_Coefficients;    use Multprec_Natural64_Coefficients;

--with text_io,integer_io; use text_io,integer_io;
--with Multprec_Integer_Numbers_io;  use Multprec_Integer_Numbers_io;

package body Multprec_Floating64_Numbers is

-- NOTES ON THE IMPLEMENTATION :
--   0) See the notes in the bodies of Multprec_{Natural,Integer}_Numbers.
--   1) 0.0 must have zero fraction and zero exponent.
--   2) The truncate and round operations all depend on the radix.
--   3) The function "Decimal_Places" occupies a central place.
--   4) The normalized representation has a nonzero least significant digit.

  rdx : constant natural64 := Multprec_Natural64_Coefficients.Radix;
  tol : constant double_float := 10.0**(-8);
  the_base : constant natural64 := Multprec_Natural64_Coefficients.Base;
  the_expo : constant natural32 := Multprec_Natural64_Coefficients.Exponent;
  sub_base : constant natural64 := the_base/rdx;

-- AUXILIARIES :

  function Truncate ( f : double_float ) return integer is

    i : integer := integer(f);

  begin
    if i >= 0
     then if double_float(i) > f + tol
           then i := i-1;
          end if;
--     else if double_float(i) < f - tol
--           then i := i+1;
--          end if;
    end if;
    return i;
  end Truncate;

  function Max_Size ( i1,i2 : Integer_Number ) return natural32 is

  -- DESCRIPTION :
  --   Return the maximal size of i1 and i2.

    sz1 : constant natural32 := Size(i1);
    sz2 : constant natural32 := Size(i2);

  begin
    if sz1 > sz2
     then return sz1;
     else return sz2;
    end if;
  end Max_Size;

  procedure Normalize ( f : in out Floating_Number ) is

  -- DESCRIPTION :
  --   Normalizes the representation of f so that its fraction has as
  --   rightmost digit a nonzero number.

    ns : natural32;

  begin
    Shift_Right(f.fraction,ns);
    Add(f.exponent,integer64(ns));
   -- if Size(f.fraction) /= 4
   --  then put_line("Number not having size of fraction = 4!!");
   -- end if;
  end Normalize;

  procedure Adjust_Expo ( f : in out double_float;
                          expo : in Integer_Number ) is

  -- DESCRIPTION :
  --   Multiplies f with rdx**expo.

    cnt : Integer_Number;

  begin
    Copy(expo,cnt);
    while cnt > 0 loop
      f := f*double_float(rdx);
      Sub(cnt,1);
    end loop;
    while cnt < 0 loop
      f := f/double_float(rdx);
      Add(cnt,1);
    end loop;
    Clear(cnt);
  end Adjust_Expo;

  function Coeff0 ( i : Integer_Number ) return integer64 is

  -- DESCRIPTION :
  --   Returns the 0-th coefficient of i with the same sign as i.
  --   This is useful when Size(i) = 0.

    res : constant natural64 := Coefficient(i,0);

  begin
    if Negative(i)
     then return -integer64(res);
     else return integer64(res);
    end if;
  end Coeff0;

-- CONSTRUCTORS :

  function Create ( i : integer ) return Floating_Number is
  begin
    return Create(integer64(i));
  end Create;

  function Create ( i : integer64 ) return Floating_Number is

    res : Floating_Number;

  begin
    res.fraction := Create64(i);
    res.exponent := Create64(0);
    Normalize(res);
    return res;
  end Create;

  function Create ( i : Integer_Number ) return Floating_Number is

    res : Floating_Number;

  begin
    res.fraction := +i;
    res.exponent := Create64(0);
    Normalize(res);
    return res;
  end Create;

  function Create ( fraction,exponent : integer64 ) return Floating_Number is

    res : Floating_Number;

  begin
    res.fraction := Create64(fraction);
    res.exponent := Create64(exponent);
    Normalize(res);
    return res;
  end Create;

  function Create ( fraction : Integer_Number;
                    exponent : integer64 ) return Floating_Number is

    res : Floating_Number;

  begin
    res.fraction := +fraction;
    res.exponent := Create64(exponent);
    Normalize(res);
    return res;
  end Create;

  function Create ( fraction : integer64;
                    exponent : Integer_Number ) return Floating_Number is

    res : Floating_Number;

  begin
    res.fraction := Create64(fraction);
    res.exponent := +exponent;
    Normalize(res);
    return res;
  end Create;

  function Create ( fraction,exponent : Integer_Number )
                  return Floating_Number is

    res : Floating_Number;

  begin
    res.fraction := +fraction;
    res.exponent := +exponent;
    Normalize(res);
    return res;
  end Create;

  function Create ( f : double_float ) return Floating_Number is

  -- NOTE :
  --   To avoid rounding errors, comparison operations are prefered
  --   above relying on intermediate results of calculations, although
  --   rounding errors still occur!

    res : Floating_Number;
    ten : Integer_Number;
    wrkflt,fltexp,lowflt,uppflt,incflt : double_float;
    intexp : integer;
    fraction,exponent,intinc : Integer_Number;

  begin
    if f = 0.0 then
      res := Create(integer64(0));
      return res;
    end if;
    if f > 0.0
     then wrkflt := f;
     else wrkflt := -f;
    end if;
    fltexp := LOG10(wrkflt);
    intexp := Truncate(fltexp);
   -- PATCH FOR BUG : 16 instead of 15
    if intexp < 16
     then wrkflt := wrkflt*(10.0**(16-intexp));
          exponent := Create64(integer64(intexp-16));
          intexp := 16;
     else exponent := Create64(0);
    end if;
    lowflt := 10.0**intexp;
    uppflt := 10.0**(intexp+1);
    while lowflt > wrkflt loop
      intexp := intexp-1;
      lowflt := 10.0**intexp;
      uppflt := 10.0**(intexp+1);
    end loop;
    ten := Create64(10);
    fraction := ten**natural64(intexp);
    for k in 1..16 loop
      incflt := 10.0**intexp;
      intinc := ten**natural64(intexp);
      uppflt := lowflt;
      for i in 1..10 loop
        uppflt := uppflt + incflt;
        exit when (uppflt > wrkflt);
        lowflt := uppflt;
        Add(fraction,intinc);
      end loop;
      intexp := intexp - 1;
      Clear(intinc);
    end loop;
    Clear(ten);
    if f < 0.0
     then Min(fraction);
    end if;
   -- res := Create(fraction,exponent);  -- dangerous, makes copy !!!
    res.fraction := fraction;
    res.exponent := exponent;
    Normalize(res);
    return res;
  end Create;

-- SELECTORS :

  function Fraction ( f : Floating_Number ) return Integer_Number is
  begin
    return f.fraction;
  end Fraction;

  function Exponent ( f : Floating_Number ) return Integer_Number is
  begin
    return f.exponent;
  end Exponent;

  function Size_Fraction ( f : Floating_Number ) return natural32 is
  begin
    return Size(Fraction(f));
  end Size_Fraction;

  function Size_Exponent ( f : Floating_Number ) return natural32 is
  begin
    return Size(Exponent(f));
  end Size_Exponent;

  function Decimal_Places_Fraction ( f : Floating_Number ) return natural32 is
  begin
    return Decimal_Places(Fraction(f));
  end Decimal_Places_Fraction;

  function Decimal_Places_Exponent ( f : Floating_Number ) return natural32 is
  begin
    return Decimal_Places(Exponent(f));
  end Decimal_Places_Exponent;

  function Decimal_to_Size ( deci : natural32 ) return natural32 is

    res : natural32;

  begin
    res := deci/the_expo - 1;
    if (res+1)*the_expo < deci
     then res := res+1;
    end if;
    return res;
  end Decimal_to_Size;

  function Size_to_Decimal ( size : natural32 ) return natural32 is
  begin
    return (size+1)*the_expo;
  end Size_to_Decimal;

  function AbsVal ( f : Floating_Number ) return Floating_Number is

    res : Floating_Number;

  begin
    if Negative(f.fraction)
     then res.fraction := -f.fraction;
     else res.fraction := +f.fraction;
    end if;
    res.exponent := +f.exponent;
    return res;
  end AbsVal;

  function Positive ( f : Floating_Number ) return boolean is
  begin
    return Multprec_Integer64_Numbers.Positive(f.fraction);
  end Positive;

  function Negative ( f : Floating_Number ) return boolean is
  begin
    return Multprec_Integer64_Numbers.Negative(f.fraction);
  end Negative;

  function Sign ( f : Floating_Number ) return integer32 is
  begin
    return Multprec_Integer64_Numbers.Sign(f.fraction);
  end Sign;

  function Truncate_to_Nearest_Integer
             ( f : Floating_Number ) return Integer_Number is

    res : Integer_Number;
    e64 : integer64;
    e32 : integer32;

  begin
    if Equal(f.exponent,0) then
      Copy(f.fraction,res);
    elsif f.exponent > 0 then
      Copy(f.fraction,res);
      e64 := Create(f.exponent);
      e32 := integer32(e64);
      Mul_Radix(res,natural32(e32));
    else
      e64 := Create(f.exponent);
      e32 := -integer32(e64);
      if natural32(e32) > Decimal_Places(f.fraction) then
        res := Create32(0);
      else
        Copy(f.fraction,res);
        for i in 1..e32 loop  -- not very efficient, but no Div_Radix yet
          Div(res,10);
        end loop;
      end if;
    end if;
    return res;
  end Truncate_to_Nearest_Integer;
 
-- COMPARISON AND COPYING :

  function Equal ( f1 : Floating_Number; f2 : double_float ) return boolean is

    ff2 : Floating_Number := Create(f2);
    res : constant boolean := Equal(f1,ff2);

  begin
    Clear(ff2);
    return res;
  end Equal;

  function Equal ( f1,f2 : Floating_Number ) return boolean is

    res : boolean;
    f1deci : constant natural32 := Decimal_Places(f1.fraction);
    f2deci : constant natural32 := Decimal_Places(f2.fraction);
    stf1expo,stf2expo : integer64;
    mpf1expo,mpf2expo,mulfra : Integer_Number;

  begin
    if Size(f1.exponent) = 0 and Size(f2.exponent) = 0 then
      stf1expo := Coeff0(f1.exponent) + integer64(f1deci);
      stf2expo := Coeff0(f2.exponent) + integer64(f2deci);
      if stf1expo = stf2expo then
        if f1deci = f2deci then
          res := Equal(f1.fraction,f2.fraction);
        else
          if f1deci < f2deci then
            mulfra := f1.fraction*integer64(rdx);
            for i in 1..(f2deci-f1deci-1) loop
              Mul(mulfra,integer64(rdx));
            end loop;
           -- Mul_Radix(mulfra,f2deci-f1deci-1);
            res := Equal(mulfra,f2.fraction);
          else
            mulfra := f2.fraction*integer64(rdx);
            for i in 1..(f1deci-f2deci-1) loop
              Mul(mulfra,integer64(rdx));
            end loop;
           -- Mul_Radix(mulfra,f1deci-f2deci-1);
            res := Equal(f1.fraction,mulfra);
          end if;
          Clear(mulfra);
        end if;
      else
        res := false;
      end if;
    else
      mpf1expo := f1.exponent + integer64(f1deci);
      mpf2expo := f2.exponent + integer64(f2deci);
      if Equal(mpf1expo,mpf2expo) then
        if f1deci = f2deci then
          res := Equal(f1.fraction,f2.fraction);
        else
          if f1deci < f2deci then
            mulfra := f1.fraction*integer64(rdx);
            for i in 1..(f2deci-f1deci-1) loop
              Mul(mulfra,integer64(rdx));
            end loop;
           -- Mul_Radix(mulfra,f2deci-f1deci-1);
            res := Equal(mulfra,f2.fraction);
          else
            mulfra := f2.fraction*integer64(rdx);
            for i in 1..(f1deci-f2deci-1) loop
              Mul(mulfra,integer64(rdx));
            end loop;
           -- Mul_Radix(mulfra,f1deci-f2deci-1);
            res := Equal(f1.fraction,mulfra);
          end if;
          Clear(mulfra);
        end if;
      else
        res := false;
      end if;
      Clear(mpf1expo); Clear(mpf2expo);
    end if;
    return res;
  end Equal;

  function "<" ( f1 : double_float; f2 : Floating_Number ) return boolean is

    ff1 : Floating_Number := Create(f1);
    res : constant boolean := (ff1 < f2);

  begin
    Clear(ff1);
    return res;
  end "<";

  function "<" ( f1 : Floating_Number; f2 : double_float ) return boolean is

    ff2 : Floating_Number := Create(f2);
    res : constant boolean := (f1 < ff2);

  begin
    Clear(ff2);
    return res;
  end "<";

  function "<" ( f1,f2 : Floating_Number ) return boolean is

    res : boolean;
    f1deci : constant natural32 := Decimal_Places(f1.fraction);
    f2deci : constant natural32 := Decimal_Places(f2.fraction);
    stf1expo,stf2expo : integer64;
    mpf1expo,mpf2expo,mulfra : Integer_Number;

  begin
    if Size(f1.exponent) = 0 and Size(f2.exponent) = 0 then
      stf1expo := Coeff0(f1.exponent) + integer64(f1deci);
      stf2expo := Coeff0(f2.exponent) + integer64(f2deci);
      if (stf1expo < stf2expo) then
        if Sign(f2.fraction) > 0 then
          res := true;
        elsif Sign(f2.fraction) < 0 then
          res := false;
        else
          res := (Sign(f1.fraction) < 0);
        end if;
      elsif (stf1expo > stf2expo) then
        if Sign(f1.fraction) > 0 then
          res := false;
        elsif Sign(f1.fraction) < 0 then
          res := true;
        else
          res := (Sign(f2.fraction) > 0);
        end if;
      elsif f1deci = f2deci then
        res := (f1.fraction < f2.fraction);
      elsif f1deci < f2deci then
        mulfra := f1.fraction*integer64(rdx);
        for i in 1..(f2deci-f1deci-1) loop
          Mul(mulfra,integer64(rdx));
        end loop;
       -- Mul_Radix(mulfra,f2deci-f1deci-1);
        res := (mulfra < f2.fraction);
        Clear(mulfra);
      else
        mulfra := f2.fraction*integer64(rdx);
        for i in 1..(f1deci-f2deci-1) loop
          Mul(mulfra,integer64(rdx));
        end loop;
       -- Mul_Radix(mulfra,f1deci-f2deci-1);
        res := (f1.fraction < mulfra);
        Clear(mulfra);
      end if;
    else
      mpf1expo := f1.exponent + integer64(f1deci);
      mpf2expo := f2.exponent + integer64(f2deci);
      if (mpf1expo < mpf2expo) then
        if Sign(f2.fraction) > 0 then
          res := true;
        elsif Sign(f2.fraction) < 0 then
          res := false;
        else
          res := (Sign(f1.fraction) < 0);
        end if;
      elsif (mpf1expo > mpf2expo) then
        if Sign(f1.fraction) > 0 then
          res := false;
        elsif Sign(f1.fraction) < 0 then
          res := true;
        else
          res := (Sign(f2.fraction) > 0);
        end if;
      elsif f1deci = f2deci then
        res := (f1.fraction < f2.fraction);
      elsif f1deci < f2deci then
        mulfra := f1.fraction*integer64(rdx);
        for i in 1..(f2deci-f1deci-1) loop
          Mul(mulfra,integer64(rdx));
        end loop;
       -- Mul_Radix(mulfra,f2deci-f1deci-1);
        res := (mulfra < f2.fraction);
        Clear(mulfra);
      else
        mulfra := f2.fraction*integer64(rdx);
        for i in 1..(f1deci-f2deci-1) loop
          Mul(mulfra,integer64(rdx));
        end loop;
       -- Mul_Radix(mulfra,f1deci-f2deci-1);
        res := (f1.fraction < mulfra);
        Clear(mulfra);
      end if;
      Clear(mpf1expo); Clear(mpf2expo);
    end if;
    return res;
  end "<";

  function ">" ( f1 : double_float; f2 : Floating_Number ) return boolean is

    ff1 : Floating_Number := Create(f1);
    res : constant boolean := (ff1 > f2);

  begin
    Clear(ff1);
    return res;
  end ">";

  function ">" ( f1 : Floating_Number; f2 : double_float ) return boolean is

    ff2 : Floating_Number := Create(f2);
    res : constant boolean := (f1 > ff2);

  begin
    Clear(ff2);
    return res;
  end ">";

  function ">" ( f1,f2 : Floating_Number ) return boolean is

    res : boolean;
    f1deci : constant natural32 := Decimal_Places(f1.fraction);
    f2deci : constant natural32 := Decimal_Places(f2.fraction);
    stf1expo,stf2expo : integer64;
    mpf1expo,mpf2expo,mulfra : Integer_Number;

  begin
    if Size(f1.exponent) = 0 and Size(f2.exponent) = 0 then
      stf1expo := Coeff0(f1.exponent) + integer64(f1deci);
      stf2expo := Coeff0(f2.exponent) + integer64(f2deci);
      if (stf1expo > stf2expo) then
        if Sign(f1.fraction) > 0 then
          res := true;
        elsif Sign(f1.fraction) < 0 then
          res := false;
        else
          res := (Sign(f2.fraction) < 0);
        end if;
      elsif (stf1expo < stf2expo) then
        if Sign(f2.fraction) > 0 then
          res := false;
        elsif Sign(f2.fraction) < 0 then
          res := true;
        else
          res := (Sign(f1.fraction) > 0);
        end if;
      else
        if f1deci = f2deci then
          res := (f1.fraction > f2.fraction);
        else
          if f1deci < f2deci then
            mulfra := f1.fraction*integer64(rdx);
            for i in 1..(f2deci-f1deci-1) loop
              Mul(mulfra,integer64(rdx));
            end loop;
            -- Mul_Radix(mulfra,f2deci-f1deci-1);
            res := (mulfra > f2.fraction);
          else
            mulfra := f2.fraction*integer64(rdx);
            for i in 1..(f1deci-f2deci-1) loop
              Mul(mulfra,integer64(rdx));
            end loop;
            -- Mul_Radix(mulfra,f1deci-f2deci-1);
            res := (f1.fraction > mulfra);
          end if;
          Clear(mulfra);
        end if;
      end if;
    else 
      mpf1expo := f1.exponent + integer64(f1deci);
      mpf2expo := f2.exponent + integer64(f2deci);
      if (mpf1expo > mpf2expo) then
        if Sign(f1.fraction) > 0 then
          res := true;
        elsif Sign(f1.fraction) < 0 then
          res := false;
        else
          res := (Sign(f2.fraction) < 0);
        end if;
      elsif (mpf1expo < mpf2expo) then
        if Sign(f2.fraction) > 0 then
          res := false;
        elsif Sign(f2.fraction) < 0 then
          res := true;
        else
          res := (Sign(f1.fraction) > 0);
        end if;
      else
        if f1deci = f2deci then
          res := (f1.fraction > f2.fraction);
        else 
          if f1deci < f2deci then
            mulfra := f1.fraction*integer64(rdx);
            for i in 1..(f2deci-f1deci-1) loop
              Mul(mulfra,integer64(rdx));
            end loop;
            -- Mul_Radix(mulfra,f2deci-f1deci-1);
            res := (mulfra > f2.fraction);
          else
            mulfra := f2.fraction*integer64(rdx);
            for i in 1..(f1deci-f2deci-1) loop
              Mul(mulfra,integer64(rdx));
            end loop;
            -- Mul_Radix(mulfra,f1deci-f2deci-1);
            res := (f1.fraction > mulfra);
          end if;
          Clear(mulfra);
        end if;
      end if;
      Clear(mpf1expo); Clear(mpf2expo);
    end if;
    return res;
  end ">";

  procedure Copy ( f1 : in Floating_Number; f2 : in out Floating_Number ) is
  begin
    Clear(f2);
    Copy(f1.fraction,f2.fraction);
    Copy(f1.exponent,f2.exponent);
  end Copy;

-- EXPANDING AND SHORTENING THE MANTISSA :

  function Expand ( i : Integer_Number; k : natural32 )
                  return Array_of_Naturals is

  -- DESCRIPTION :
  --   Returns the expanded coefficient vector of i.

    cff : constant Array_of_Naturals := Coefficients(Unsigned(i));
    res : Array_of_Naturals(cff'first..cff'last+k);

  begin
    res(cff'range) := cff(cff'range);
    res(cff'last+1..res'last) := (cff'last+1..res'last => 0);
    return res;
  end Expand;

  procedure Truncate ( i : in Integer_Number; k : in natural32;
                       trc : out Array_of_Naturals; exp : out integer32 ) is

  -- DESCRIPTION :
  --   On return is the truncated coefficient vector of i and the number
  --   of decimal places that are lost, or gained when exp is positive.

    icff : Array_of_Naturals(0..Size(i)) := Coefficients(i);
    ns : natural32;          -- #shifts = #decimal places gained

  begin
    Shift_Left(icff,ns);
    exp := integer32(k*the_expo) - integer32(ns);
    for j in trc'range loop
      trc(j) := icff(j+k);
    end loop;
  end Truncate;

 -- procedure Truncate ( i : in out Integer_Number; k : in natural;
 --                      exp : out integer ) is
 --
 --   neg : constant boolean := Negative(i);
 --   ns : natural;
 --   newint : Integer_Number;
 --   icff : Array_of_Naturals(0..Size(i)) := Coefficients(i);
 --
 -- begin
 --   Shift_Left(icff,ns);
 --   exp := k*the_expo - ns;
 --   declare
 --     trc : Array_of_Naturals(icff'first..icff'last-k);
 --   begin
 --     for j in trc'range loop
 --       trc(j) := icff(j+k);
 --     end loop;
 --     newint := Create(trc);
 --     Clear(i);
 --     i := newint;
 --     if neg
 --      then Min(i);
 --     end if;
 --   end;
 -- end Truncate;

  procedure Round ( i : in Integer_Number; k : in natural32;
                    trc : out Array_of_Naturals; exp : out integer32 ) is

  -- DESCRIPTION :
  --   On return is the rounded coefficient vector of i and the number
  --   of decimal places that are lost, or gained when exp is positive.

    icff : Array_of_Naturals(0..Size(i)) := Coefficients(i);
    ns : natural32;          -- #shifts = #decimal places gained

  begin
    Shift_Left(icff,ns);
    exp := integer32(k*the_expo) - integer32(ns);
    for j in trc'range loop
      trc(j) := icff(j+k);
    end loop;
    if k > 0 then
      if icff(k-1) >= the_base/2
       then trc(0) := trc(0)+1;
      end if;
    end if;
  end Round;

  function Expand ( f : Floating_Number; k : natural32 )
                  return Floating_Number is

    res : Floating_Number;
    cff : constant Array_of_Naturals := Expand(f.fraction,k);

  begin
    res.fraction := Create(Create(cff));
    if Negative(f.fraction)
     then Min(res.fraction);
    end if;
    Copy(f.exponent,res.exponent);
    Normalize(res);
    return res;
  end Expand;

  procedure Expand ( f : in out Floating_Number; k : in natural32 ) is

    cff : constant Array_of_Naturals := Expand(f.fraction,k);
    neg : constant boolean := Negative(f.fraction);

  begin
    Clear(f.fraction);
    f.fraction := Create(cff);
    if neg
     then Min(f.fraction);
    end if;
    Normalize(f);
  end Expand;

  function Round ( f : Floating_Number ) return double_float is

    res : double_float;
    szf,ns : natural32;
    n : Natural_Number;
    expo : Integer_Number;

  begin
    if (Empty(f.fraction) or else Equal(f.fraction,0)) then
      res := 0.0;
    else
      Copy(Unsigned(f.fraction),n);
      Shift_Left(n,ns);
      expo := f.exponent - integer64(ns);
      szf := Size(n);
      while Coefficient(n,szf) = 0 and szf > 0 loop
        szf := szf - 1;
      end loop;
      res := double_float(Coefficient(n,szf));
      if szf >= 1 then
        res := res*double_float(the_base)
                 + double_float(Coefficient(n,szf-1));
        if szf >= 2 then
          if Coefficient(n,szf-2) > the_base/2
           then res := res+1.0;
          end if;
          ns := (szf - 1)*the_expo;
          Add(expo,integer64(ns));
        end if;
      end if;
      Adjust_Expo(res,expo);
      Clear(n); Clear(expo);
      if Negative(f.fraction)
       then res := -res;
      end if;
    end if;
    return res;
  end Round;

  function Round ( f : Floating_Number; k : natural32 )
                 return Floating_Number is

    res : Floating_Number;
    cff : Array_of_Naturals(0..(Size(f.fraction)-k));
    exp : integer32;

  begin
    Round(f.fraction,k,cff,exp);
    res.fraction := Create(Create(cff));
    if Negative(f.fraction)
     then Min(res.fraction);
    end if;
    Copy(f.exponent,res.exponent);
    Add(res.exponent,integer64(exp));
    Normalize(res);
    return res;
  end Round;

  procedure Round ( f : in out Floating_Number; k : in natural32 ) is

    exp : integer32;
    neg : constant boolean := Negative(f.fraction);
    ns : natural32;
    icff : Array_of_Naturals(0..Size(f.fraction)) := Coefficients(f.fraction);

  begin
    Shift_Left(icff,ns);
    exp := integer32(k*the_expo) - integer32(ns);
    declare
      trc : Array_of_Naturals(icff'first..icff'last-k);
    begin
      for i in trc'range loop
        trc(i) := icff(i+k);
      end loop;
      if k > 0
       then if icff(k-1) >= the_base/2
             then trc(0) := trc(0)+1;
            end if;
      end if;
      Clear(f.fraction); 
      f.fraction := Create(trc);
      if neg
       then Min(f.fraction);
      end if;
    end;
    Add(f.exponent,integer64(exp));
    Normalize(f);
  end Round;

  function Trunc ( f : Floating_Number ) return double_float is

    res : double_float;
    szf,ns : natural32;
    n : Natural_Number;
    expo : Integer_Number;

  begin
    if (Empty(f.fraction) or else Equal(f.fraction,0)) then
      res := 0.0;
    else
      Copy(Unsigned(f.fraction),n);
      Shift_Left(n,ns);
      expo := f.exponent - integer64(ns);
      szf := Size(n);
      while Coefficient(n,szf) = 0 and szf > 0 loop
        szf := szf - 1;
      end loop;
      res := double_float(Coefficient(f.fraction,szf));
      if szf >= 1 then
        res := res*double_float(the_base)
                 + double_float(Coefficient(f.fraction,szf-1));
        if szf >= 2 then
          ns := (szf - 1)*the_expo;
          Add(expo,integer64(ns));
        end if;
      end if;
      Adjust_Expo(res,expo);
      Clear(n); Clear(expo);
      if Negative(f.fraction)
       then res := -res;
      end if;
    end if;
    return res;
  end Trunc;

  function Trunc ( f : Floating_Number; k : natural32 )
                 return Floating_Number is
 
    res : Floating_Number;
    cff : Array_of_Naturals(0..(Size(f.fraction)-k));
    exp : integer32;

  begin
    Truncate(f.fraction,k,cff,exp);
    res.fraction := Create(Create(cff));
    if Negative(f.fraction)
     then Min(res.fraction);
    end if;
    Copy(f.exponent,res.exponent);
    Add(res.exponent,integer64(exp));
    Normalize(res);
    return res;
  end Trunc;

  procedure Trunc ( f : in out Floating_Number; k : in natural32 ) is

    exp : integer32;
    neg : constant boolean := Negative(f.fraction);
    ns : natural32;

  begin
    Shift_Left(f.fraction,ns);
    exp := integer32(k*the_expo) - integer32(ns);
    declare
      cff : constant Array_of_Naturals := Coefficients(f.fraction);
      trc : Array_of_Naturals(cff'first..cff'last-k);
    begin
      for i in trc'range loop
        trc(i) := cff(i+k);
      end loop;
      Clear(f.fraction);
      f.fraction := Create(trc);
      if neg
       then Min(f.fraction);
      end if;
    end;
    Add(f.exponent,integer64(exp));
  end Trunc;

  procedure Set_Size ( f : in out Floating_Number; k : in natural32 ) is

    sz : constant natural32 := Size_Fraction(f);

  begin
    if sz > k then
      Round(f,sz-k);
    elsif sz < k then
      Expand(f,k-sz);
    end if;
  end Set_Size;

-- ARITHMETIC OPERATIONS as functions (no data sharing) :

  function "+" ( f1 : Floating_Number; f2 : double_float )
               return Floating_Number is

    ff2 : Floating_Number := Create(f2);
    res : constant Floating_Number := f1+ff2;

  begin
    Clear(ff2);
    return res;
  end "+";

  function "+" ( f1 : double_float; f2 : Floating_Number )
               return Floating_Number is
  begin
    return f2+f1;
  end "+";

  function "+" ( f1,f2 : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    maxsz,max_diff : natural32;
    diff_size : integer32;
    diff_expo : Integer_Number;
    cnt : integer64;

  begin
    if (Empty(f1.fraction) or else Equal(f1.fraction,0)) then
      Copy(f2,res);
    elsif (Empty(f2.fraction) or else Equal(f2.fraction,0)) then
      Copy(f1,res);
    else
      maxsz := Max_Size(f1.fraction,f2.fraction) + 1;
      if Equal(f1.exponent,f2.exponent) then
        res.fraction := f1.fraction + f2.fraction;
        Copy(f1.exponent,res.exponent);
      else
        diff_expo := f1.exponent - f2.exponent;
        max_diff := maxsz*the_expo;
        if diff_expo > 0 then
          Copy(f1.fraction,res.fraction);
          if diff_expo < integer64(2*max_diff) then
            cnt := Create(diff_expo);
            for i in 1..cnt loop
              Mul(res.fraction,integer64(rdx));
            end loop;
           -- Mul_Radix(res.fraction,cnt);
            Add(res.fraction,f2.fraction);
            Copy(f2.exponent,res.exponent);
          else
            Copy(f1.exponent,res.exponent);
          end if;
        else
          Copy(f2.fraction,res.fraction);
          Min(diff_expo);
          if diff_expo < integer64(2*max_diff) then
            cnt := Create(diff_expo);
            for i in 1..cnt loop
              Mul(res.fraction,integer64(rdx));
            end loop;
           -- Mul_Radix(res.fraction,cnt);
            Add(res.fraction,f1.fraction);
            Copy(f1.exponent,res.exponent);
          else
            Copy(f2.exponent,res.exponent);
          end if;
        end if;
        Clear(diff_expo);
      end if;
      diff_size := integer32(Size(res.fraction)) + 1 - integer32(maxsz);
      if diff_size > 0
       then Round(res,natural32(diff_size));
       else Normalize(res);
      end if;
     -- if Size_Fraction(res) /= maxsz - 1
     --  then put_line("in + result not of max size");
     -- end if;
    end if;
    return res;
  end "+";

  function "+" ( f : Floating_Number ) return Floating_Number is

    res : Floating_Number;

  begin
    Copy(f,res);
    return res;
  end "+";

  function "-" ( f : Floating_Number ) return Floating_Number is

    res : Floating_Number;

  begin
    Copy(f,res);
    Min(res.fraction);
    return res;
  end "-";

  function "-" ( f1 : Floating_Number; f2 : double_float )
               return Floating_Number is

    minf2 : constant double_float := -f2;
    res : constant Floating_Number := f1+minf2;

  begin
    return res;
  end "-";

  function "-" ( f1 : double_float; f2 : Floating_Number )
               return Floating_Number is

    ff1 : Floating_Number := Create(f1);
    res : constant Floating_Number := ff1 - f2;

  begin
    Clear(ff1);
    return res;
  end "-";

  function "-" ( f1,f2 : Floating_Number ) return Floating_Number is

    res,minf2 : Floating_Number;

  begin
    minf2.exponent := f2.exponent;
    minf2.fraction := -f2.fraction;
    res := f1 + minf2;
    Clear(minf2.fraction);
    return res;
  end "-";

  function "*" ( f1 : Floating_Number; f2 : double_float )
               return Floating_Number is

    ff2 : Floating_Number := Create(f2);
    res : constant Floating_Number := f1*ff2;

  begin
    Clear(ff2);
    return res;
  end "*";

  function "*" ( f1 : double_float; f2 : Floating_Number )
               return Floating_Number is

    ff1 : Floating_Number := Create(f1);
    res : constant Floating_Number := ff1*f2;

  begin
    Clear(ff1);
    return res;
  end "*";

  function "*" ( f1,f2 : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    diffsize : integer32;

  begin
    if not (Empty(f1.fraction) or else Equal(f1.fraction,0)) then
      if not (Empty(f2.fraction) or else Equal(f2.fraction,0)) then
        res.fraction := f1.fraction*f2.fraction;
        res.exponent := f1.exponent+f2.exponent;
        diffsize := integer32(Size(res.fraction))
                     - integer32(Max_Size(f1.fraction,f2.fraction));
        if diffsize > 0
         then Round(res,natural32(diffsize));
        end if;
        Normalize(res);
      else
        res := Create(integer64(0));
      end if;
    else
      res := Create(integer64(0));
    end if;
    return res;
  end "*";

  function "**" ( f : double_float; n : Natural_Number )
                return Floating_Number is

    ff : Floating_Number := Create(f);
    res : constant Floating_Number := ff**n;

  begin
    Clear(ff);
    return res;
  end "**";

  function "**" ( f : Floating_Number; n : Natural_Number )
                return Floating_Number is

    res : Floating_Number;
    cnt : Natural_Number;

  begin
    if (Empty(n) or else Equal(n,0)) then
      res := Create(integer64(1));
    else
      Copy(f,res);
      cnt := Create(1);
      while cnt < n loop
        Mul(res,f);
        Add(cnt,1);
      end loop;
      Clear(cnt);
      Normalize(res);
    end if;
    return res;
  end "**";

  function "**" ( f : Floating_Number; i : integer32 ) return Floating_Number is

    res : Floating_Number;

  begin
    if i = 0 then
      res := Create(1.0);
    else
      if i > 0 then
        Copy(f,res);
        for j in 1..(i-1) loop
          Mul(res,f);
        end loop;
      else
        res := Create(integer64(1));
        for j in 1..(-i) loop
          Div(res,f);
        end loop;
      end if;
      Normalize(res);
    end if;
    return res;
  end "**";

  function "**" ( f : double_float; i : Integer_Number )
                return Floating_Number is

    ff : Floating_Number := Create(f);
    res : constant Floating_Number := ff**i;

  begin
    Clear(ff);
    return res;
  end "**";

  function "**" ( f : Floating_Number; i : Integer_Number )
                return Floating_Number is

    res : Floating_Number;
    cnt,n : Natural_Number;

  begin
    if (Empty(i) or else Equal(i,0)) then
      res := Create(integer64(1));
    else
      n := Unsigned(i);
      if i > 0 then
        Copy(f,res);
        cnt := Create(1);
        while cnt < n loop
          Mul(res,f);
          Add(cnt,1);
        end loop;
      else
        res := Create(integer64(1));
        cnt := Create(0);
        while cnt < n loop
          Div(res,f);
          Add(cnt,1);
        end loop;
      end if;
      Clear(cnt);
      Normalize(res);
    end if;
    return res;
  end "**";

  function "/" ( f1 : Floating_Number; f2 : double_float )
               return Floating_Number is

    ff2 : Floating_Number := Create(f2);
    res : constant Floating_Number := f1/ff2;

  begin
    Clear(ff2);
    return res;
  end "/";

  function "/" ( f1 : double_float; f2 : Floating_Number )
               return Floating_Number is

    ff1 : Floating_Number := Create(f1);
    res : constant Floating_Number := ff1/f2;

  begin
    Clear(ff1);
    return res;
  end "/";

  function Pos_Div ( f1,f2 : Floating_Number ) return Floating_Number is

  -- DESCRIPTION :
  --   This division only applies when f1 and f2 are both positive.

    res : Floating_Number;
    diffsize : integer32;
    maxsz : constant natural32 := Max_Size(f1.fraction,f2.fraction);
    f1frac : constant Array_of_Naturals(0..Size(f1.fraction))
           := Coefficients(f1.fraction);
    f2frac : constant Array_of_Naturals(0..Size(f2.fraction))
           := Coefficients(f2.fraction);
    resfrac,backup : Array_of_Naturals(0..maxsz+1);
    rmdfrac,new_rmdfrac : Array_of_Naturals(0..f2frac'last+1);
    f2fraclast,sizeresfrac,cnt : natural32;
    iszero : boolean;

  begin
   -- put_line("in function Pos_Div :");
   -- put("f1.frac : "); put(f1.fraction); new_line;
   -- put("f2.frac : "); put(f2.fraction); new_line;
    res.exponent := f1.exponent - f2.exponent;
    f2fraclast := 0;
    for i in reverse f2frac'range loop
      if f2frac(i) /= 0
       then f2fraclast := i; exit;
      end if;
    end loop;
    resfrac := (resfrac'range => 0);
    if f2frac > f1frac then
      sizeresfrac := 0;
      for i in rmdfrac'range loop
        if i <= f1frac'last
         then rmdfrac(i) := f1frac(i);
         else rmdfrac(i) := 0;
        end if;
      end loop;
    else
      if f2fraclast = 0 then
        if f2frac(0) <= rdx
         then Small_Div(f1frac,f2frac(0),resfrac,rmdfrac(0));
         else Big_Div(f1frac,f2frac(0),resfrac,rmdfrac(0));
        end if;
        for i in 1..rmdfrac'last loop
          rmdfrac(i) := 0;
        end loop;
      else
        declare
          rest : Array_of_Naturals(0..f2fraclast);
        begin
         -- put_line("Calling Div from function Pos_Div");
         -- put("f1.fraction : "); put(f1.fraction); new_line;
         -- put("f1.exponent : "); put(f1.exponent); new_line;
         -- put("f2.fraction : "); put(f2.fraction); new_line;
         -- put("f2.exponent : "); put(f2.exponent); new_line;
          Div(f1frac,f2frac(0..f2fraclast),resfrac,rest);
          rmdfrac(rest'range) := rest;
          for i in f2fraclast+1..rmdfrac'last loop
            rmdfrac(i) := 0;
          end loop;
        end;
      end if;
      sizeresfrac := 0;
      for i in reverse resfrac'range loop
        if resfrac(i) /= 0
         then sizeresfrac := i; exit;
        end if;
      end loop;
    end if;
    diffsize := integer32(sizeresfrac) - integer32(maxsz);
    iszero := true;
    for i in rmdfrac'range loop
      if rmdfrac(i) /= 0
       then iszero := false; exit;
      end if;
    end loop;
    while not iszero and (diffsize <= 0) loop
      cnt := 0;
     -- put_line("resfrac before the shift : ");
     -- for i in resfrac'range loop
     --   put(resfrac(i));
     -- end loop;
     -- new_line;
      backup := resfrac;
      while rmdfrac < f2frac loop
       -- Mul_Fact(rmdfrac,rdx);
       -- Mul_Fact(resfrac,rdx);
        if resfrac(resfrac'last) >= sub_base then
          iszero := true;
          -- put_line("we must leave the loop with resfrac :");
          -- for i in resfrac'range loop
          --   put(resfrac(i));
          -- end loop;
          -- new_line;
          resfrac := backup;
        end if;  -- but here rmdfrac >= f2frac ??????????
        exit when iszero;
        Mul_Fact(rmdfrac,rdx);
        Mul_Fact(resfrac,rdx);
        cnt := cnt+1;
      end loop;
      exit when iszero;
     -- put("resfrac after the shift with cnt = "); put(cnt,1);
     -- put_line(" : ");
     -- for i in resfrac'range loop
     --   put(resfrac(i));
     -- end loop;
     -- new_line;
      Sub(res.exponent,integer64(cnt));
      new_rmdfrac := (new_rmdfrac'range => 0);
      if f2fraclast = 0 then
        if f2frac(0) <= rdx
         then Small_Div(rmdfrac,f2frac(0),new_rmdfrac(0));
         else Big_Div(rmdfrac,f2frac(0),new_rmdfrac(0));
        end if;
      else
        declare
          rest : Array_of_Naturals(0..f2fraclast);
        begin
         -- put_line("Calling Div from function Pos_Div");
         -- put("f1.fraction : "); put(f1.fraction); new_line;
         -- put("f1.exponent : "); put(f1.exponent); new_line;
         -- put("f2.fraction : "); put(f2.fraction); new_line;
         -- put("f2.exponent : "); put(f2.exponent); new_line;
          Div(rmdfrac,f2frac(0..f2fraclast),rest);
          new_rmdfrac(rest'range) := rest;
          for i in f2fraclast+1..new_rmdfrac'last loop
            new_rmdfrac(i) := 0;
          end loop;
        end;
      end if;
      Add(resfrac,rmdfrac);
      sizeresfrac := 0;
      for i in reverse resfrac'range loop
        if resfrac(i) /= 0
         then sizeresfrac := i; exit;
        end if;
      end loop;
      diffsize := integer32(sizeresfrac) - integer32(maxsz);
      rmdfrac := new_rmdfrac;
      iszero := true;
      for i in rmdfrac'range loop
        if rmdfrac(i) /= 0
         then iszero := false; exit;
        end if;
      end loop;
    end loop;
    if sizeresfrac <= maxsz
     then res.fraction := Create(resfrac(0..maxsz));
     else res.fraction := Create(resfrac);
    end if;
    if diffsize > 0
     then Round(res,natural32(diffsize));
    end if;
    Normalize(res);
    return res;
  end Pos_Div;

  function "/" ( f1,f2 : Floating_Number ) return Floating_Number is

    res,minf1,minf2 : Floating_Number;

  begin
    if not (Empty(f1.fraction) or else Equal(f1.fraction,0)) then
      if not (Empty(f2.fraction) or else Equal(f2.fraction,0)) then
        if Multprec_Integer64_Numbers.Positive(f1.fraction) then
          if Multprec_Integer64_Numbers.Positive(f2.fraction) then
            res := Pos_Div(f1,f2);
          else
            minf2.fraction := -f2.fraction;
            minf2.exponent := f2.exponent;
            res := Pos_Div(f1,minf2);
            Clear(minf2.fraction);
            Min(res);
          end if;
        else
          minf1.fraction := -f1.fraction;
          minf1.exponent := f1.exponent;
          if Multprec_Integer64_Numbers.Positive(f2.fraction) then
            res := Pos_Div(minf1,f2);
            Min(res);
          else
            minf2.fraction := -f2.fraction;
            minf2.exponent := f2.exponent;
            res := Pos_Div(minf1,minf2);
            Clear(minf2.fraction);
          end if;
          Clear(minf1.fraction);
        end if;
      else
        raise NUMERIC_ERROR;
      end if;
    else
      res := Create(integer64(0));
    end if;
    return res;
  end "/";

-- ARITHMETIC OPERATIONS as procedures for memory management :

  procedure Add ( f1 : in out Floating_Number; f2 : in double_float ) is

    ff2 : Floating_Number := Create(f2);

  begin
    Add(f1,ff2);
    Clear(ff2);
  end Add;

  procedure Add ( f1 : in out Floating_Number; f2 : in Floating_Number ) is

    res : Floating_Number;
    maxsz,max_diff : natural32;
    diff_size : integer32;
    diff_expo : Integer_Number;
   -- szf1,szf2 : natural32;
    cnt : integer64;

  begin
    if (Empty(f1.fraction) or else Equal(f1.fraction,0)) then
      Copy(f1 => f2,f2 => f1);
    elsif (Empty(f2.fraction) or else Equal(f2.fraction,0)) then
      null;
    else
      maxsz := Max_Size(f1.fraction,f2.fraction) + 1;
     -- szf1 := Size(f1.fraction);
     -- szf2 := Size(f2.fraction);
      if Equal(f1.exponent,f2.exponent) then
        Add(f1.fraction,f2.fraction);
      else
        diff_expo := f1.exponent - f2.exponent;
        max_diff := maxsz*the_expo;
        if diff_expo > 0 then
          if diff_expo < integer64(2*max_diff) then
            cnt := Create(diff_expo);
            for i in 1..cnt loop
              Mul(f1.fraction,integer64(rdx));
            end loop;
           -- Mul_Radix(f1.fraction,cnt);
            Add(f1.fraction,f2.fraction);
            Copy(f2.exponent,f1.exponent);
          else
            null;
          end if;
        else
          Copy(f2.fraction,res.fraction);
          Min(diff_expo);
          if diff_expo < integer64(2*max_diff) then
            cnt := Create(diff_expo);
            for i in 1..cnt loop
              Mul(res.fraction,integer64(rdx));
            end loop;
           -- Mul_Radix(res.fraction,cnt);
            Add(f1.fraction,res.fraction);
            Clear(res.fraction);
          else
            Copy(f2.exponent,f1.exponent);
            Clear(f1.fraction);
            f1.fraction := res.fraction;
          end if;
        end if;
        Clear(diff_expo);
      end if;
      diff_size := integer32(Size(f1.fraction)) + 1 - integer32(maxsz);
      if diff_size > 0
       then Round(f1,natural32(diff_size));
       else Normalize(f1);
      end if;
     -- if Size_Fraction(f1) /= maxsz - 1
     --  then put_line("in Add f1, has not size = max size!!!");
     --       put("Size f1 : "); put(Size_Fraction(f1),1);
     --       put("  diff_size : "); put(diff_size,1); new_line;
     --       put("f1.fraction : "); put(f1.fraction); new_line;
     --       put("Original Size f1 : "); put(szf1,1);
     --       put("  Original Size f2 : "); put(szf2,1); new_line;
     -- end if;
    end if;
  end Add;

  procedure Sub ( f1 : in out Floating_Number; f2 : in double_float ) is

    minf2 : constant double_float := -f2;

  begin
    Add(f1,minf2);
  end Sub;

  procedure Sub ( f1 : in out Floating_Number; f2 : in Floating_Number ) is

    minf2 : Floating_Number;

  begin
    if (not Empty(f2.fraction) and then not Equal(f2.fraction,0))
     then minf2.exponent := f2.exponent;
          minf2.fraction := -f2.fraction;
          Add(f1,minf2);
          Clear(minf2.fraction);
    end if;
  end Sub;

  procedure Min ( f : in out Floating_Number ) is
  begin
    Min(f.fraction);
  end Min;

  procedure Mul ( f1 : in out Floating_Number; f2 : in double_float ) is

    ff2 : Floating_Number := Create(f2);

  begin
    Mul(f1,ff2);
    Clear(ff2);
  end Mul;

  procedure Mul ( f1 : in out Floating_Number; f2 : in Floating_Number ) is

    diffsize : integer32;
    maxsize : natural32;

  begin
    if Empty(f2.fraction) then
      Clear(f1); f1 := Create(integer64(0));
    elsif Equal(f2.fraction,0) then
      Clear(f1); f1 := Create(integer64(0));
    elsif Empty(f1.fraction) then
      null;
    elsif Equal(f1.fraction,0) then
      null;
    else
      maxsize := Max_Size(f1.fraction,f2.fraction);
      Mul(f1.fraction,f2.fraction);
      Add(f1.exponent,f2.exponent);
      diffsize := integer32(Size(f1.fraction)) - integer32(maxsize);
      if diffsize > 0
       then Round(f1,natural32(diffsize));
      end if;
      Normalize(f1);
    end if;
  end Mul;

  procedure Div ( f1 : in out Floating_Number; f2 : in double_float ) is

    ff2 : Floating_Number := Create(f2);

  begin
    Div(f1,ff2);
    Clear(ff2);
  end Div;

  procedure Pos_Div ( f1 : in out Floating_Number;
                      f2 : in Floating_Number ) is

  -- DESCRIPTION :
  --   This division only applies when f1 and f2 are both positive.

    diffsize : integer32;
    maxsz : constant natural32 := Max_Size(f1.fraction,f2.fraction);
    f1frac : constant Array_of_Naturals(0..Size(f1.fraction))
           := Coefficients(f1.fraction);
    f2frac : constant Array_of_Naturals(0..Size(f2.fraction))
           := Coefficients(f2.fraction);
    resfrac,backup : Array_of_Naturals(0..maxsz+1);
    rmdfrac,new_rmdfrac : Array_of_Naturals(0..f2frac'last+1);
    f2fraclast,sizeresfrac,cnt : natural32;
    iszero : boolean;

  begin
   -- put_line("in procedure Pos_Div 2 : ");
   -- put("f1.frac : "); put(f1.fraction); new_line;
   -- put("f2.frac : "); put(f2.fraction); new_line;
    Sub(f1.exponent,f2.exponent);
    f2fraclast := 0;
    for i in reverse f2frac'range loop
      if f2frac(i) /= 0
       then f2fraclast := i; exit;
      end if;
    end loop;
    resfrac := (resfrac'range => 0);
    if f2frac > f1frac then
      sizeresfrac := 0;
      for i in rmdfrac'range loop
        if i <= f1frac'last
         then rmdfrac(i) := f1frac(i);
         else rmdfrac(i) := 0;
        end if;
      end loop;
    else
      if f2fraclast = 0 then
        if f2frac(0) <= rdx
         then Small_Div(f1frac,f2frac(0),resfrac,rmdfrac(0));
         else Big_Div(f1frac,f2frac(0),resfrac,rmdfrac(0));
        end if;
        for i in 1..rmdfrac'last loop
          rmdfrac(i) := 0;
        end loop;
      else
        declare
          rest : Array_of_Naturals(0..f2fraclast);
        begin
        -- put_line("Calling Div from procedure Pos_Div");
          Div(f1frac,f2frac(0..f2fraclast),resfrac,rest);
          rmdfrac(rest'range) := rest;
          for i in f2fraclast+1..rmdfrac'last loop
            rmdfrac(i) := 0;
          end loop;
        end;
      end if;
      sizeresfrac := 0;
      for i in reverse resfrac'range loop
        if resfrac(i) /= 0
         then sizeresfrac := i; exit;
        end if;
      end loop;
    end if;
    diffsize := integer32(sizeresfrac) - integer32(maxsz);
    iszero := true;
    for i in rmdfrac'range loop
      if rmdfrac(i) /= 0
       then iszero := false; exit;
      end if;
    end loop;
    while not iszero and (diffsize <= 0) loop
      cnt := 0;
      backup := resfrac;
      while rmdfrac < f2frac loop
        if resfrac(resfrac'last) >= sub_base
         then iszero := true;
              resfrac := backup;
        end if;
        exit when iszero;
        Mul_Fact(rmdfrac,rdx);
        Mul_Fact(resfrac,rdx);
        cnt := cnt+1;
      end loop;
      exit when iszero;
      Sub(f1.exponent,integer64(cnt));
      new_rmdfrac := (new_rmdfrac'range => 0);
      if f2fraclast = 0 then
        if f2frac(0) <= rdx
         then Small_Div(rmdfrac,f2frac(0),new_rmdfrac(0));
         else Big_Div(rmdfrac,f2frac(0),new_rmdfrac(0));
        end if;
      else
        declare
          rest : Array_of_Naturals(0..f2fraclast);
        begin
         -- put_line("Calling Div from procedure Pos_Div");
          Div(rmdfrac,f2frac(0..f2fraclast),rest);
          new_rmdfrac(rest'range) := rest;
          for i in f2fraclast+1..new_rmdfrac'last loop
            new_rmdfrac(i) := 0;
          end loop;
        end;
      end if;
      Add(resfrac,rmdfrac);
      sizeresfrac := 0;
      for i in reverse resfrac'range loop
        if resfrac(i) /= 0
         then sizeresfrac := i; exit;
        end if;
      end loop;
      diffsize := integer32(sizeresfrac) - integer32(maxsz);
      rmdfrac := new_rmdfrac;
      iszero := true;
      for i in rmdfrac'range loop
        if rmdfrac(i) /= 0
         then iszero := false; exit;
        end if;
      end loop;
    end loop;
    Clear(f1.fraction);
    if sizeresfrac <= maxsz
     then f1.fraction := Create(resfrac(0..maxsz));
     else f1.fraction := Create(resfrac);
    end if;
    if diffsize > 0
     then Round(f1,natural32(diffsize));
    end if;
    Normalize(f1);
  end Pos_Div;

  procedure Div ( f1 : in out Floating_Number; f2 : in Floating_Number ) is

    minf1,minf2 : Floating_Number;

  begin
    if not (Empty(f1.fraction) or else Equal(f1.fraction,0)) then
      if not (Empty(f2.fraction) or else Equal(f2.fraction,0)) then
        if Multprec_Integer64_Numbers.Positive(f1.fraction) then
          if Multprec_Integer64_Numbers.Positive(f2.fraction) then
            Pos_Div(f1,f2);
          else
            minf2.fraction := -f2.fraction;
            minf2.exponent := f2.exponent;
            Pos_Div(f1,minf2);
            Clear(minf2.fraction); Min(f1);
          end if;
        else
          minf1 := -f1;
          if Multprec_Integer64_Numbers.Positive(f2.fraction) then
            Pos_Div(minf1,f2);
            Clear(f1); f1 := minf1; Min(f1);
          else
            minf2.fraction := -f2.fraction;
            minf2.exponent := f2.exponent;
            Pos_Div(minf1,minf2);
            Clear(minf2.fraction);
            Clear(f1); f1 := minf1;
          end if;
        end if;
      else 
        raise NUMERIC_ERROR;
      end if;
    end if;
  end Div;

-- DESTRUCTOR :

  procedure Clear ( f : in out Floating_Number ) is
  begin
    Clear(f.fraction);
    Clear(f.exponent);
  end Clear;

end Multprec_Floating64_Numbers;
