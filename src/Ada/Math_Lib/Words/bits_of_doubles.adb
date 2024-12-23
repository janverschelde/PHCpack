with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Mask_Bits_of_Doubles;

package body Bits_of_Doubles is

  procedure expand_52bits
              ( bits : out Standard_Natural_Vectors.Vector;
                nbr : in integer64 ) is

    temp : integer64 := nbr;
    idx : integer32;

  begin
    for i in 0..51 loop
      idx := integer32(51 - i);
      bits(idx) := natural32(temp mod 2);
      temp := temp/2;
    end loop;
  end expand_52bits;

  procedure write_52bits
              ( bits : in Standard_Natural_Vectors.Vector ) is

  begin
    for i in 0..51 loop
      if i > 0 then
        if i mod 4 = 0
         then put(" ");
        end if;
      end if;
      put(integer32(bits(integer32(i))));
    end loop;
  end write_52bits;

  procedure Fraction_Exponent
              ( x : in double_float;
                f : out integer64; e : out integer32 ) is

    r : constant double_float := double_float'fraction(x);
    s : constant double_float := double_float'compose(r, 52);

  begin
    e := integer32(double_float'exponent(x));
    f := integer64(double_float'truncation(s));
  end Fraction_Exponent;

  procedure write_52bits_expo ( x : in double_float ) is

    e : integer32;
    m : integer64;
    b : Standard_Natural_Vectors.Vector(0..51);

  begin
    Fraction_Exponent(x,m,e);
    expand_52bits(b,m);
    if x < 0.0
     then put("-");
     else put("+");
    end if;
    write_52bits(b);
    put(", "); put(e,1);
  end write_52bits_expo;

  function Bit_Equal ( x,y : double_float ) return boolean is

    ex,ey : integer32;
    mx,my : integer64;

  begin
    if x < 0.0 and y > 0.0 then
      return false;
    elsif x > 0.0 and y < 0.0 then
      return false;
    else
      Fraction_Exponent(x,mx,ex);
      Fraction_Exponent(y,my,ey);
      if ex /= ey then
        return false;
      else
        return mx = my;
      end if;
    end if;
  end Bit_Equal;

  procedure write_fraction_bits ( nbr : in double_float ) is

    e : integer32;
    mfr : integer64;

  begin
    Fraction_Exponent(nbr,mfr,e);
    put(mfr,1,b=>2); new_line;
  end write_fraction_bits;

  function value_52bits
             ( bits : Standard_Natural_Vectors.Vector ) return integer64 is

    result : integer64 := 0;
 
  begin
    for i in 0..51 loop
      result := 2*result + integer64(bits(integer32(i)));
    end loop;
    return result;
  end value_52bits;

  function chop_last_bits
             ( nbr : double_float; lastbits : natural32 )
             return double_float is

    frc : constant double_float := double_float'fraction(nbr);
    exn : constant integer32 := integer32(double_float'exponent(nbr));
    sfr : constant double_float := double_float'compose(frc, 52);
    mfr : constant integer64 := integer64(double_float'truncation(sfr));
    bits : Standard_Natural_Vectors.Vector(0..51);
    value : integer64;

  begin
    expand_52bits(bits,mfr);
    for i in 0..(lastbits-1) loop
      bits(integer32(51-i)) := 0;
    end loop;
    value := value_52bits(bits);
    return double_float'compose(double_float(value), exn);
  end chop_last_bits;

  procedure chop_last_bits
             ( nbr : in out double_float; lastbits : in natural32;
               headbits : out Standard_Natural_Vectors.Vector;
               tailbits : out Standard_Natural_Vectors.Vector ) is

    frc : constant double_float := double_float'fraction(nbr);
    exn : constant integer32 := integer32(double_float'exponent(nbr));
    sfr : constant double_float := double_float'compose(frc, 52);
    mfr : constant integer64 := integer64(double_float'truncation(sfr));
    value : integer64;

  begin
    expand_52bits(headbits,mfr);
    for i in 0..51-lastbits loop
      tailbits(integer32(i)) := 0;
    end loop;
    for i in 0..(lastbits-1) loop
      tailbits(integer32(51-i)) := headbits(integer32(51-i));
      headbits(integer32(51-i)) := 0;
    end loop;
    value := value_52bits(headbits);
    nbr := double_float'compose(double_float(value), exn);
  end chop_last_bits;

  procedure insert_first_bits
              ( bits : in out Standard_Natural_Vectors.Vector;
                firstbits : in natural32;
                headbits : in Standard_Natural_Vectors.Vector ) is
  begin
    for i in 0..firstbits-1 loop
      bits(integer32(51-i)) := bits(integer32(51-i-firstbits));
      bits(integer32(51-i-firstbits)) := 0;
    end loop;
    for i in 0..firstbits-1 loop
      bits(integer32(i)) := headbits(integer32(i));
    end loop;
  end insert_first_bits;

  procedure insert_first_bits
              ( nbr : in out double_float;
                firstbits : in natural32;
                headbits : in Standard_Natural_Vectors.Vector ) is

    fnbr : constant double_float := double_float'fraction(nbr);
    enbr : constant integer32 := integer32(double_float'exponent(nbr));
    snbr : constant double_float := double_float'compose(fnbr, 52);
    mnbr : constant integer64 := integer64(double_float'truncation(snbr));
    bits : Standard_Natural_Vectors.Vector(0..51);
    value : integer64;
    ifirst : constant integer32 := integer32(firstbits);

  begin
    expand_52bits(bits, mnbr);
    insert_first_bits(bits, firstbits, headbits);
    value := value_52bits(bits);
    nbr := double_float'compose(double_float(value), enbr + ifirst);
  end insert_first_bits;

  function insert_first_bits
             ( nbr : double_float;
               firstbits : natural32;
               headbits : Standard_Natural_Vectors.Vector )
             return double_float is

    res : double_float := nbr;
  
  begin
    insert_first_bits(res,firstbits,headbits);
    return res;
  end insert_first_bits;

  procedure Mod_Split ( x : in double_float;
                        xhi,xlo : out double_float ) is

    e : constant integer32 := integer32(double_float'exponent(x));
    f : constant double_float := double_float'fraction(x);
    s : constant double_float := double_float'compose(f, 52);
    m : constant integer64 := integer64(double_float'truncation(s));
    mlast,mchop : integer64;
    cnt : integer32 := 0;

  begin
    mlast := Mask_Bits_of_Doubles.last_bits(m,26);
    mchop := m - mlast;
    put("m     : "); put(m,1,b=>2); new_line;
    put("mchop : "); put(mchop,1,b=>2); new_line;
    xhi := double_float'compose(double_float(mchop),e);
    if mlast /= 0 then
      while mlast < 2**natural(25-cnt) loop
        cnt := cnt + 1;
      end loop;
      put("mlast : "); put(mlast,1,b=>2); new_line;
      put("cnt : "); put(cnt,1); new_line;
    end if;
    xlo := double_float'compose(double_float(mlast),e-26-cnt);
  end Mod_Split;

  procedure Vec_Split ( x : in double_float;
                        xhi,xlo : out double_float ) is

    head,tail : Standard_Natural_Vectors.Vector(0..51);
    val : integer64;
    expo : constant integer32 := integer32(double_float'exponent(x));
    cnt : integer32 := 0;

  begin
    xhi := x;
    chop_last_bits(xhi,26,head,tail);
    val := value_52bits(tail);
    if val /= 0 then
      while tail(cnt+26) = 0 loop
        cnt := cnt + 1;
      end loop;
    end if;
    xlo := double_float'compose(double_float(val),expo-26-cnt);
  end Vec_Split;

  procedure Split ( x : in double_float;
                    x0,x1,x2,x3 : out double_float ) is

    f : constant double_float := double_float'fraction(x);
    e : constant integer32 := integer32(double_float'exponent(x));
    s : constant double_float := double_float'compose(f, 52);
    m : constant integer64 := integer64(double_float'truncation(s));
    thebits,xp0,xp1,xp2,xp3 : Standard_Natural_Vectors.Vector(0..51);
    part : constant integer32 := 52/4;
    valbits : integer64;
    idx,cnt : integer32;
  
  begin
    expand_52bits(thebits,m);
   -- put("the bits : "); write_52bits(thebits); new_line;
    for i in 0..(part-1) loop
      xp0(i) := thebits(i);
      xp1(i) := 0; xp2(i) := 0; xp3(i) := 0;
    end loop;
    for i in part..(2*part-1) loop
      xp1(i) := thebits(i);
      xp0(i) := 0; xp2(i) := 0; xp3(i) := 0;
    end loop;
    for i in 2*part..(3*part-1) loop
      xp2(i) := thebits(i);
      xp0(i) := 0; xp1(i) := 0; xp3(i) := 0;
    end loop;
    for i in 3*part..51 loop
      xp3(i) := thebits(i);
      xp0(i) := 0; xp1(i) := 0; xp2(i) := 0;
    end loop;
   -- put("1st part : "); write_52bits(xp0); new_line;
   -- put("2nd part : "); write_52bits(xp1); new_line;
   -- put("3rd part : "); write_52bits(xp2); new_line;
   -- put("4th part : "); write_52bits(xp3); new_line;
    valbits := value_52bits(xp0);
    x0 := double_float'compose(double_float(valbits),e);
    if x < 0.0
     then x0 := -x0;
    end if;
    valbits := value_52bits(xp1);
    cnt := 0;
    if valbits /= 0 then
      idx := part;
      while xp1(idx) = 0 loop
        idx := idx + 1;
        cnt := cnt + 1;
      end loop;
    end if;
    x1 := double_float'compose(double_float(valbits),e - part - cnt);
    if x < 0.0
     then x1 := -x1;
    end if;
    valbits := value_52bits(xp2);
    cnt := 0;
    if valbits /= 0 then
      idx := 2*part;
      while xp2(idx) = 0 loop
        idx := idx + 1;
        cnt := cnt + 1;
      end loop;
    end if;
    x2 := double_float'compose(double_float(valbits),e - 2*part - cnt);
    if x < 0.0
     then x2 := -x2;
    end if;
    valbits := value_52bits(xp3);
    cnt := 0;
    if valbits /= 0 then
      idx := 3*part;
      while xp3(idx) = 0 loop
        idx := idx + 1;
        cnt := cnt + 1;
      end loop;
    end if;
   -- x3 := double_float'compose(double_float(valbits),e - 3*part - cnt);
   -- if x < 0.0
   --  then x3 := -x3;
   -- end if;
    x3 := x - (x0 + x1 + x2);
  end Split;

  procedure Sign_Balance ( hi,lo : in out double_float;
                           verbose : in boolean := true ) is

    mhi,mlo : double_float;

  begin
    if hi < 0.0 then
      mhi := -hi;
      mlo := -lo;
      Sign_Balance(mhi,mlo,verbose);
      hi := -mhi;
      lo := -mlo;
    else
      for k in 1..52 loop
        mhi := hi - chop_last_bits(hi,natural32(k));
        if verbose then
          put("b hi : "); write_fraction_bits(hi);
          put("  hi : "); put(hi); new_line;
          put("  lo : "); put(lo); new_line;
          put("last bit : "); put(mhi); new_line;
        end if;
        exit when (mhi > 0.0);
      end loop;
      hi := hi - mhi;
      lo := lo + mhi;
      if verbose then
        put("  hi : "); put(hi); new_line;
        put("  lo : "); put(lo); new_line;
      end if;
    end if;
  end Sign_Balance;

  procedure Sign_Balance ( x : in out double_double;
                           verbose : in boolean := true ) is

    hi : double_float := hi_part(x);
    lo : double_float := lo_part(x);
    prd : double_float;

  begin
    prd := hi*lo;
    if prd < 0.0 then
      Sign_Balance(hi,lo,verbose);
      x := create(hi,lo);
    end if;
  end Sign_Balance;

  function Is_Sign_Balanced ( x : double_double ) return boolean is

    hi : constant double_float := hi_part(x);
    lo : constant double_float := lo_part(x);
    prd : constant double_float := hi*lo;

  begin
    return (prd >= 0.0);
  end Is_Sign_Balanced;

end Bits_of_Doubles;
