with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;

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
      if i mod 4 = 0
       then put(" ");
      end if;
      put(integer32(bits(integer32(i))));
    end loop;
  end write_52bits;

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

end Bits_of_Doubles;
