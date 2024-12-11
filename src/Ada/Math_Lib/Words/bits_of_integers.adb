with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package body Bits_of_Integers is

  function Bit_Size ( x : integer64 ) return integer32 is

    pwr : natural32;
    bnd : unsigned_integer64;

  begin
    if x <= 0 then
      return -1;
    else
      pwr := 63;
      bnd := 2**natural(pwr);
      while unsigned_integer64(x) < bnd loop
        bnd := bnd / 2;
        pwr := pwr - 1;
      end loop;
      return integer32(pwr)+1;
    end if;
  end Bit_Size;

  procedure Split ( nbr : in unsigned_integer64;
                    high,low : out unsigned_integer64 ) is

    mask : constant unsigned_integer64 := 4294967295;

  begin
    low := nbr and mask;
    high := nbr - low;
  end Split;

  procedure Split ( nbr : in integer64; high,low : out integer64 ) is

    mask : constant unsigned_integer64 := 4294967295;
    unbr : constant unsigned_integer64 := unsigned_integer64(nbr);

  begin
    low := integer64(unbr and mask);
    high := nbr - low;
  end Split;

  procedure Split_30_Bits ( nbr : in unsigned_integer64;
                            high,low,carry : out unsigned_integer64 ) is

    mask30 : constant unsigned_integer64 := 1073741823;
    mask60 : constant unsigned_integer64 := 1152921504606846975;

  begin
    low := nbr and mask30;
    high := (nbr - low) and mask60;
    carry := nbr - high - low;
  end Split_30_Bits;

  procedure Split_30_Bits ( nbr : in integer64;
                            high,low,carry : out integer64 ) is

    unbr : constant unsigned_integer64 := unsigned_integer64(nbr);
    uhigh,ulow,ucarry : unsigned_integer64;

  begin
    Split_30_Bits(unbr,uhigh,ulow,ucarry);
    low := integer64(ulow);
    high := integer64(uhigh);
    carry := integer64(ucarry);
  end Split_30_Bits;

  procedure Split_30_Bit_Words ( nbr : in unsigned_integer64;
                                 high,low,carry : out integer32 ) is

    uhigh,ulow,ucarry : unsigned_integer64;
    base30 : constant unsigned_integer64 := 2**30;
    base60 : constant unsigned_integer64 := 2**60;

  begin
    Split_30_Bits(nbr,uhigh,ulow,ucarry);
    low := integer32(ulow);
    high := integer32(uhigh/base30);
    carry := integer32(ucarry/base60);
  end Split_30_Bit_Words;

  procedure Split_30_Bit_Words ( nbr : in integer64;
                                 high,low,carry : out integer32 ) is

    unbr : constant unsigned_integer64 := unsigned_integer64(nbr);

  begin
    Split_30_Bit_Words(unbr,high,low,carry);
  end Split_30_Bit_Words;

  procedure Quarter ( nbr : in unsigned_integer64;
                      hihi,lohi,hilo,lolo : out unsigned_integer64 ) is

    high,low : unsigned_integer64;
    lowmask : constant unsigned_integer64 := 65535;
    highmask : constant unsigned_integer64 := 281474976710655;

  begin
    Split(nbr,high,low);
    lolo := low and lowmask;
    hilo := low - lolo;
    lohi := high and highmask;
    hihi := high - lohi;
  end Quarter;

  procedure Quarter ( nbr : in integer64;
                      hihi,lohi,hilo,lolo : out integer64 ) is

    high,low : integer64;
    lowmask : constant unsigned_integer64 := 65535;
    highmask : constant unsigned_integer64 := 281474976710655;
    ulow,uhigh : unsigned_integer64;

  begin
    Split(nbr,high,low);
    ulow := unsigned_integer64(low);
    lolo := integer64(ulow and lowmask);
    hilo := low - lolo;
    uhigh := unsigned_integer64(high);
    lohi := integer64(uhigh and highmask);
    hihi := high - lohi;
  end Quarter;

  procedure Split_Bits ( nbr : in integer64; expbits : integer32;
                         high,low : out integer64 ) is

    expo : natural;
    modulus : integer64;

  begin
    if expbits < 0
     then expo := natural(-expbits);
     else expo := natural(expbits);
    end if;
    modulus := 2**expo;
    high := nbr/modulus;
    low := nbr mod modulus;
  end Split_Bits;

end Bits_of_Integers;
