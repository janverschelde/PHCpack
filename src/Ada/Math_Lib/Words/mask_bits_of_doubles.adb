package body Mask_Bits_of_Doubles is

  function last_mask ( nbrbits : natural32 ) return unsigned_integer64 is

    res : unsigned_integer64 := 1;
    pwr : unsigned_integer64 := 2;

  begin
    for i in 1..nbrbits-1 loop
      res := res + pwr;
      pwr := 2*pwr;
    end loop;
    return res;
  end last_mask;

  function last_bits ( nbr : unsigned_integer64;
                       nbrbits : natural32 ) return unsigned_integer64 is

    mask : constant unsigned_integer64 := last_mask(nbrbits);

  begin
    return nbr and mask;
  end last_bits;

  function last_bits ( nbr : integer64;
                       nbrbits : natural32 ) return integer64 is
  begin
    return integer64(last_bits(unsigned_integer64(nbr),nbrbits));
  end last_bits;

  function first_bits ( nbr : unsigned_integer64;
                        nbrbits : natural32 ) return unsigned_integer64 is
  begin
     return (nbr - last_bits(nbr, 52 - nbrbits));
  end first_bits;

  function first_bits ( nbr : integer64;
                        nbrbits : natural32 ) return integer64 is

    unbr : constant unsigned_integer64 := unsigned_integer64(nbr);

  begin
     return integer64(unbr - last_bits(unbr, 52 - nbrbits));
  end first_bits;

end Mask_Bits_of_Doubles;
