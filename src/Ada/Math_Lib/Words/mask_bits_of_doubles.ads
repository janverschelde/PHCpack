with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Mask_Bits_of_Doubles is

-- DESCRIPTION :
--   Provides function to extract the bits of the 52-bit fraction
--   of a 64-bit double precision floating point number.

  type unsigned_integer64 is mod 2**integer64'size;

  function last_mask ( nbrbits : natural32 ) return unsigned_integer64;

  -- DESCRIPTION :
  --   Returns 1 + 2 + .. + 2^n, where n = nbrbits - 1,
  --   to serve as a mask to extract the last nbrbits
  --   of an unsigned 64-bit number with the "and" operator.

  -- REQUIRED : nbrbits > 0.

  function last_bits ( nbr : unsigned_integer64;
                       nbrbits : natural32 ) return unsigned_integer64;

  -- DESCRIPTION :
  --   Returns the last nbr of bits (value in nbrbits) of the number nbr.

  function last_bits ( nbr : integer64;
                       nbrbits : natural32 ) return integer64;

  -- DESCRIPTION :
  --   Returns the last nbr of bits (value in nbrbits) of the number nbr.

  function first_bits ( nbr : unsigned_integer64;
                        nbrbits : natural32 ) return unsigned_integer64;

  -- DESCRIPTION :
  --   Returns the first nbr of bits (value in nbrbits) of the number nbr,
  --   when viewed as a 52-bit integer number, that is:
  --   the number on return has 52 - nbrbits bits as zeros appended.

  function first_bits ( nbr : integer64;
                        nbrbits : natural32 ) return integer64;

  -- DESCRIPTION :
  --   Returns the first nbr of bits (value in nbrbits) of the number nbr,
  --   when viewed as a 52-bit integer number, that is:
  --   the number on return has 52 - nbrbits bits as zeros appended.

end Mask_Bits_of_Doubles;
