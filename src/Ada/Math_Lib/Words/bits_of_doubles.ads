with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;

package Bits_of_Doubles is

-- DESCRIPTION :
--   Provides functions to manipulate the bits of the fraction
--   of a double precision floating point number.

  procedure expand_52bits
              ( bits : out Standard_Natural_Vectors.Vector;
                nbr : in integer64 );

  -- DESCRIPTION :
  --   Given in bits space for 52 bits, in the vector bits of range 0..51,
  --   fills the bits with the binary expansin of the number nbr.

  procedure write_52bits
              ( bits : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes 52 bits in groups of 4.

  function value_52bits
             ( bits : Standard_Natural_Vectors.Vector ) return integer64;

  -- DESCRIPTION :
  --   Given 52 bits in the vector bits, or range 0..51,
  --   returns its value as a 64-bit integer.

  function chop_last_bits
             ( nbr : double_float; lastbits : natural32 )
             return double_float;

  -- DESCRIPTION :
  --   Returns the number nbr with the last bits, equal as lastbits,
  --   chopped off, that is: set to zero.

  -- REQUIRED : lastbits < 52.

  procedure chop_last_bits
             ( nbr : in out double_float; lastbits : in natural32;
               headbits : out Standard_Natural_Vectors.Vector;
               tailbits : out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Sets the last bits of the number nbr to zero.

  -- REQUIRED : lastbits < 52.

  -- ON ENTRY :
  --   nbr     a 64-bit double floating-point number;
  --   lastbits is the number of last bits to set to zero.

  -- ON RETURN :
  --   nbr     the number with the last bits removed;
  --   headbits are the leading bits which remain unchanged;
  --   tailbits are the last bits of the original number nbr.

end Bits_of_Doubles;
