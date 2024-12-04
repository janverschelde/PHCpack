with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Bits_of_Integers is

-- DESCRIPTION :
--   Provides operations for multiword integer arithmetic.

  procedure Split ( nbr : in integer64; high,low : out integer64 );

  -- DESCRIPTION :
  --   Splits the number nbr into a high and a low word.
  --   The high word stores the first 32 bits of the number.
  --   The low word stores the last 32 bits of the number.

  procedure Quarter ( nbr : in integer64;
                      hihi,lohi,hilo,lolo : out integer64 );

  -- DESCRIPTION :
  --   Quarters the number nbr into four equal words.
  --   The hihi stores the first 16 bits of the number.
  --   The lohi stores the second 16 bits of the number.
  --   The hilo stores the third 16 bits of the number.
  --   The lolo stores the last 16 bits of the number.

end Bits_of_Integers;
