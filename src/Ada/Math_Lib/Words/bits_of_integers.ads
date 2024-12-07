with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Bits_of_Integers is

-- DESCRIPTION :
--   Provides operations for multiword integer arithmetic.

  type unsigned_integer64 is mod 2**integer64'size;

  function Bit_Size ( x : integer64 ) return integer32;

  -- DESCRIPTION :
  --   Returns -1 if x is zero or negative; 
  --   otherwise returns the size of x,
  --   defined as n+1 where 2^n > x >= 2^(n-1).

  procedure Split ( nbr : in unsigned_integer64;
                    high,low : out unsigned_integer64 );
  procedure Split ( nbr : in integer64; high,low : out integer64 );

  -- DESCRIPTION :
  --   Splits the number nbr into a high and a low word.
  --   If not unsigned, nbr must not be negative.
  --   The high word stores the first 32 bits of the number.
  --   The low word stores the last 32 bits of the number.
  --   On return, nbr = high + low.

  procedure Split_30_Bits ( nbr : in unsigned_integer64;
                            high,low,carry : out unsigned_integer64 );
  procedure Split_30_Bits ( nbr : in integer64;
                            high,low,carry : out integer64 );

  -- DESCRIPTION :
  --   Splits the number in a high word, low word, and a carry over.
  --   If not unsigned, nbr must not be negative.
  --   The carry over stores the first four highest bits of the number.
  --   The high word stores the next 30 bits of the number.
  --   The low word stores the last 30 bits of the number.
  --   On return, nbr = carry + high + low.

  procedure Split_30_Bit_Words ( nbr : in unsigned_integer64;
                                 high,low,carry : out integer32 );
  procedure Split_30_Bit_Words ( nbr : in integer64;
                                 high,low,carry : out integer32 );

  -- DESCRIPTION :
  --   Splits the number in three 32-bit integers,
  --   If not unsigned, nbr must not be negative.
  --   where carry stores the first 4 bits,
  --   high the next 30 bits, and low has the last 30 bits.
  --   On return, nbr = carry*2^60 + high*2^30 + low.

  procedure Quarter ( nbr : in unsigned_integer64;
                      hihi,lohi,hilo,lolo : out unsigned_integer64 );
  procedure Quarter ( nbr : in integer64;
                      hihi,lohi,hilo,lolo : out integer64 );

  -- DESCRIPTION :
  --   Quarters the number nbr into four equal words.
  --   If not unsigned, nbr must not be negative.
  --   The hihi stores the first 16 bits of the number.
  --   The lohi stores the second 16 bits of the number.
  --   The hilo stores the third 16 bits of the number.
  --   The lolo stores the last 16 bits of the number.
  --   On return, nbr = hihi + lohi + hilo + lolo.

end Bits_of_Integers;
