package Test_Bits_of_Integers is

-- DESCRIPTION :
--   Tests the multiword arithmetic with 64-bit integers.

  procedure Test_Signed_Split;

  -- DESCRIPTION :
  --   Tests the splitting of a 64-bit integer into a high and low word.

  procedure Test_Unsigned_Split;

  -- DESCRIPTION :
  --   Tests the splitting of an unsigned 64-bit integer
  --   into a high and low word.

  procedure Test_30_Bit_Split;

  -- DESCRIPTION :
  --   Tests the splitting of a 64-bit integer into a low word
  --   with the last 30 bits, a high word with the next 30 bits,
  --   and a carry over, with the 4 highest bits.

  procedure Test_Signed_Quarter;

  -- DESCRIPTION :
  --   Tests the quartering of a 64-bit integer into four words.

  procedure Test_Bit_Split;

  -- DESCRIPTION :
  --   Generates a random 64-bit integer, prompts for the number of bits,
  --   and then splits the number.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests.

end Test_Bits_of_Integers;
