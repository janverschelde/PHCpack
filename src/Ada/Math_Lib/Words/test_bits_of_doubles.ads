package Test_Bits_of_Doubles is

-- DESCRIPTION :
--   Test to get the fraction and the exponent from a double,
--   the equivalent constructions of frexp() and ldexp() in C.

  procedure Eat_Last_Bits_of_Pi;

  -- DESCRIPTION :
  --   Removes the last 26 bits of the fraction of the 64-bit
  --   double representation of pi.

  procedure Add_First_Bits_of_Pi;

  -- DESCRIPTION :
  --   Adds the first 26 bits of the fraction of the 64-bit
  --   double representation of pi to the last 26 bits.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests.

end Test_Bits_of_Doubles;
