package Test_Bits_of_Doubles is

-- DESCRIPTION :
--   Test to get the fraction and the exponent from a double,
--   the equivalent constructions of frexp() and ldexp() in C.

  procedure Eat_Last_Bits_of_Pi;

  -- DESCRIPTION :
  --   Removes the last 26 bits of the fraction of the 64-bit
  --   double representation of pi, with a vector of natural numbers.

  procedure Mod_Last_Bits_of_Pi;

  -- DESCRIPTION :
  --   Removes the last 26 bits of the fraction of the 64-bit
  --   double representation of pi, with modular operations,
  --   without the use of a vector of natural numbers.

  procedure Add_First_Bits_of_Pi;

  -- DESCRIPTION :
  --   Adds the first 26 bits of the fraction of the 64-bit
  --   double representation of pi to the last 26 bits,
  --   with a vector of natural numbers.

  procedure Mod_First_Bits_of_Pi;

  -- DESCRIPTION :
  --   Adds the first 26 bits of the fraction of the 64-bit
  --   double representation of pi to the last 26 bits,
  --   with modular operations, without a vector of natural numbers.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests.

end Test_Bits_of_Doubles;
