with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;

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

  procedure Test_Mod_Mask_Bits;

  -- DESCRIPTION :
  --   Runs the above four test procedures.

  procedure Test_Last_Bits ( lastbits : in natural32 );

  -- DESCRIPTION :
  --   Tests the extraction of the last bits,
  --   equal to the value of lastbits.

  procedure Test_First_Bits ( firstbits : in natural32 );

  -- DESCRITPION :
  --   Tests the extraction of the first bits,
  --   equal to the value of first bits.

  procedure Test_Mask_Bits;

  -- DESCRIPTION :
  --   Runs the two above test procedures.

  procedure Test_Bit_Split ( x : in double_float );

  -- DESCRIPTION :
  --   Splits the fraction of x in four equal parts,
  --   using the bit representation of the fraction.

  procedure to_Double_Double ( s : in string; x : out double_double;
                               verbose : in boolean := true );

  -- DESCRIPTION :
  --   Reads the string into a double double number
  --   to test for representation errors.

  procedure Test_Particular_Representations;

  -- DESCRIPTION :
  --   Tests the representations of two numbers,
  --   defined in decimal representation of double floats,
  --   which result in double doubles with a nonzero low double,
  --   because the decimal representation leads to a representation
  --   error when converted to binary.
  --   For the particular numbers in this test, multiplying only the
  --   high doubles leads to an error equal to machine precision,
  --   a representation error caused by omitting the low doubles.

  procedure Test_Product ( x,y : in double_double );

  -- DESCRIPTION :
  --   Tests the product of two double doubles via splitting.

  procedure Test_Particular_Product;

  -- DESCRIPTION :
  --   Tests the product of two particular doubles.

  procedure Test_Split_Product;

  -- DESCRIPTION :
  --   Tests the splitting of a 64-bit floating-point number
  --   in four equal parts, for use in making the product.
  --   Two issues matter:
  --   (1) A decimal representation may not have a finite binary expansion.
  --   (2) A bit equal split does not suffice to give accurate results.

  procedure Test_Sign_Balance ( nbr : in double_double );

  -- DESCRIPTION :
  --   Balances the number nbr so the sign of the high double 
  --   is the same as the sign of the low double.

  procedure Test_Sign_Balance;

  -- DESCRIPTION :
  --   Generates two random double doubles and balances
  --   so the sign of the high double is the same
  --   as the sign of the low double.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests.

end Test_Bits_of_Doubles;
