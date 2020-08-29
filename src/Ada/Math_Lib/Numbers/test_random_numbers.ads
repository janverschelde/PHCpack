with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;

package Test_Random_Numbers is

-- DESCRIPTION :
--   Tests the generation of random numbers.

  function Modulus ( c : Standard_Complex_Numbers.Complex_Number ) 
                   return double_float;

  -- DESCRIPTION :
  --   Returns the modulus of the complex number c.

  procedure Random_Standard_Integer;

  -- DESCRIPTION :
  --   Prompts for lower and upper bounds and generates a random integer
  --   number between the given lower and upper bounds.

  procedure Random_Standard_Complex;

  -- DESCRIPTION :
  --   Prompts for the magnitude and generates a random complex number
  --   with the given magnitude.

  procedure Random_Multprec_Natural;

  -- DESCRIPTION :
  --   Prompts for the size of a multiprecision natural number
  --   and generates a random multiprecision natural number of
  --   the given size.

  procedure Random_Multprec_Integer;

  -- DESCRIPTION :
  --   Prompts for the size of a multiprecision integer number
  --   and generates a random multiprecision integer number of
  --   the given size.

  procedure Random_Multprec_Integer64;

  -- DESCRIPTION :
  --   Prompts for the size of a multiprecision integer number
  --   and generates a random multiprecision integer number of
  --   the given size, with 64-bit words.

  procedure Random_Multprec_Floating;

  -- DESCRIPTION :
  --   Prompts for the size of a multiprecision floating-point number
  --   and generates a random multiprecision floating-point number of
  --   the given size.

  procedure Random_Multprec_Complex;

  -- DESCRIPTION :
  --   Prompts for the size of a multiprecision complex number
  --   and generates a random multiprecision complex number of
  --   the given size.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Random_Numbers;
