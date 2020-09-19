package Test_Multprec_OctoDobl_Casts is

-- DESCRIPTION :
--   Tests type casts between multiprecision floating numbers
--   and octo double numbers.

  procedure Multprec_to_Octo_Double;

  -- DESCRIPTION :
  --   Starts with a given multiprecision float, converts to
  --   octo double and then back to multiprecision float.

  procedure Octo_Double_to_Multprec;

  -- DESCRIPTION :
  --   Starts with a given octo double number, converts to
  --   multiprecision float and then back to octo double.

  procedure Complex_to_Octo_Double;
 
  -- DESCRIPTION :
  --   Prompts for a complex number in standard double precision
  --   and converts this number into octo double precision.
    
  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Multprec_OctoDobl_Casts;
