package Test_Multprec_DoblDobl_Casts is

-- DESCRIPTION :
--   Tests type casts between multiprecision floating numbers
--   and double double numbers.

  procedure Multprec_to_Double_Double;

  -- DESCRIPTION :
  --   Starts with a given multiprecision float, converts to
  --   double double and then back to multiprecision float.

  procedure Double_Double_to_Multprec;

  -- DESCRIPTION :
  --   Starts with a given double double number, converts to
  --   multiprecision float and then back to double double.

  procedure Complex_to_Double_Double;
 
  -- DESCRIPTION :
  --   Prompts for a complex number in standard double precision
  --   and converts this number into double double precision.
    
  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Multprec_DoblDobl_Casts;
