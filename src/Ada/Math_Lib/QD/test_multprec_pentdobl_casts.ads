package Test_Multprec_PentDobl_Casts is

-- DESCRIPTION :
--   Tests type casts between multiprecision floating numbers
--   and penta double numbers.

  procedure Multprec_to_Penta_Double;

  -- DESCRIPTION :
  --   Starts with a given multiprecision float, converts to
  --   penta double and then back to multiprecision float.

  procedure Penta_Double_to_Multprec;

  -- DESCRIPTION :
  --   Starts with a given penta double number, converts to
  --   multiprecision float and then back to penta double.

  procedure Complex_to_Penta_Double;
 
  -- DESCRIPTION :
  --   Prompts for a complex number in standard double precision
  --   and converts this number into penta double precision.
    
  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Multprec_PentDobl_Casts;
