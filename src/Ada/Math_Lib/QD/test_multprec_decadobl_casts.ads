package Test_Multprec_DecaDobl_Casts is

-- DESCRIPTION :
--   Tests type casts between multiprecision floating numbers
--   and deca double numbers.

  procedure Multprec_to_Deca_Double;

  -- DESCRIPTION :
  --   Starts with a given multiprecision float, converts to
  --   deca double and then back to multiprecision float.

  procedure Deca_Double_to_Multprec;

  -- DESCRIPTION :
  --   Starts with a given deca double number, converts to
  --   multiprecision float and then back to deca double.

  procedure Complex_to_Deca_Double;
 
  -- DESCRIPTION :
  --   Prompts for a complex number in standard double precision
  --   and converts this number into deca double precision.
    
  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Multprec_DecaDobl_Casts;
