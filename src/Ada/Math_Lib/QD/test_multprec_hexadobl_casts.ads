package Test_Multprec_HexaDobl_Casts is

-- DESCRIPTION :
--   Tests type casts between multiprecision floating numbers
--   and hexa double numbers.

  procedure Multprec_to_Hexa_Double;

  -- DESCRIPTION :
  --   Starts with a given multiprecision float,
  --   converts to hexa double and then back to multiprecision float.

  procedure Hexa_Double_to_Multprec;

  -- DESCRIPTION :
  --   Starts with a given hexa double number,
  --   converts to multiprecision float and then back to hexa double.

  procedure Complex_to_Hexa_Double;
 
  -- DESCRIPTION :
  --   Prompts for a complex number in standard double precision
  --   and converts this number into hexa double precision.
    
  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Multprec_HexaDobl_Casts;
