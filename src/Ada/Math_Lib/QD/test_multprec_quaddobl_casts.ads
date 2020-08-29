package Test_Multprec_QuadDobl_Casts is

-- DESCRIPTION :
--   Tests type casts between multiprecision floating numbers
--   and double double numbers.

  procedure Multprec_to_Quad_Double;

  -- DESCRIPTION :
  --   Starts with a given multiprecision float, converts to
  --   quad double and then back to multiprecision float.

  procedure Quad_Double_to_Multprec;

  -- DESCRIPTION :
  --   Starts with a given quad double number, converts to
  --   multiprecision float and then back to quad double.

  procedure Complex_to_Quad_Double;

  -- DESCRIPTION :
  --   Prompts for a complex number in standard double precision
  --   and converts the number into quad double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Multprec_QuadDobl_Casts;
