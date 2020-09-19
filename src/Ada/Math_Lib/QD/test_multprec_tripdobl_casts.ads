package Test_Multprec_TripDobl_Casts is

-- DESCRIPTION :
--   Tests type casts between multiprecision floating numbers
--   and triple double numbers.

  procedure Multprec_to_Triple_Double;

  -- DESCRIPTION :
  --   Starts with a given multiprecision float, converts to
  --   triple double and then back to multiprecision float.

  procedure Triple_Double_to_Multprec;

  -- DESCRIPTION :
  --   Starts with a given triple double number, converts to
  --   multiprecision float and then back to triple double.

  procedure Complex_to_Triple_Double;
 
  -- DESCRIPTION :
  --   Prompts for a complex number in standard double precision
  --   and converts this number into triple double precision.
    
  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Multprec_TripDobl_Casts;
