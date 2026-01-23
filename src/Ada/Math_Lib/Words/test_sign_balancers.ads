with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;

package Test_Sign_Balancers is

-- DESCRIPTION :
--   Tests the procedures to balance multiple double numbers.

  procedure Test_One_Last_Bit;

  -- DESCRIPTION :
  --   The one last bit of a double float is a number where the fraction
  --   is one bit with an exponent computed so that subtracting this one
  --   last bit from the double results in a fraction that is one bit less.

  procedure Test_Equalize_Signs ( nbr : in double_double );
  procedure Test_Equalize_Signs ( nbr : in quad_double );

  -- DESCRIPTION :
  --   Equalizes the parts of the multiple double nbr
  --   so they all have the same sign.

  procedure Test_DD_Equalize_Signs;

  -- DESCRIPTION :
  --   Generates two random double doubles and equalizes their signs.

  procedure Test_QD_Equalize_Signs;

  -- DESCRIPTION :
  --   Generates two random quad doubles and equalizes their signs.

  procedure Test_Sign_Balance ( nbr : in double_double );
  procedure Test_Sign_Balance ( nbr : in quad_double );
  procedure Test_Sign_Balance ( nbr : in octo_double );
  procedure Test_Sign_Balance ( nbr : in hexa_double );

  -- DESCRIPTION :
  --   Balances the number nbr so all parts have the same sign.

  procedure Test_Sign_DD_Balance;

  -- DESCRIPTION :
  --   Generates two random double doubles and balances
  --   so the sign of the high double is the same
  --   as the sign of the low double.

  procedure Test_Sign_QD_Balance;

  -- DESCRIPTION :
  --   Generates two random quad doubles and balances
  --   so that all parts have the same sign.

  procedure Test_Sign_OD_Balance;

  -- DESCRIPTION :
  --   Generates two random octo doubles and balances
  --   so that all parts have the same sign.

  procedure Test_Sign_HD_Balance;

  -- DESCRIPTION :
  --   Generates two random hexa doubles and balances
  --   so that all parts have the same sign.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests.

end Test_Sign_Balancers;
