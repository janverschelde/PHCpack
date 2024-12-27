with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;

package Test_Sign_Balancers is

-- DESCRIPTION :
--   Tests the procedures to balance multiple double numbers.

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
