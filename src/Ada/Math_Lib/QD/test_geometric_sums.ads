with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Octo_Double_Numbers;               use Octo_Double_Numbers;
with Hexa_Double_Numbers;               use Hexa_Double_Numbers;

package Test_Geometric_Sums is

-- DESCRIPTION :
--   Provides test on computing geometric sums in several precisions.

  procedure Test_Double_Sum
              ( dim : in integer32; ratio : in double_float );

  -- DESCRIPTION :
  --   Sums the geometric sum defined by the ratio and of size dim
  --   and computes the error of this sum.

  procedure Test_Double_Double_Sum
              ( dim : in integer32; ratio : in double_double );

  -- DESCRIPTION :
  --   Sums the geometric sum defined by the ratio and of size dim,
  --   computes the error of this sum, all using double double arithmetic.

  procedure Test_Quad_Double_Sum
              ( dim : in integer32; ratio : in quad_double );

  -- DESCRIPTION :
  --   Sums the geometric sum defined by the ratio and of size dim,
  --   computes the error of this sum, all using quad double arithmetic.

  procedure Test_Octo_Double_Sum
              ( dim : in integer32; ratio : in octo_double );

  -- DESCRIPTION :
  --   Sums the geometric sum defined by the ratio and of size dim,
  --   computes the error of this sum, all using octo double arithmetic.

  procedure Test_Hexa_Double_Sum
              ( dim : in integer32; ratio : in Hexa_double );

  -- DESCRIPTION :
  --   Sums the geometric sum defined by the ratio and of size dim,
  --   computes the error of this sum, all using hexa double arithmetic.

  procedure Test;

  -- DESCRIPTION :
  --   Runs tests in all precisions.

end Test_Geometric_Sums;
