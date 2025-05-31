with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Octo_Double_Numbers;               use Octo_Double_Numbers;
with Hexa_Double_Numbers;               use Hexa_Double_Numbers;

package Geometric_Sums is

-- DESCRIPTION :
--   A geometric sum is the sum of a sequence of powers of a fixed ratio.

  function Double_Sum
             ( dim : integer32; ratio : double_float ) return double_float;

  -- DESCRIPTION :
  --   Sums the geometric sum defined by the ratio and of size dim.

  function Double_Double_Sum
             ( dim : integer32; ratio : double_double ) return double_double;

  -- DESCRIPTION :
  --   Sums the geometric sum defined by the ratio and of size dim,
  --   using double double arithmetic.

  function Quad_Double_Sum
             ( dim : integer32; ratio : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Sums the geometric sum defined by the ratio and of size dim,
  --   using quad double arithmetic.

  function Octo_Double_Sum
             ( dim : integer32; ratio : octo_double ) return octo_double;

  -- DESCRIPTION :
  --   Sums the geometric sum defined by the ratio and of size dim,
  --   using octo double arithmetic.

  function Hexa_Double_Sum
             ( dim : integer32; ratio : hexa_double ) return hexa_double;

  -- DESCRIPTION :
  --   Sums the geometric sum defined by the ratio and of size dim,
  --   using hexa double arithmetic.

end Geometric_Sums;
