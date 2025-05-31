with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;

package Geometric_Inner_Products is

-- DESCRIPTION :
--   A geometric inner product is the sum of the componentwise
--   product of two geometric sequences.

  function Inner_Product ( dim : integer32; rtx,rty : double_float ) 
                         return double_float;
  function Inner_Product ( dim : integer32; rtx,rty : double_double ) 
                         return double_double;
  function Inner_Product ( dim : integer32; rtx,rty : quad_double ) 
                         return quad_double;
  function Inner_Product ( dim : integer32; rtx,rty : octo_double ) 
                         return octo_double;
  function Inner_Product ( dim : integer32; rtx,rty : hexa_double ) 
                         return hexa_double;

  -- DESCRIPTION :
  --   Returns the inner product of two geometric sequences of length dim,
  --   where the ratios of the sequences are rtx and rty.

end Geometric_Inner_Products;
