with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Triple_Double_Numbers;             use Triple_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Penta_Double_Numbers;              use Penta_Double_Numbers;
with Octo_Double_Numbers;               use Octo_Double_Numbers;
with Deca_Double_Numbers;               use Deca_Double_Numbers;
with Hexa_Double_Numbers;               use Hexa_Double_Numbers;

package Binomial_Coefficients is

-- DESCRIPTION :
--   Exports a function to compute binomial coefficients.

  function binomial ( n,k : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the binomial coefficient n choose k.

  function binomial ( n,k : integer32 ) return double_float;

  -- DESCRIPTION :
  --   Returns the binomial coefficient n choose k,
  --   computed with floating-point arithmetic to avoid
  --   integer arithmetic overflow.

  function binomial ( n,k : integer32 ) return double_double;

  -- DESCRIPTION :
  --   Returns the binomial coefficient n choose k,
  --   computed with double double arithmetic for more accuracy
  --   and to avoid integer arithmetic overflow.

  function binomial ( n,k : integer32 ) return triple_double;

  -- DESCRIPTION :
  --   Returns the binomial coefficient n choose k,
  --   computed with triple double arithmetic for more accuracy
  --   and to avoid integer arithmetic overflow.

  function binomial ( n,k : integer32 ) return quad_double;

  -- DESCRIPTION :
  --   Returns the binomial coefficient n choose k,
  --   computed with quad double arithmetic for even more accuracy
  --   and to avoid integer arithmetic overflow.

  function binomial ( n,k : integer32 ) return penta_double;

  -- DESCRIPTION :
  --   Returns the binomial coefficient n choose k,
  --   computed with penta double arithmetic for even more accuracy
  --   and to avoid integer arithmetic overflow.

  function binomial ( n,k : integer32 ) return octo_double;

  -- DESCRIPTION :
  --   Returns the binomial coefficient n choose k,
  --   computed with octo double arithmetic for even more accuracy
  --   and to avoid integer arithmetic overflow.

  function binomial ( n,k : integer32 ) return deca_double;

  -- DESCRIPTION :
  --   Returns the binomial coefficient n choose k,
  --   computed with deca double arithmetic for even more accuracy
  --   and to avoid integer arithmetic overflow.

  function binomial ( n,k : integer32 ) return hexa_double;

  -- DESCRIPTION :
  --   Returns the binomial coefficient n choose k,
  --   computed with hexa double arithmetic for even more accuracy
  --   and to avoid integer arithmetic overflow.

end Binomial_Coefficients;
