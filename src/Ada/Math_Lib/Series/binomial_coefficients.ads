with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Binomial_Coefficients is

-- DESCRIPTION :
--   Exports a function to compute binomial coefficients.

  function binomial ( n,k : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the binomial coefficient n choose k.

end Binomial_Coefficients;
