with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Test_Leading_Evaluations is

-- DESCRIPTION :
--   Tests the evaluation of Laurent monomials at leading terms
--   of series with real powers.

  function Leading_Power 
             ( deg : Standard_Integer_Vectors.Vector;
               pwr : Standard_Floating_Vectors.Vector;
               difidx : integer32 := 0 ) return double_float;

  -- DESCRIPTION :
  --   Returns the evaluation of the leading powers of a series with 
  --   powers in pwr at a monomial with exponents in deg.
  --   If difidx = 0, then the value of monomial is returned,
  --   otherwise, if difidx in deg'range, then the value of
  --   the derivative with respect to difidx is returned.

  -- REQUIRED : deg'range = pwr'range = 1..nvr,
  --   where nvr equals the number of variables,
  --   as pwr(i) is the power for the variable raised to deg(i).

  function Leading_Coefficient
             ( deg : Standard_Integer_Vectors.Vector;
               cff : Standard_Complex_Vectors.Vector;
               difidx : integer32 := 0 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the monomial with exponents in deg
  --   at the leading coefficients in cff of a series.
  --   If difidx = 0, then the value of monomial is returned,
  --   otherwise, if difidx in deg'range, then the value of
  --   the derivative with respect to difidx is returned.

  -- REQUIRED : deg'range = cff'range, and
  --   if any of the deg(i)'s are negative, then cff(i) is nonzero,
  --   which is implied by cff being leading coefficients.

  function Random_Monomial
             ( dim,low,upp : integer32 )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of exponents of dimension dim,
  --   with values randomly generated between low and upp.

  function Random_Leading_Powers
             ( dim : integer32 ) return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..dim with random powers of series
  --   with real positive powers. 

  procedure Test_Monomial_Derivative
              ( deg : in Standard_Integer_Vectors.Vector;
                pwr : in Standard_Floating_Vectors.Vector;
                cff : in Standard_Complex_Vectors.Vector;
                idx : in integer32; err : out double_float );

  -- DESCRIPTION :
  --   Tests the derivative of the monomial with exponents in deg
  --   at the series with leading powers in pwr and coefficients in cff,
  --   with respect to the variable index idx.
  --   On return in err is the magnitude of the error of a random point test.

  procedure Test_Monomial ( dim : in integer32 );
               
  -- DESCRIPTION :
  --   Tests monomial evaluation and differentation in dim many variables
  --    at the leading terms of a series with real positive powers.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and then launches a test.

end Test_Leading_Evaluations;
