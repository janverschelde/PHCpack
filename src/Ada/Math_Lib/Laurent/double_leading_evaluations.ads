with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Double_Leading_Evaluations is

-- DESCRIPTION :
--   Evaluates Laurent monomials and polynomials at leading terms
--   of series with real powers.

  function Leading_Power 
             ( deg : Standard_Integer_Vectors.Vector;
               pwr : Standard_Floating_Vectors.Vector;
               difidx : integer32 := 0; vrblvl : integer32 := 0 )
             return double_float;

  -- DESCRIPTION :
  --   Returns the evaluation of the leading powers of a series with 
  --   powers in pwr at a monomial with exponents in deg.
  --   If difidx = 0, then the value of monomial is returned,
  --   otherwise, if difidx in deg'range, then the value of
  --   the derivative with respect to difidx is returned.

  -- REQUIRED : deg'range = pwr'range = 1..nvr,
  --   where nvr equals the number of variables,
  --   as pwr(i) is the power for the variable raised to deg(i).

  procedure Leading_Power 
              ( deg : in Standard_Integer_VecVecs.VecVec;
                pwr : in Standard_Floating_Vectors.Vector;
                val : out double_float; idx : out integer32;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Returns the evaluation of the leading powers of a series with 
  --   powers in pwr at a polynomial with exponents in deg.

  -- REQUIRED :
  --   for all i in deg'range: deg(i)'range = pwr'range = 1..nvr,
  --   where nvr equals the number of variables,
  --   as pwr(j) is the power for the variable raised to deg(i)(j).

  -- ON ENTRY :
  --   deg      exponents of the monomials in a Laurent polynomial;
  --   pwr      leading positive powers of a series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   val      value of the minimum power over all monomials;
  --   idx      index in deg'range where the minimum happened.

  function Leading_Coefficient
             ( deg : Standard_Integer_Vectors.Vector;
               cff : Standard_Complex_Vectors.Vector;
               difidx : integer32 := 0; vrblvl : integer32 := 0 )
             return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the monomial with exponents in deg
  --   at the leading coefficients in cff of a series.
  --   If difidx = 0, then the value of monomial is returned,
  --   otherwise, if difidx in deg'range, then the value of
  --   the derivative with respect to difidx is returned.

  -- REQUIRED : deg'range = cff'range, and
  --   if any of the deg(i)'s are negative, then cff(i) is nonzero,
  --   which is implied by cff being leading coefficients.

end Double_Leading_Evaluations;
