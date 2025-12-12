with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;

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

  procedure Sort ( x : in out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Sorts the numbers in x in increasing order.

  procedure Sort ( x : in out Standard_Floating_Vectors.Vector;
                   y : in out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Sorts the numbers in x in increasing order
  --   and swaps the corresponding numbers in y accordingly.

  procedure Evaluate_Powers
              ( deg : in Standard_Integer_VecVecs.VecVec;
                pwr : in Standard_Floating_Vectors.Vector;
                val : out Standard_Floating_Vectors.Vector;
                idx : out integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Returns the value of the Laurent monomials at real powers
  --   in a sorted vector, also returning the index where the
  --   minimum occurred in the original degrees.

  -- REQUIRED :
  --   for all i in deg'range: deg(i)'range = pwr'range = 1..nvr,
  --   where nvr equals the number of variables,
  --   as pwr(j) is the power for the variable raised to deg(i)(j);
  --   and val'range = deg'range.

  -- ON ENTRY :
  --   deg      exponents of the monomials in a Laurent polynomial;
  --   pwr      leading positive powers of a series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   val      values of the minimum power over all monomials,
  --            sorted in increasing order;
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

  procedure Evaluate_Polynomial
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                xcf : in Standard_Complex_Vectors.Vector;
                xdg : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                idx : out integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Returns the value of the Laurent polynomial p at real powers
  --   in a sorted vector, also returning the index where the
  --   minimum occurred in the original degrees.

  -- REQUIRED :
  --   for all i in pdg'range: pdg(i)'range = xdg'range = 1..nvr,
  --   where nvr equals the number of variables,
  --   as xdg(j) is the power for the variable raised to pdg(i)(j);
  --   and ydg'range = pdg'range.

  -- ON ENTRY :
  --   pcf      coefficients of the Laurent monomials in p;
  --   pct      powers of t with the coefficients of the monomials of p;
  --   pdg      exponents of the Laurent monomials in p;
  --   xcf      leading coefficients of the series;
  --   xdg      leading exponents of the series;
  --   vrblvl   is the verbose level.

   -- ON RETURN :
  --   ycf      coefficients of the evaluated series,
  --            sorted according to the powers in ydg;
  --   ydg      values of the minimum power over all monomials,
  --            sorted in increasing order;
  --   idx      index in pdg'range where the minimum happened.

  procedure Evaluate_System
              ( pcf : in Standard_Complex_VecVecs.VecVec;
                pct : in Standard_Floating_VecVecs.VecVec;
                pdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                xcf : in Standard_Complex_Vectors.Vector;
                xdg : in Standard_Floating_Vectors.Vector;
                ycf : in Standard_Complex_VecVecs.VecVec;
                ydg : in Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Evaluates a system of Laurent polynomials at the leading term
  --   of a series with real powers and complex coefficients.

  -- REQUIRED :
  --   All ranges conform and ycf and ydg are properly allocated.

  -- ON ENTRY :
  --   pcf      pcf(i) has the coefficients of the i-th polynomial;
  --   pct      pct(i) has the power of t in the coefficients;
  --   pdg      pdg(i) has the exponents of the i-th polynomial;
  --   xcf      leading coefficients of the series;
  --   xdg      leading exponents of the series;
  --   vrblvl   is the verbose level.

   -- ON RETURN :
  --   ycf      coefficients of the series, evaluated at i-th polynomial,
  --            sorted according to the powers in ydg;
  --   ydg      values of the minimum power over all monomials,
  --            sorted in increasing order;

end Double_Leading_Evaluations;
