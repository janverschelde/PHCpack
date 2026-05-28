with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;

package Double_rpSeries_Operations is

-- DESCRIPTION :
--   A real powered series is defined by a support and a coefficient
--   vector of the same size as the support, plus one for the constant.
--   The support is an increasing sequence of floating-point numbers
--   and the coefficients are complex numbers.
--   The coefficient vectors start at index zero, at the constant term,
--   whereas the positive powers are stored starting at index one,
--   stored in increasing order.

  procedure Sort ( x : in out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Sorts the numbers in x in increasing order.

  procedure Sort ( x : in out Standard_Floating_Vectors.Vector;
                   y : in out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Sorts the numbers in x in increasing order
  --   and swaps the corresponding numbers in y accordingly.

  procedure Shift_Zeros ( cff : in out Standard_Complex_Vectors.Vector;
                          pwt : in out Standard_Floating_Vectors.Vector;
                          tol : in double_float := 1.0E-12 );

  -- DESCRIPTION :
  --   Coefficients which are smaller than the tolerance tol
  --   are shifted towards the end and their corresponding power
  --   is set to zero.

  procedure Normalize ( cff : in out Standard_Complex_Vectors.Vector;
                        pwt : in out Standard_Floating_Vectors.Vector;
                        tol : in double_float := 1.0E-12 );

  -- DESCRIPTION :
  --   Given a sorted sequence of real powers in pwt,
  --   adds the coefficients in cff which correspond to equal powers. 
  --   The tolerance is used to decide on the equality of coefficients.

  procedure Normalize ( cf : in Standard_Complex_VecVecs.VecVec;
                        dg : in Standard_Floating_VecVecs.VecVec;
                        tol : in double_float := 1.0E-12 );

  -- DESCRIPTION :
  --   Applies the normalization to all vectors in (cf, dg).

  function Equal ( acf,bcf : Standard_Complex_Vectors.Vector;
                   apw,bpw : Standard_Floating_Vectors.Vector;
                   tol : double_float := 1.0E-12 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the sum of the componentwise errors of
  --   the coefficients and the powers is less than the tolerance.
  --   Both series are expected to have the same size.

-- BASIC ARITHMETIC :

  procedure Add ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  apw,bpw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in (cff, pwt) the sum of (acf, apw) and (bcf, bpw).
  --   The addition of two random series is defined by the
  --   sum of their constants and the merge sort of the terms.
  --   For the result to be correct, the size of (cff, pwt)
  --   must be the sum of the sizes of (acf, apw) and (bcf, bpw).

  procedure Sub ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  apw,bpw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in (cff, pwt) the result of (acf, apw) - (bcf, bpw).
  --   For the result to be correct, the size of (cff, pwt)
  --   must be the sum of the sizes of (acf, apw) and (bcf, bpw).

  procedure Mul ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  apw,bpw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in (cff, pwt) the product of (acf, apw) and (bcf, bpw).
  --   For the result to be correct, the size of (cff, pwt)
  --   must be the product of the one plus the sizes of (acf, apw)
  --   and (bcf, bpw) minus one (for the constant).

  procedure Inv ( acf : in Standard_Complex_Vectors.Vector;
                  apw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in (cff, pwt) the multiplicative inverse of the
  --   series given in (acf, apw).

  -- REQUIRED : acf(0) /= 0.

  procedure Div ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  apw,bpw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in (cff, pwt) the series (acf, apw) divided by (bcf, bpw),
  --   via multiplication of the inverse of (bcf, bpw).

  -- REQUIRED : bcf(0) /= 0.

-- USEFUL FUNCTIONS :

  function Positive_Minimum_Index
             ( c : Standard_Complex_Vectors.Vector;
               v : Standard_Floating_Vectors.Vector;
               tol : double_float := 1.0E-12 ) return integer32;

  -- DESCRIPTION :
  --   Returns index of the smallest positive number in v,
  --   skipping the entries from which the corresponding c is zero.
  --   The tolerance tol is used to decide if a number is zero.

  function Positive_Minimum
             ( v : Standard_Floating_Vectors.Vector;
               tol : double_float := 1.0E-12 ) return double_float;

  -- DESCRIPTION :
  --   Returns the smallest positive number in v,
  --   using the tolerance tol to decide if a number is zero.

  function Positive_Minimum
             ( c : Standard_Complex_Vectors.Vector;
               v : Standard_Floating_Vectors.Vector;
               tol : double_float := 1.0E-12 ) return double_float;

  -- DESCRIPTION :
  --   Returns the smallest positive number in v,
  --   skipping the entries from which the corresponding c is zero.
  --   The tolerance tol is used to decide if a number is zero.

  function Positive_Minima
             ( c : Standard_Complex_VecVecs.VecVec;
               v : Standard_Floating_VecVecs.VecVec;
               tol : double_float := 1.0E-12 )
             return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the positive minima for each (c(i), v(i)).

  function Coefficient ( c : Standard_Complex_Vectors.Vector;
                         e : Standard_Floating_Vectors.Vector;
                         p : double_float; tol : double_float := 1.0E-12 )
                       return Complex_Number;

  -- DESCRIPTION :
  --   Returns the coefficient c(i) for which e(i) = p,
  --   decided by the test |e(i) - p| < tol.

  function Coefficients ( c : Standard_Complex_VecVecs.VecVec;
                          e : Standard_Floating_VecVecs.VecVec;
                          p : Standard_Floating_Vectors.Vector;
                          tol : double_float := 1.0E-12 )
                       return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the coefficients in each (c(i), e(i))
  --   corresponding to the power p(i).

  function Random_Leading_Powers
             ( dim : integer32 ) return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..dim with random powers of series
  --   with real positive powers. 

  procedure Random_Power_Series
              ( dim : in integer32;
                nbt : in Standard_Integer_Vectors.Vector;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Generates dim power series with random coefficients and real powers.
  --   The next term has power less than twice the power of the previous term.

  -- REQUIRED :
  --   cff'range = nbt'range = pwr'range = 1..dim.

  -- ON ENTRY :
  --   dim      number of power series;
  --   nbt      nbt(i) equals the number of terms in the i-th series;
  --   cff      coefficients of the power series;
  --   pwr      powers of the series.

  function Evaluate_Series
             ( cff : Standard_Complex_VecVecs.VecVec;
               pwr : Standard_Floating_VecVecs.VecVec; tpt : double_float )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Given in cff and pwr are the coefficients and real powers of a series,
  --   and in tpt a value for t.  Returns the value of the series.

end Double_rpSeries_Operations;
