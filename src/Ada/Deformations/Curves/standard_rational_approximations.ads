with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;

package Standard_Rational_Approximations is

-- DESCRIPTION :
--   Given a power series, computes a rational approximation for
--   the function represented by the series.

  procedure Denominator_System
              ( numdeg,dendeg : in integer32; 
                cff : in Standard_Complex_Vectors.Vector;
                mat : out Standard_Complex_Matrices.Matrix;
                rhs : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Defines the coefficient matrix and the right hand side vector
  --   for the denominator of a rational approximation.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   cff      coefficients of the series, of range 0..numdeg+dendeg.

  -- ON RETURN :
  --   mat      square coefficient matrix of range 1..dendeg.
  --   rhs      right hand side vector of range 1..dendeg.

  function Numerator_Coefficients
              ( numdeg,dendeg : integer32;
                dencff,sercff : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the coefficients of the numerator in the rational
  --   approximation as a vector of range 0..numdeg.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   dencff   coefficients of the denominator;
  --   sercff   coefficients of the power series.

  procedure Pade
              ( numdeg,dendeg : in integer32;
                cff : in Standard_Complex_Vectors.Vector;
                numcff,dencff : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   for a series with coefficients given in cff.
  --   The simplest textbook definition is applied in the computation.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   cff      coefficients of the power series.

  -- ON RETURN :
  --   numcff   coefficients of the numerator;
  --   dencff   coefficients of the denominator.

  -- REQUIRED : cff'range = 0..numdeg+dendeg.

  function Evaluate
              ( p : Standard_Complex_Vectors.Vector;
                x : Complex_Number ) return Complex_Number;

  -- DESCRIPTION :
  --   Evaluates the polynomial with coefficients in p at x.

  function Evaluate
              ( num,den : Standard_Complex_Vectors.Vector;
                x : Complex_Number ) return Complex_Number;

  -- DESCRIPTION :
  --   Given the coefficients of the numerator in num
  --   and the coefficients of the denominator in den,
  --   return the value of the rational approximation at x.

end Standard_Rational_Approximations;
