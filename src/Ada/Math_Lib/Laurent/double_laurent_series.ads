with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;

package Double_Laurent_Series is

-- DESCRIPTION :
--   Defines operations on Laurent series with complex coefficients,
--   in double precision.
--   A Laurent series is defined by a leading exponent and a coefficient
--   vector, with a typical range from 0 to some degree d, d > 0.
--   In this intermediate proof-of-concept package,
--   the data are not encapsulated by types.

  procedure Write ( e : in integer32;
                    c : in Standard_Complex_Vectors.Vector );
  procedure Write ( file : in file_type; e : in integer32;
                    c : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the power series with leading exponent e
  --   and coefficients in c, to standard output or to file.

  procedure Multiply ( d,xe,ye : in integer32;
                       xc,yc : in Standard_Complex_Vectors.Vector;
                       ze : out integer32;
                       zc : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Multiplies two Laurent series.

  -- REQUIRED :
  --   All coefficient vectors range between 0 and d,
  --   where d is the same constant for all series.

  -- ON ENTRY :
  --   d       only coefficients in the range 0 to d are considered;
  --   xe      leading exponent of the first series x;
  --   ye      leading exponent of the second series y;
  --   xc      coefficient vector of the first series x;
  --   yc      coefficient vector of the second series y.

  -- ON RETURN :
  --   ze      leading exponent of the product of x with y;
  --   zc      coefficient vector of the product of x with y.

  procedure Inverse ( d,xe : in integer32;
                      xc : in Standard_Complex_Vectors.Vector;
                      ye : out integer32;
                      yc : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the inverse of a Laurent series.
  --   Any nonzero Laurent series has a nonzero leading coefficient.

  -- REQUIRED : xc'first = 0 = yc'first and xc'last = yc'last.

  -- ON ENTRY :
  --   d       truncation degree of the series;
  --   xe      leading exponent of a nonzero series x;
  --   xc      coefficient vector of the series x.

  -- ON RETURN :
  --   ye      leading exponent of the inverse of x;
  --   yc      coefficient vector of the inverse of x.

  procedure Divide ( d,xe,ye : in integer32;
                     xc,yc : in Standard_Complex_Vectors.Vector;
                     ze : out integer32;
                     zc : out Standard_Complex_Vectors.Vector;
                     iyc : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Divides two Laurent series.

  -- REQUIRED :
  --   All coefficient vectors range between 0 and d,
  --   where d is the same constant for all series.
  --   The second series y should be nonzero.

  -- ON ENTRY :
  --   d       only coefficients in the range 0 to d are considered;
  --   xe      leading exponent of the first series x;
  --   ye      leading exponent of the second series y;
  --   xc      coefficient vector of the first series x;
  --   yc      coefficient vector of the second series y;

  -- ON RETURN :
  --   ze      leading exponent of the product of x with y;
  --   zc      coefficient vector of the product of x with y;
  --   iyc     coefficient vector of the inverse of y, as work space.

  function Is_Zero ( d : integer32;
                     c : Standard_Complex_Vectors.Vector;
                     tol : double_float := 1.0E-15 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if for all i in 0..d,
  --   the magnitude of c(i) is less than the tolerance tol.

  procedure Normalize ( d : in integer32; e : in out integer32;
                        c : in out Standard_Complex_Vectors.Vector;
                        tol : in double_float := 1.0E-15 );

  -- DESCRIPTION :
  --   Normalizes the representation of the Laurent series
  --   so that the leading coefficient is larger than tol.

  -- ON ENTRY :
  --   d       only coefficients in the range 0 to d are considered;
  --   e       leading exponent of the series;
  --   c       coefficient vector of the series;
  --   tol     tolerance to decide the leading degree
  --           when the first coefficient becomes too small.

  -- ON RETURN :
  --   e       augmented exponent for each shift,
  --           if all coefficients are zero, then e is zero as well;
  --   c       shifted coefficient vector so the leading coefficient
  --           is nonzero, or zero for a zero series.

  function Exponent_Gap ( a,b : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the gap beween a and b,
  --   viewed as the leading exponents of two Laurent power series.

  procedure Add ( d,xe,ye : in integer32;
                  xc,yc : in Standard_Complex_Vectors.Vector;
                  ze : out integer32;
                  zc : out Standard_Complex_Vectors.Vector;
                  tol : in double_float := 1.0E-15 );

  -- DESCRIPTION :
  --   Adds two Laurent series, with a tolerance to determine
  --   the leading degree when the leading coefficient is too small.

  -- REQUIRED :
  --   All coefficient vectors range between 0 and d,
  --   where d is the same constant for all series.

  -- ON ENTRY :
  --   d       only coefficients in the range 0 to d are considered;
  --   xe      leading exponent of the first series x;
  --   ye      leading exponent of the second series y;
  --   xc      coefficient vector of the first series x;
  --   yc      coefficient vector of the second series y;
  --   tol     tolerance to decide the leading degree
  --           when the first coefficient becomes too small.

  -- ON RETURN :
  --   ze      leading exponent of the sum of x and y;
  --   zc      coefficient vector of the sum of x and y.

  procedure Subtract ( d,xe,ye : in integer32;
                       xc,yc : in Standard_Complex_Vectors.Vector;
                       ze : out integer32;
                       zc : out Standard_Complex_Vectors.Vector;
                       tol : in double_float := 1.0E-15 );

  -- DESCRIPTION :
  --   Subtracts two Laurent series, with a tolerance to determine
  --   the leading degree when the leading coefficient is too small.

  -- REQUIRED :
  --   All coefficient vectors range between 0 and d,
  --   where d is the same constant for all series.

  -- ON ENTRY :
  --   d       only coefficients in the range 0 to d are considered;
  --   xe      leading exponent of the first series x;
  --   ye      leading exponent of the second series y;
  --   xc      coefficient vector of the first series x;
  --   yc      coefficient vector of the second series y;
  --   tol     tolerance to decide the leading degree
  --           when the first coefficient becomes too small.

  -- ON RETURN :
  --   ze      leading exponent of the difference x - y;
  --   zc      coefficient vector of the difference x - y.

end Double_Laurent_Series;
