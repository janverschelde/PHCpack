with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;

package Standard_Dense_Series is

-- DESCRIPTION :
--   A standard dense series is defined by a coefficient vector x,
--   which stores the coefficients of a truncated power series,
--   such as x(0) + x(1)*t + x(2)*t^2 + .. + x(x'last)*t^x'last.
--   Because the series is dense, the powers of the parameter t
--   are stored implicitly as the indices in the coefficient vector.
--   With the definitions and operations in this package,
--   the Standard_Dense_Series_Ring is defined.

  max_order : constant integer32 := 16;

-- Because the type of the number in a ring must be a definite type,
-- we cannot allow truncated series of arbitrary length, but must fix
-- the order of the series at compile time to some number.

  type Series is record
    order : integer32; 
    -- the last exponent in the series, the error is O(t^(order+1))
    cff : Standard_Complex_Vectors.Vector(0..max_order);
    -- only the coefficients in the range 0..order are used
  end record;

-- CREATORS :

  function Create ( i : integer ) return Series;
  function Create ( f : double_float ) return Series;
  function Create ( c : Complex_Number ) return Series;

  -- DESCRIPTION :
  --   Returns a series with one element, the value of i, f, or c.

  function Create ( i : integer; order : integer32 ) return Series;
  function Create ( f : double_float; order : integer32 ) return Series;
  function Create ( c : Complex_Number; order : integer32 ) return Series;

  -- DESCRIPTION :
  --   Returns a series of the given order where the leading element
  --   at position 0 has the value of i, f, or c.

  function Create ( s : Series; order : integer32 ) return Series;

  -- DESCRIPTION :
  --   Returns a series of the given order with coefficients in s.
  --   If order > s'last, then the series on return has all
  --   coefficients of s, padded with zeros.
  --   If order < s'last, then the series on return is a truncation
  --   of the series in s, with only the coefficients in s(0..order).

-- EQUALITY AND COPY :

  function Equal ( s,t : Series ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all coefficients are the same.
  --   If the orders of s and t are different,
  --   then the extra coefficients must be zero for equality to hold.

  procedure Copy ( s : in Series; t : in out Series );

  -- DESCRIPTION :
  --   Copies the coefficients of s to t.
  --   If the order of t is larger than s,
  --   the extra coefficients of t are set to zero.

  -- REQUIRED : s'last <= t'last.

-- ARITHMETICAL OPERATORS :

  function "+" ( s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns a copy of s.

  function "+" ( s,t : Series ) return Series;

  -- DESCRIPTION :
  --   Adds the series s to t.  The order of the series
  --   on return is the maximum of s.order and t.order.
 
  procedure Add ( s : in out Series; t : in Series );

  -- DESCRIPTION :
  --   Adds the series t to s.  If t.order > s.order,
  --   then the terms of t of index > s.order will be ignored.

  function "-" ( s : Series ) return Series;

  -- DESCRIPTION :
  --   The coefficients of series on return has all signs flipped.

  procedure Min ( s : in out Series );

  -- DESCRIPTION :
  --   This is equivalent to s := -s.

  function "-" ( s,t : Series ) return Series;

  -- DESCRIPTION :
  --   Subtracts the series t from s.  The order of the series
  --   on return is the maximum of s.order and t.order.

  procedure Sub ( s : in out Series; t : in Series );

  -- DESCRIPTION :
  --   Subtracts the series t from s.  If t.order > s.order,
  --   then the terms of t of index > s.order will be ignored.

  function "*" ( s,t : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the multiplication of the series s with t.
  --   The order of the series on return is the maximum
  --   of s.order and t.order.
 
  procedure Mul ( s : in out Series; t : in Series );

  -- DESCRIPTION :
  --   Multiplies the series s with t.  If t.order > s.order,
  --   then terms of t with index > s.order will be ignored.

  function Inverse ( s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the inverse of the series with coefficients in s.

  -- REQUIRED : s(0) /= 0.

  function "/" ( s,t : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the series c = s/t.  The order of the series
  --   on return is the maximum of s.order and t.order.

  procedure Div ( s : in out Series; t : in Series );

  -- DESCRIPTION :
  --   Divides the series s by t.  If t.order > s.order,
  --   then terms in t with index > s.order will be ignored.

  -- REQUIRED : t(0) /= 0.

-- DESTRUCTOR :

  procedure Clear ( s : in out Series );

  -- DESCRIPTION :
  --   All coefficients of s are set to zero.
  --   Also the order of s is set to zero.

end Standard_Dense_Series;
