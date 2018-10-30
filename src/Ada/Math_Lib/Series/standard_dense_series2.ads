with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;

package Standard_Dense_Series2 is

-- DESCRIPTION :
--   A standard dense series is defined by a coefficient vector x,
--   which stores the coefficients of a truncated power series,
--   such as x(0) + x(1)*t + x(2)*t^2 + .. + x(x'last)*t^x'last,
--   in standard double hardware precision.
--   Because the series is dense, the powers of the parameter t
--   are stored implicitly as the indices in the coefficient vector.
--   With the definitions and operations in this package,
--   the Standard_Dense_Series_Ring is defined.

-- This is the second design of the series.
-- In the first design, a fixed constant max_deg limited the degrees
-- of all series and max_deg could only be changed at compile time.
-- In this second design we can work with truncated series of fixed
-- degree via Series and series with variable size via Link_to_Series.

  type Series ( deg : integer32 ) is record
    cff : Standard_Complex_Vectors.Vector(0..deg);
  end record;
  type Link_to_Series is access Series;

-- CONSTRUCTORS :

  function Create ( i : integer ) return Series;
  function Create ( f : double_float ) return Series;
  function Create ( c : Complex_Number ) return Series;

  -- DESCRIPTION :
  --   Returns a series with one element, the value of i, f, or c.
  --   The degree of the series on return equals zero.

  function Create ( i : integer ) return Link_to_Series;
  function Create ( f : double_float ) return Link_to_Series;
  function Create ( c : Complex_Number ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a series with one element, the value of i, f, or c.
  --   The degree of the series on return equals zero.

  function Create ( i : integer; deg : integer32 ) return Series;
  function Create ( f : double_float; deg : integer32 ) return Series;
  function Create ( c : Complex_Number; deg : integer32 ) return Series;

  -- DESCRIPTION :
  --   Returns a series of the given degree where the leading element
  --   at position 0 has the value of i, f, or c.

  function Create ( i : integer; deg : integer32 )
                  return Link_to_Series;
  function Create ( f : double_float; deg : integer32 )
                  return Link_to_Series;
  function Create ( c : Complex_Number; deg : integer32 )
                  return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a series of the given degree where the leading element
  --   at position 0 has the value of i, f, or c.

  function Create ( c : Standard_Complex_Vectors.Vector ) return Series;

  -- DESCRIPTION :
  --   Returns a series with coefficients in c,
  --   with degree equal to c'last.

  -- REQUIRED : c'first = 0.

  function Create ( c : Standard_Complex_Vectors.Vector )
                  return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a series with coefficients in c,
  --   with degree equal to c'last.

  -- REQUIRED : c'first = 0.

  function Create ( s : Series; deg : integer32 ) return Series;

  -- DESCRIPTION :
  --   Returns a series of the given degree with coefficients in s.
  --   If deg > s'last, then the series on return has all
  --   coefficients of s, padded with zeros.
  --   If deg < s'last, then the series on return is a truncation
  --   of the series in s, with only the coefficients in s(0..deg).

  function Create ( s : Series; deg : integer32 )
                  return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a series of the given degree with coefficients in s.
  --   If deg > s'last, then the series on return has all
  --   coefficients of s, padded with zeros.
  --   If deg < s'last, then the series on return is a truncation
  --   of the series in s, with only the coefficients in s(0..deg).

  procedure Set_Degree ( s : in out Link_to_Series; deg : in integer32 );

  -- DESCRIPTION :
  --   Sets the degree of the series s to deg,
  --   either by truncating its coefficients if deg < s.deg, or
  --   otherwise by padding the coefficients with deg - s.deg zeroes.

-- EQUALITY AND COPY :

  function Equal ( s,t : Series ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all coefficients are the same.
  --   If the degrees of s and t are different,
  --   then the extra coefficients must be zero for equality to hold.

  function Equal ( s,t : Link_to_Series ) return boolean;

  -- DESCRIPTION :
  --   Wraps the test on equality for the content of s and t.
  --   Returns False if s is empty and t is not, or vice versa.

  procedure Copy ( s : in Series; t : in out Series );

  -- DESCRIPTION :
  --   Copies the coefficients of s to t.
  --   All coefficients of s are copied if s.deg <= t.deg,
  --   otherwise, only the first t.deg coefficients of s
  --   are copied to t.

  procedure Copy ( s : in Link_to_Series; t : in out Link_to_Series );

  -- DESCRIPTION :
  --   Clears t and assign to t a series with the same degree
  --   and the same coefficients as s.
  --   Compared to the Copy on the fixed degree series,
  --   all coefficients of s are copied to t.

-- ORDER :

  function Order ( s : Series; tol : double_float := 0.0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the smallest integer k in the range 0..s.deg
  --   for which AbsVal(s.cff(k)) > tol.
  --   If all coefficients are less than tol, then s.deg+1 is returned.

  function Order ( s : Link_to_Series; tol : double_float := 0.0 )
                 return integer32;

  -- DESCRIPTION :
  --   Returns the smallest integer k in the range 0..s.deg
  --   for which AbsVal(s.cff(k)) > tol.
  --   If all coefficients are less than tol, then s.deg+1 is returned.
  --   Returns -1 if s is null.

-- COMPLEX CONJUGATE :

  function Conjugate ( s : Series ) return Series;

  -- DESCRIPTION :
  --   The complex conjugate of a series s has as coefficients
  --   the complex conjugates of the coefficients of s.

  function Conjugate ( s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns s if s is null, otherwise returns the complex
  --   conjugate of all coefficients in s.

-- ARITHMETICAL OPERATORS :

  function "+" ( s : Series; c : Complex_Number ) return Series;

  -- DESCRIPTION :
  --   Returns a series which is a copy of s,
  --   but with c added to the constant term of s.

  function "+" ( s : Link_to_Series;
                 c : Complex_Number ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns the series s with the constant c added to s.

  function "+" ( c : Complex_Number; s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns a series which is a copy of s,
  --   but with c added to the constant term of s.

  function "+" ( c : Complex_Number;
                 s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns the series s with the constant c added to s.

  procedure Add ( s : in out Series; c : in Complex_Number );

  -- DESCRIPTION :
  --   This is equivalent to s := s + c.

  procedure Add ( s : in out Link_to_Series;
                  c : in Complex_Number );

  -- DESCRIPTION :
  --   Does s := s + c without the creation of an extra object.
  --   For memory efficiency, this Add(s,c) is preferred over s := s+c.

  function "+" ( s : Series ) return Series;
  function "+" ( s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a copy of s.

  function "+" ( s,t : Series ) return Series;

  -- DESCRIPTION :
  --   Adds the series s to t.  The degree of the series
  --   on return is the maximum of s.deg and t.deg.

  function "+" ( s,t : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Adds the series s to t.  The degree of the series
  --   on return is the maximum of s.deg and t.deg.
 
  procedure Add ( s : in out Series; t : in Series );

  -- DESCRIPTION :
  --   Adds the series t to s.  If t.deg > s.deg,
  --   then the coefficients of t corresponding to degrees
  --   higher than s.deg are ignored.
  --
  procedure Add ( s : in out Link_to_Series; t : in Link_to_Series );

  -- DESCRIPTION :
  --   Adds the series t to s.  If t.deg > s.deg,
  --   then the degree of s is extended to the degree of t.

  function "-" ( s : Series; c : Complex_Number ) return Series;

  -- DESCRIPTION :
  --   Returns the series s - c, subtracting c from the constant
  --   coefficient of s.

  function "-" ( s : Link_to_Series;
                 c : Complex_Number ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns the series s - c, subtracting c from the constant
  --   coefficient of s.

  function "-" ( c : Complex_Number; s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the series c - s.

  function "-" ( c : Complex_Number;
                 s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns the series c - s.

  procedure Sub ( s : in out Series; c : in Complex_Number );

  -- DESCRIPTION :
  --   This is equivalent to s := s - c.

  procedure Sub ( s : in out Link_to_Series;
                  c : in Complex_Number );

  -- DESCRIPTION :
  --   Does s := s - c without the creation of an extra object.
  --   For memory efficiency, this Sub(s,c) is preferred over s := s-c.

  function "-" ( s : Series ) return Series;

  -- DESCRIPTION :
  --   The coefficients of series on return has all signs flipped.

  function "-" ( s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns null if s is null or otherwise returns a new series
  --   with the same coefficients as s, but with the signs flipped.

  procedure Min ( s : in out Series );

  -- DESCRIPTION :
  --   This is equivalent to s := -s.

  procedure Min ( s : in out Link_to_Series );

  -- DESCRIPTION :
  --   Flips the signs of the coefficients of s.
  --   For memory efficiency, this is better than s := -s,
  --   if s is of type Link_to_Series.

  function "-" ( s,t : Series ) return Series;

  -- DESCRIPTION :
  --   Subtracts the series t from s.  The degree of the series
  --   on return is the maximum of s.deg and t.deg.

  function "-" ( s,t : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns s - t as a variable degree power series.

  procedure Sub ( s : in out Series; t : in Series );

  -- DESCRIPTION :
  --   Subtracts the series t from s.  If t.deg > s.deg,
  --   then the coefficients of t with degrees higher than t
  --   will be ignored.

  procedure Sub ( s : in out Link_to_Series;
                  t : in Link_to_Series );

  -- DESCRIPTION :
  --   Subtracts t from s.  In case t.deg > s.deg,
  --   the degree of s will be extended to t.deg
  --   and the coefficients of t will be subtracted as well.

  function "*" ( s : Series; c : Complex_Number ) return Series;

  -- DESCRIPTION :
  --   Returns s*c.

  function "*" ( s : Link_to_Series;
                 c : Complex_Number ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns s*c.  Note that, if c is zero, then the series
  --   with same degree as s and zero coefficients is returned.

  function "*" ( c : Complex_Number; s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns c*s.

  function "*" ( c : Complex_Number;
                 s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns c*s.  Note that, if c is zero, then the series
  --   with same degree as s and zero coefficients is returned.

  procedure Mul ( s : in out Series; c : in Complex_Number );

  -- DESCRIPTION :
  --   Is equivalent to s := s*c.

  procedure Mul ( s : in out Link_to_Series; c : in Complex_Number );

  -- DESCRIPTION :
  --   Does s := s*c, but for memory management,
  --   this procedure should be applied instead of s := s*c,
  --   as the latter will create another copy of s.

  function "*" ( s,t : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the multiplication of the series s with t.
  --   The degree of the series on return is the maximum of s.deg and t.deg.

  function "*" ( s,t : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns the multiplication of the series s with t.
  --   The degree of the series on return is the maximum of s.deg and t.deg.
  --   If either s or t is null, then null is returned.
 
  procedure Mul ( s : in out Series; t : in Series );

  -- DESCRIPTION :
  --   Multiplies the series s with t.

  -- REQUIRED : s.deg = t.deg.

  procedure Mul ( s : in out Link_to_Series; t : in Link_to_Series );

  -- DESCRIPTION :
  --   Multiplies t with s and unlike the fixed degree series,
  --   the precision of s will be raised to t.deg if t.deg > s.deg.
  --   This procedure is equivalent to s := s*t, but Mul(s,t) should
  --   be used because of memory management.

  function Inverse ( s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the inverse of the series with coefficients in s.

  -- REQUIRED : s.cff(0) /= 0.

  function Inverse ( s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns null if s is null, otherwise returns the series
  --   with content Inverse(s.all).

  -- REQUIRED : s.cff(0) /= 0.

  function "/" ( s : Series; c : Complex_Number ) return Series;

  -- DESCRIPTION :
  --   Returns the series s/c, where every coefficient of s
  --   is divided by c.

  -- REQUIRED : c /= 0.

  function "/" ( s : Link_to_Series;
                 c : Complex_Number ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns null if s is null, or otherwise the series s/c
  --   where every coefficient of s is divided by c.

  -- REQUIRED : c /= 0.

  function "/" ( c : Complex_Number; s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns c/s, obtained by multiplying all coefficients of
  --   the inverse of s with c.

  -- REQUIRED : s.cff(0) /= 0.

  function "/" ( c : Complex_Number;
                 s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns null if s is null (the result is undefined anyway),
  --   or otherwise returns c/s.

  -- REQUIRED : s.cff(0) /= 0.

  procedure Div ( s : in out Series; c : in Complex_Number );

  -- DESCRIPTION :
  --   This is equivalent to s := s/c.

  -- REQUIRED : c /= 0.

  procedure Div ( s : in out Link_to_Series; c : in Complex_Number );

  -- DESCRIPTION :
  --   This is equivalent to s := s/c, but for memory management,
  --   Div(s,c) should be used over s := s/c on a Link_to_Series.

  -- REQUIRED : c /= 0.

  function "/" ( s,t : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the series c = s/t.  The degree of the series
  --   on return is the maximum of s.deg and t.deg.

  -- REQUIRED : t.cff(0) /= 0.

  function "/" ( s,t : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns null if s or t are null, otherwise,
  --   returns the series c = s/t.  The degree of the series
  --   on return is the maximum of s.deg and t.deg.

  -- REQUIRED : t.cff(0) /= 0.

  procedure Div ( s : in out Series; t : in Series );

  -- DESCRIPTION :
  --   Divides the series s by t.  If t.deg > s.deg,
  --   then the degree of the series of s will be upgraded to t.deg.

  -- REQUIRED : t.cff(0) /= 0.

  procedure Div ( s : in out Link_to_Series; t : in Link_to_Series );

  -- DESCRIPTION :
  --   Divides the series s by t.  If t.deg > s.deg,
  --   then the degree of the series of s will be upgraded to t.deg.

  -- REQUIRED : t.cff(0) /= 0.

  function "**" ( s : Series; p : integer ) return Series;
  function "**" ( s : Series; p : natural32 ) return Series;

  -- DESCRIPTION :
  --   Returns s**p, s to the power p.

  -- REQUIRED : if p < 0, then s.cff(0) /= 0.

  function "**" ( s : Link_to_Series;
                  p : integer ) return Link_to_Series;
  function "**" ( s : Link_to_Series;
                  p : natural32 ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns null if s is null, otherwise,
  --   returns s**p, s to the power p.

  -- REQUIRED : if p < 0, then s.cff(0) /= 0.

-- EVALUATORS :

  function Eval ( s : Series; t : double_float ) return Complex_Number;
  function Eval ( s : Series; t : Complex_Number ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value c(0) + c(1)*t + .. + c(s.deg)*t**s.deg,
  --   where c abbreviates the coefficient vector s.cff.

  function Eval ( s : Link_to_Series;
                  t : double_float ) return Complex_Number;
  function Eval ( s : Link_to_Series;
                  t : Complex_Number ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns zero if s is null, otherwise,
  --   returns the value c(0) + c(1)*t + .. + c(s.deg)*t**s.deg,
  --   where c abbreviates the coefficient vector s.cff.

  function Eval ( s : Series; t : double_float;
                  a,b : integer32 ) return Complex_Number;
  function Eval ( s : Series; t : Complex_Number;
                  a,b : integer32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Evaluates the series using a as the numerator and b as the
  --   numerator of the power for t, so the series starts with
  --   c(0)*t**(a/b) + c(1)*t**((1+a)/b) + ...

  -- REQUIRED : b /= 0 and t /= 0.0 if a < 0.

  function Eval ( s : Link_to_Series; t : double_float;
                  a,b : integer32 ) return Complex_Number;
  function Eval ( s : Link_to_Series; t : Complex_Number;
                  a,b : integer32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns zero if s is null, otherwise,
  --   evaluates the series using a as the numerator and b as the
  --   numerator of the power for t, so the series starts with
  --   c(0)*t**(a/b) + c(1)*t**((1+a)/b) + ...

  -- REQUIRED : b /= 0 and t /= 0.0 if a < 0.

  procedure Filter ( s : in out Series; tol : in double_float );
  procedure Filter ( s : in out Link_to_Series; tol : in double_float );

  -- DESCRIPTION :
  --   All coefficients of s that are less than tol in magnitude 
  --   are set to zero.

-- SHIFT OPERATORS :

  function Shift ( s : Series; c : double_float ) return Series;
  function Shift ( s : Series; c : Complex_Number ) return Series;

  -- DESCRIPTION :
  --   The series on return has the coefficients of the series s,
  --   where the series parameter is replaced by t-c.

  function Shift ( s : Link_to_Series;
                   c : double_float ) return Link_to_Series;
  function Shift ( s : Link_to_Series;
                   c : Complex_Number ) return Link_to_Series;

  -- DESCRIPTION :
  --   if s is null, then null is returned, otherwise,
  --   the series on return has the coefficients of the series s,
  --   where the series parameter is replaced by t-c.

  procedure Shift ( s : in out Series; c : in double_float );
  procedure Shift ( s : in out Series; c : in Complex_Number );

  -- DESCRIPTION :
  --   On return, s = Shift(s,c).

  procedure Shift ( s : in out Link_to_Series; c : in double_float );
  procedure Shift ( s : in out Link_to_Series; c : in Complex_Number );

  -- DESCRIPTION :
  --   On return, s = Shift(s,c).

-- DESTRUCTOR :

  procedure Clear ( s : in out Series );

  -- DESCRIPTION :
  --   All coefficients of s are set to zero.
  --   Also the degree of s is set to zero.

  procedure Clear ( s : in out Link_to_Series );

  -- DESCRIPTION :
  --   The space allocated for s is freed.

end Standard_Dense_Series2;
