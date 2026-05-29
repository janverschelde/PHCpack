with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Double_Real_Power_Series is

-- DESCRIPTION :
--   A double real power series consists of a complex vector of
--   coefficients c and a real vector p of powers of t,
--   symbolically represented as c(0) + c(1)*t^p(1) + .. + c(d)*t^p(d),
--   where d is the index of truncation.
--   The powers are positive and sorted in increasing order.
--   If c(0) /= 0, then the multiplicative inverse is well defined
--   and the arithmetical operations + and * define a field.
--   The precision of the numbers is the double precision.

  type Series ( tdx : integer32 ) is record
    cff : Standard_Complex_Vectors.Vector(0..tdx);
    pwt : Standard_Floating_Vectors.Vector(1..tdx);
  end record;
  type Link_to_Series is access Series;

-- CONSTRUCTORS :

  function Make ( x : integer32; tdx : integer32 ) return Series;
  function Make ( x : integer64; tdx : integer32 ) return Series;
  function Make ( x : double_float; tdx : integer32 ) return Series;
  function Make ( x : complex_number; tdx : integer32 ) return Series;

  -- DESCRIPTION :
  --   Given a number x and the truncation index tdx,
  --   returns a series s where s.cff(0) and all else is zero.

  function Make ( x : integer32; tdx : integer32 ) return Link_to_Series;
  function Make ( x : integer64; tdx : integer32 ) return Link_to_Series;
  function Make ( x : double_float; tdx : integer32 ) return Link_to_Series;
  function Make ( x : complex_number; tdx : integer32 ) return Link_to_Series;

  -- DESCRIPTION :
  --   Given a number x and the truncation index tdx,
  --   returns a pointer to a series s where s.cff(0) and all else is zero.

  function Make ( cff : Standard_Complex_Vectors.Vector;
                  pwt : Standard_Floating_Vectors.Vector ) return Series;
  function Make ( cff : Standard_Complex_Vectors.Vector;
                  pwt : Standard_Floating_Vectors.Vector )
                return Link_to_Series;

  -- DESCRIPTION :
  --   Given coefficients in cff and powers of t in pwt,
  --   returns the encapsulation as a series or pointer to a series.

  -- REQUIRED : cff'last = pwt'last, cff'first = 0 and pwt'first = 1.

-- EQUALITY AND COPY :

  function Equal ( a,b : Series;
                   tol : double_float := 1.0E-12 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all coefficients are the same,
  --   using the tolerance tol for a numerical test.
  --   If the truncation indices of a and b are different,
  --   then the extra coefficients must be zero for equality to hold.

  function Equal ( a,b : Link_to_Series;
                   tol : double_float := 1.0E-12 ) return boolean;

  -- DESCRIPTION :
  --   Wraps the test on equality for the content of a and b.
  --   Returns False if a is empty and b is not, or vice versa.

  procedure Copy ( a : in Series; b : in out Series );

  -- DESCRIPTION :
  --   Copies the coefficients and powers of a to b.
  --   All coefficients and powers of a are copied if a.tdx <= b.tdx,
  --   otherwise, only the first b.tdx coefficients and powers of a
  --   are copied to b.

  procedure Copy ( a : in Link_to_Series; b : in out Link_to_Series );

  -- DESCRIPTION :
  --   Clears b and assign to b a series with the same truncation index 
  --   and the same coefficients and powers as a, thus: a deep copy.

-- ARITHMETICAL OPERATORS :

  function "+" ( s : Series; c : complex_number ) return Series;

  -- DESCRIPTION :
  --   Returns a series which is a copy of s,
  --   but with c added to the constant term of s.

  function "+" ( s : Link_to_Series;
                 c : complex_number ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns the series s with the constant c added to s.

  function "+" ( c : complex_number; s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns a series which is a copy of s,
  --   but with c added to the constant term of s.

  function "+" ( c : complex_number;
                 s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns the series s with the constant c added to s.

  procedure Add ( s : in out Series; c : in complex_number );

  -- DESCRIPTION :
  --   This is equivalent to s := s + c.

  procedure Add ( s : in out Link_to_Series;
                  c : in complex_number );

  -- DESCRIPTION :
  --   Does s := s + c without the creation of an extra object.
  --   For memory efficiency, this Add(s,c) is preferred over s := s+c.

  function "+" ( s : Series ) return Series;
  function "+" ( s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a copy of s.

  function "+" ( a,b : Series ) return Series;

  -- DESCRIPTION :
  --   Adds the series a to b.  The truncation index of the series
  --   on return can be as large as the sum of a.tdx and b.tdx.

  function "+" ( a,b : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Adds the series a to b.  The degree of the series
  --   on return can be as large as the sum of a.tdx and b.tdx.
 
  procedure Add ( a : in out Series; b : in Series );

  -- DESCRIPTION :
  --   Adds the series b to a, keeping a.tdx invariant.
  --
  procedure Add ( a : in out Link_to_Series; b : in Link_to_Series );

  -- DESCRIPTION :
  --   Adds the series b to a, keeping the truncation degree of a.

  function "-" ( s : Series; c : complex_number ) return Series;

  -- DESCRIPTION :
  --   Returns the series s - c, subtracting c from the constant
  --   coefficient of s.

  function "-" ( s : Link_to_Series;
                 c : complex_number ) return Link_to_Series; 

  -- DESCRIPTION :
  --   Returns the series s - c, subtracting c from the constant
  --   coefficient of s.

  function "-" ( c : complex_number; s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the series c - s.

  function "-" ( c : complex_number;
                 s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns the series c - s.

  procedure Sub ( s : in out Series; c : in complex_number );

  -- DESCRIPTION :
  --   This is equivalent to s := s - c.

  procedure Sub ( s : in out Link_to_Series;
                  c : in complex_number );

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

  function "-" ( a,b : Series ) return Series;

  -- DESCRIPTION :
  --   Subtracts the series b from a.  The truncation index of the series
  --   on return is the sum of a.tdx and b.tdx.

  function "-" ( a,b : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns a - b.

  procedure Sub ( a : in out Series; b : in Series );

  -- DESCRIPTION :
  --   Subtracts the series b from a.
  --   In contrast to the "-" function, the truncation index of a
  --   does not change, which may lead to a loss of accuracy,
  --   but which could also prevent expression swell.

  procedure Sub ( a : in out Link_to_Series; b : in Link_to_Series );

  -- DESCRIPTION :
  --   Subtracts b from a.
  --   In contrast to the "-" function, the truncation index of a
  --   does not change, which may lead to a loss of accuracy,
  --   but which could also prevent expression swell.

  function "*" ( s : Series; c : complex_number ) return Series;

  -- DESCRIPTION :
  --   Returns s*c.

  function "*" ( s : Link_to_Series;
                 c : complex_number ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns s*c.  Note that, if c is zero, then the series
  --   with same degree as s and zero coefficients is returned.

  function "*" ( c : complex_number; s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns c*s.

  function "*" ( c : complex_number;
                 s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns c*s.  Note that, if c is zero, then the series
  --   with same degree as s and zero coefficients is returned.

  procedure Mul ( s : in out Series; c : in complex_number );

  -- DESCRIPTION :
  --   Is equivalent to s := s*c.

  procedure Mul ( s : in out Link_to_Series; c : in complex_number );

  -- DESCRIPTION :
  --   Does s := s*c, but for memory management,
  --   this procedure should be applied instead of s := s*c,
  --   as the latter will create another copy of s.

  function "*" ( a,b : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the multiplication of the series a with b.
  --   The truncation degree of the series on return
  --   is (a.tdx+1)*(b.tdx+1) - 1.

  function "*" ( a,b : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns the multiplication of the series a with b.
 
  procedure Mul ( a : in out Series; b : in Series );

  -- DESCRIPTION :
  --   Multiplies the series a with b.
  --   The truncation index of a remains invariant,
  --   which may therefore be less accurate,
  --   but also less prone to expression swell.

  procedure Mul ( a : in out Link_to_Series; b : in Link_to_Series );

  -- DESCRIPTION :
  --   Multiplies a with b, with preservation of a.tdx,
  --   except when b = null, then a is cleared.

  function Inverse ( s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the inverse of the series with coefficients in s.

  -- REQUIRED : s.cff(0) /= 0.

  function Inverse ( s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns null if s is null, otherwise returns the series
  --   with content Inverse(s.all).

  -- REQUIRED : s.cff(0) /= 0.

  function "/" ( s : Series; c : complex_number ) return Series;

  -- DESCRIPTION :
  --   Returns the series s/c, where every coefficient of s
  --   is divided by c.

  -- REQUIRED : c /= 0.

  function "/" ( s : Link_to_Series;
                 c : complex_number ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns null if s is null, or otherwise the series s/c
  --   where every coefficient of s is divided by c.

  -- REQUIRED : c /= 0.

  function "/" ( c : complex_number; s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns c/s, obtained by multiplying all coefficients of
  --   the inverse of s with c.

  -- REQUIRED : s.cff(0) /= 0.

  function "/" ( c : complex_number;
                 s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns null if s is null (the result is undefined anyway),
  --   or otherwise returns c/s.

  -- REQUIRED : s.cff(0) /= 0.

  procedure Div ( s : in out Series; c : in complex_number );

  -- DESCRIPTION :
  --   This is equivalent to s := s/c.

  -- REQUIRED : c /= 0.

  procedure Div ( s : in out Link_to_Series; c : in complex_number );

  -- DESCRIPTION :
  --   This is equivalent to s := s/c, but for memory management,
  --   Div(s,c) should be used over s := s/c on a Link_to_Series.

  -- REQUIRED : c /= 0.

  function "/" ( a,b : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the series a/b.  The truncation index of the series
  --   on return is (a.tdx+1)*(b.tdx+1)-1.

  -- REQUIRED : b.cff(0) /= 0.

  function "/" ( a,b : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns null if s or t are null, otherwise,
  --   returns the series s/t, wrapping the other division.

  -- REQUIRED : b.cff(0) /= 0.

  procedure Div ( a : in out Series; b : in Series );

  -- DESCRIPTION :
  --   Divides the series a by b, keeping a.tdx invariant.

  -- REQUIRED : b.cff(0) /= 0.

  procedure Div ( a : in out Link_to_Series; b : in Link_to_Series );

  -- DESCRIPTION :
  --   Divides the series a by b.

  -- REQUIRED : b.cff(0) /= 0.

-- DESTRUCTORS :

  procedure Clear ( s : in out Series );

  -- DESCRIPTION :
  --   Sets all coefficients and all powers to zero.

  procedure Clear ( s : in out Link_to_Series );

  -- DESCRIPTION :
  --   Deallocates the space occupied by s.

end Double_Real_Power_Series;
