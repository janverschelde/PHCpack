with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;

package Quad_Double_Numbers is

-- DESCRIPTION :
--   This package exports working with quad doubles to quadruple the
--   precision of the standard (hardware) doubles.
--   The code is based on QD-2.3.9 (Y. Hida, X.S. Li, and D.H. Bailey).

  type quad_double is private;

-- CONSTRUCTORS :

  function Create ( i : integer ) return quad_double;

  -- DESCRIPTION :
  --   Returns the quad double representation of the standard integer i.

  function Create ( n : natural32 ) return quad_double;
  function Create ( n : natural64 ) return quad_double;

  -- DESCRIPTION :
  --   Returns the quad double representation of the 32-bit or 64-bit
  --   natural number n.

  function Create ( i : integer32 ) return quad_double;
  function Create ( i : integer64 ) return quad_double;

  -- DESCRIPTION :
  --   Returns the quad double representation of the 32-bit or 64-bit
  --   integer number i.

  function Create ( f : double_float ) return quad_double;

  -- DESCRIPTION :
  --   Returns the quad double representation of the double float f.

  function Create ( d : double_double ) return quad_double;

  -- DESCRIPTION :
  --   The high part of the quad double on return equals
  --   the double double d given on entry.

  function Create ( hi,lo : double_double ) return quad_double;

  -- DESCRIPTION :
  --   Returns the quad double with high word equal to hi
  --   and low word equal to lo.

  function Create ( hihi,lohi,hilo,lolo : double_float ) return quad_double;

  -- DESCRIPTION :
  --   Uses the four given numbers in order of significance
  --   to define a quad double.

-- SELECTORS :

  function hi_part ( q : quad_double ) return double_double;

  -- DESCRIPTION :
  --   Returns the most significant or high word of the quad double.

  function hihi_part ( q : quad_double ) return double_float;

  -- DESCRIPTION :
  --   Returns the most significant part of the quad double,
  --   or equivalently: returns hi_part(hi_part(q)).

  function lohi_part ( q : quad_double ) return double_float;

  -- DESCRIPTION :
  --   Returns the second most significant part of the quad double,
  --   or equivalently: returns lo_part(hi_part(q)).

  function lo_part ( q : quad_double ) return double_double;

  -- DESCRIPTION :
  --   Returns the less significant or low word of the quad double.

  function hilo_part ( q : quad_double ) return double_float;

  -- DESCRIPTION :
  --   Returns the high part of the least significant part of q,
  --   or equivalently: returns hi_part(lo_part(q)).

  function lolo_part ( q : quad_double ) return double_float;

  -- DESCRIPTION :
  --   Returns the least significant part of the quad double,
  --   or equivalently: returns lo_part(lo_part(q)).

-- COMPARISON and COPYING :

  function is_zero ( q : quad_double ) return boolean;

  -- DESCRIPTION :
  --   Returns true if q is zero, returns false otherwise.

  function is_one ( q : quad_double ) return boolean;

  -- DESCRIPTION :
  --   Returns true if q is one, returns false otherwise.

  function is_positive ( q : quad_double ) return boolean;

  -- DESCRIPTION : 
  --   Returns true if q is positive, returns false otherwise.

  function is_negative ( q : quad_double ) return boolean;

  -- DESCRIPTION : 
  --   Returns true if q is negative, returns false otherwise.

  function equal ( x,y : quad_double ) return boolean;
  function equal ( x : quad_double; y : double_double ) return boolean;
  function equal ( x : quad_double; y : double_float ) return boolean;

  function "<" ( x,y : quad_double ) return boolean;
  function "<" ( x : quad_double; y : double_double ) return boolean;
  function "<" ( x : double_double; y : quad_double ) return boolean;
  function "<" ( x : quad_double; y : double_float ) return boolean;
  function "<" ( x : double_float; y : quad_double ) return boolean;
  function "<=" ( x,y : quad_double ) return boolean;
  function "<=" ( x : quad_double; y : double_double ) return boolean;
  function "<=" ( x : double_double; y : quad_double ) return boolean;
  function "<=" ( x : quad_double; y : double_float ) return boolean;
  function "<=" ( x : double_float; y : quad_double ) return boolean;

  function ">" ( x,y : quad_double ) return boolean;
  function ">" ( x : quad_double; y : double_double ) return boolean;
  function ">" ( x : double_double; y : quad_double ) return boolean;
  function ">" ( x : quad_double; y : double_float ) return boolean;
  function ">" ( x : double_float; y : quad_double ) return boolean;
  function ">=" ( x,y : quad_double ) return boolean;
  function ">=" ( x : quad_double; y : double_double ) return boolean;
  function ">=" ( x : double_double; y : quad_double ) return boolean;
  function ">=" ( x : quad_double; y : double_float ) return boolean;
  function ">=" ( x : double_float; y : quad_double ) return boolean;

  procedure copy ( x : in quad_double; y : in out quad_double );

-- Absolute value and type casts :

  function to_int ( x : quad_double ) return integer32;
  function to_double ( x : quad_double ) return double_float;
  function to_double_double ( x : quad_double ) return double_double;
  function to_triple_double ( x : quad_double ) return triple_double;
  function "abs" ( x : quad_double ) return quad_double;
  function AbsVal ( x : quad_double ) return quad_double; -- same as abs
  function floor ( x : quad_double ) return quad_double;
  function nint ( x : quad_double ) return quad_double;
    -- nearest integer to x

-- ARITHMETICAL OPERATIONS :

  function "+" ( x,y : quad_double ) return quad_double;
  function "+" ( x : quad_double; y : double_double ) return quad_double;
  function "+" ( x : double_double; y : quad_double ) return quad_double;
  function "+" ( x : quad_double; y : double_float ) return quad_double;
  function "+" ( x : double_float; y : quad_double ) return quad_double;
  function "+" ( x,y : double_double ) return quad_double;
  function "+" ( x : double_double; y : double_float ) return quad_double;
  function "+" ( x : double_float; y : double_double ) return quad_double;
  function "+" ( x,y : double_float ) return quad_double;
  function "+" ( x : quad_double ) return quad_double; -- returns +x
  procedure Add ( x : in out quad_double; y : in quad_double );   -- x += y
  procedure Add ( x : in out quad_double; y : in double_double ); -- x += y
  procedure Add ( x : in out quad_double; y : in double_float );  -- x += y

  function "-" ( x,y : quad_double ) return quad_double;
  function "-" ( x : quad_double; y : double_double ) return quad_double;
  function "-" ( x : double_double; y : quad_double ) return quad_double;
  function "-" ( x : quad_double; y : double_float ) return quad_double;
  function "-" ( x : double_float; y : quad_double ) return quad_double;
  function "-" ( x,y : double_double ) return quad_double;
  function "-" ( x : double_double; y : double_float ) return quad_double;
  function "-" ( x : double_float; y : double_double ) return quad_double;
  function "-" ( x,y : double_float ) return quad_double;
  function "-" ( x : quad_double ) return quad_double; -- returns -x
  procedure Min ( x : in out quad_double ); -- x := -x
  procedure Sub ( x : in out quad_double; y : in quad_double );   -- x -= y
  procedure Sub ( x : in out quad_double; y : in double_double ); -- x -= y
  procedure Sub ( x : in out quad_double; y : in double_float );  -- x -= y

  function "*" ( x,y : quad_double ) return quad_double;
  function "*" ( x : quad_double; y : double_double ) return quad_double;
  function "*" ( x : double_double; y : quad_double ) return quad_double;
  function "*" ( x : quad_double; y : double_float ) return quad_double;
  function "*" ( x : double_float; y : quad_double ) return quad_double;
  function "*" ( x,y : double_double ) return quad_double;
  function "*" ( x : double_double; y : double_float ) return quad_double;
  function "*" ( x : double_float; y : double_double ) return quad_double;
  function "*" ( x,y : double_float ) return quad_double;
  procedure Mul ( x : in out quad_double; y : in quad_double );   -- x *= y
  procedure Mul ( x : in out quad_double; y : in double_double ); -- x *= y
  procedure Mul ( x : in out quad_double; y : in double_float );  -- x *= y
  function Mul_pwr2 ( x : quad_double; y : double_float ) -- y = 2^k
                    return quad_double;
  procedure Mul_pwr2 ( x : in out quad_double; y : in double_float );
     -- multiplies x with y, where y is a power of 2

  function "/" ( x,y : quad_double ) return quad_double;
  function "/" ( x : quad_double; y : double_double ) return quad_double;
  function "/" ( x : double_double; y : quad_double ) return quad_double;
  function "/" ( x : quad_double; y : double_float ) return quad_double;
  function "/" ( x : double_float; y : quad_double ) return quad_double;
  function "/" ( x,y : double_double ) return quad_double;
  function "/" ( x : double_double; y : double_float ) return quad_double;
  function "/" ( x : double_float; y : double_double ) return quad_double;
  function "/" ( x,y : double_float ) return quad_double;
  procedure Div ( x : in out quad_double; y : in quad_double );   -- x /= y
  procedure Div ( x : in out quad_double; y : in double_double ); -- x /= y
  procedure Div ( x : in out quad_double; y : in double_float );  -- x /= y

  function sqr ( x : double_float ) return quad_double;  -- x^2
  function sqr ( x : double_double ) return quad_double; -- x^2
  function sqr ( x : quad_double ) return quad_double;   -- x^2
  function "**" ( x : quad_double; n : integer ) return quad_double; -- x^n
  function "**" ( x : quad_double; n : integer32 ) return quad_double; -- x^n
  function "**" ( x : quad_double; n : integer64 ) return quad_double; -- x^n

  function ldexp ( x : quad_double; n : integer ) return quad_double;
     -- multiplies the quad double in x with 2^n
  function "**" ( x,y : quad_double ) return quad_double; -- x^y
  function "**" ( x : quad_double; y : double_double ) return quad_double;
  function "**" ( x : quad_double; y : double_float ) return quad_double;
 
  function exp ( x : quad_double ) return quad_double;   -- returns exp(x)
  function log ( x : quad_double ) return quad_double;   -- natural log
  function log10 ( x : quad_double ) return quad_double; -- decimal log

-- DESTRUCTOR :

  procedure clear ( q : in out quad_double ); -- sets d to zero

private

-- We can view a quad double as two double doubles with a high and low part,
-- just as the double_double type is defined as a high and a low double.
-- But using two doubles would lead to cumbersome double references,
-- e.g. as q.hi.lo for a quad double q.  Instead of using two double doubles,
-- we use four doubles to represent a quad double.

  type quad_double is record
    hihi : double_float;  -- highest word, most significant part
    lohi : double_float;  -- second highest word, lo_part of hi_part
    hilo : double_float;  -- second lowest word, hi_part of lo_part
    lolo : double_float;  -- lowest word, least significant part
  end record;

end Quad_Double_Numbers;
