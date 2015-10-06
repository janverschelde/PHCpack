with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Double_Double_Numbers is

-- DESCRIPTION :
--   This package exports working with double doubles to double the
--   precision of the standard (hardware) doubles.
--   The code is based on QD-2.3.9 (Y. Hida, X.S. Li, and D.H. Bailey).

  type double_double is private;

-- CONSTRUCTORS :

  function Create ( i : integer ) return double_double;

  -- DECRIPTION :
  --   Returns the double double representation of the standard integer i.

  function Create ( n : natural32 ) return double_double;
  function Create ( n : natural64 ) return double_double;

  -- DESCRIPTION :
  --   Returns the double double representation of the 32-bit or 64-bit
  --   machine natural number n.

  function Create ( i : integer32 ) return double_double;
  function Create ( i : integer64 ) return double_double;

  -- DESCRIPTION :
  --   Returns the double double representation of the 32-bit or 64-bit
  --   machine integer number i.

  function Create ( f : double_float ) return double_double;

  -- DESCRIPTION :
  --   Returns the double double representation of the double float f.

  function Create ( hi,lo : double_float ) return double_double;

  -- DESCRIPTION :
  --   Returns the double double with high word equal to hi
  --   and low word equal to lo.

-- SELECTORS :

  function hi_part ( d : double_double ) return double_float;

  -- DESCRIPTION :
  --   Returns the most significant or high word of the double double.

  function lo_part ( d : double_double ) return double_float;

  -- DESCRIPTION :
  --   Returns the less significant or low word of the double double.

-- COMPARISON and COPYING :

  function is_zero ( d : double_double ) return boolean;

  -- DESCRIPTION :
  --   Returns true if d is zero, returns false otherwise.

  function is_one ( d : double_double ) return boolean;

  -- DESCRIPTION :
  --   Returns true if d is one, returns false otherwise.

  function is_positive ( d : double_double ) return boolean;

  -- DESCRIPTION : 
  --   Returns true if d is positive, returns false otherwise.

  function is_negative ( d : double_double ) return boolean;

  -- DESCRIPTION : 
  --   Returns true if d is negative, returns false otherwise.

  function equal ( x,y : double_double ) return boolean;
  function equal ( x : double_double; y : double_float ) return boolean;

  function "<" ( x,y : double_double ) return boolean;
  function "<" ( x : double_double; y : double_float ) return boolean;
  function "<" ( x : double_float; y : double_double ) return boolean;
  function "<=" ( x,y : double_double ) return boolean;
  function "<=" ( x : double_double; y : double_float ) return boolean;
  function "<=" ( x : double_float; y : double_double ) return boolean;

  function ">" ( x,y : double_double ) return boolean;
  function ">" ( x : double_double; y : double_float ) return boolean;
  function ">" ( x : double_float; y : double_double ) return boolean;
  function ">=" ( x,y : double_double ) return boolean;
  function ">=" ( x : double_double; y : double_float ) return boolean;
  function ">=" ( x : double_float; y : double_double ) return boolean;

  procedure copy ( x : in double_double; y : in out double_double );

-- Absolute value and type casts :

  function to_int ( x : double_double ) return integer32;
  function to_double ( x : double_double ) return double_float;
  function "abs" ( x : double_double ) return double_double;
  function AbsVal ( x : double_double ) return double_double; -- same as abs
  function floor ( x : double_double ) return double_double;
  function nint ( x : double_float ) return double_float;   -- nearest integer
  function nint ( x : double_double ) return double_double; -- to x

-- ARITHMETICAL OPERATIONS :

  function "+" ( x,y : double_double ) return double_double;
  function "+" ( x : double_double; y : double_float ) return double_double;
  function "+" ( x : double_float; y : double_double ) return double_double;
  function "+" ( x,y : double_float ) return double_double;
  function "+" ( x : double_double ) return double_double; -- returns +x
  procedure Add ( x : in out double_double; y : in double_double ); -- x += y
  procedure Add ( x : in out double_double; y : in double_float );  -- x += y

  function "-" ( x,y : double_double ) return double_double;
  function "-" ( x : double_double; y : double_float ) return double_double;
  function "-" ( x : double_float; y : double_double ) return double_double;
  function "-" ( x,y : double_float ) return double_double;
  function "-" ( x : double_double ) return double_double; -- returns -x
  procedure Min ( x : in out double_double ); -- x := -x
  procedure Sub ( x : in out double_double; y : in double_double ); -- x -= y
  procedure Sub ( x : in out double_double; y : in double_float );  -- x -= y

  function "*" ( x,y : double_double ) return double_double;
  function "*" ( x : double_double; y : double_float ) return double_double;
  function "*" ( x : double_float; y : double_double ) return double_double;
  function "*" ( x,y : double_float ) return double_double;
  procedure Mul ( x : in out double_double; y : in double_double ); -- x *= y
  procedure Mul ( x : in out double_double; y : in double_float );  -- x *= y
  function Mul_pwr2 ( x : double_double; y : double_float ) -- y = 2^k
                    return double_double;
  procedure Mul_pwr2 ( x : in out double_double; y : in double_float );
     -- multiplies x with y, where y is a power of 2

  function "/" ( x,y : double_double ) return double_double;
  function "/" ( x : double_double; y : double_float ) return double_double;
  function "/" ( x : double_float; y : double_double ) return double_double;
  function "/" ( x,y : double_float ) return double_double;
  procedure Div ( x : in out double_double; y : in double_double ); -- x /= y
  procedure Div ( x : in out double_double; y : in double_float );  -- x /= y

  function sqr ( x : double_float ) return double_double;  -- x^2
  function sqr ( x : double_double ) return double_double; -- x^2
  function "**" ( x : double_double; n : integer ) return double_double; -- x^n
  function "**" ( x : double_double; n : integer32 ) return double_double;
  function "**" ( x : double_double; n : integer64 ) return double_double;

  function ldexp ( x : double_double; n : integer ) return double_double;
     -- multiplies the double double in x with 2^n
  function "**" ( x,y : double_double ) return double_double; -- x^y
  function "**" ( x : double_double; y : double_float ) return double_double;
 
  function exp ( x : double_double ) return double_double;   -- returns exp(x)
  function log ( x : double_double ) return double_double;   -- natural log
  function log10 ( x : double_double ) return double_double; -- decimal log

-- DESTRUCTOR :

  procedure clear ( d : in out double_double ); -- sets d to zero

private

  type double_double is record
    hi : double_float;  -- high word, most significant part
    lo : double_float;  -- low word, least significant part
  end record;

end Double_Double_Numbers;
