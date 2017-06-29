with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Standard_Complex_Numbers is

-- DESCRIPTION :
--   Offers a Cartesian view on the complex numbers, where the type of the
--   reals are standard hardware double floats.

  type Complex_Number is private;

-- CREATORS :

  function Create ( n : natural32 ) return Complex_Number;
  function Create ( i : integer32 ) return Complex_Number;
  function Create ( i : integer ) return Complex_Number;
  function Create ( f : double_float ) return Complex_Number;
  function Create ( re,im : double_float ) return Complex_Number;
  function Conjugate ( c : Complex_Number ) return Complex_Number;

-- SELECTORS :

  function REAL_PART ( x : Complex_Number ) return double_float;
  function IMAG_PART ( x : Complex_Number ) return double_float;

  function AbsVal ( x : Complex_Number ) return double_float;
    -- x=a+bi, |x|=|a|+|b|
  function AbsVal ( x : Complex_Number ) return Complex_Number;

-- COMPARISON/COPYING :

  function Equal ( x,y : Complex_Number ) return boolean;
  procedure Copy ( x : in Complex_Number; y : in out Complex_Number );

  function "<" ( x,y : Complex_Number ) return boolean;  -- return |x|<|y|
  function ">" ( x,y : Complex_Number ) return boolean;  -- return |x|>|y|

-- ARITHMETHIC OPERATIONS AS FUNCTIONS :

  function "+" ( x : Complex_Number; y : double_float ) return Complex_Number;
  function "-" ( x : Complex_Number; y : double_float ) return Complex_Number;
  function "*" ( x : Complex_Number; y : double_float ) return Complex_Number;
  function "/" ( x : Complex_Number; y : double_float ) return Complex_Number;

  function "+" ( x : double_float; y : Complex_Number ) return Complex_Number;
  function "-" ( x : double_float; y : Complex_Number ) return Complex_Number;
  function "*" ( x : double_float; y : Complex_Number ) return Complex_Number;
  function "/" ( x : double_float; y : Complex_Number ) return Complex_Number;

  function "+" ( x,y : Complex_Number ) return Complex_Number;
  function "+" ( x : Complex_Number )   return Complex_Number;  -- copies x
  function "-" ( x,y : Complex_Number ) return Complex_Number;
  function "-" ( x : Complex_Number )   return Complex_Number;
  function "*" ( x,y : Complex_Number ) return Complex_Number;
  function "/" ( x,y : Complex_Number ) return Complex_Number;

  function "**" ( x : Complex_Number; m : integer ) return Complex_Number;

-- ARITHMETIC OPERATIONS AS PROCEDURES :

  procedure Add ( x : in out Complex_Number; y : in double_float );
  procedure Sub ( x : in out Complex_Number; y : in double_float );
  procedure Mul ( x : in out Complex_Number; y : in double_float );
  procedure Div ( x : in out Complex_Number; y : in double_float );

  procedure Add ( x : in out Complex_Number; y : in Complex_Number );
  procedure Sub ( x : in out Complex_Number; y : in Complex_Number );
  procedure Min ( x : in out Complex_Number );
  procedure Mul ( x : in out Complex_Number; y : in Complex_Number );
  procedure Div ( x : in out Complex_Number; y : in Complex_Number );

-- CHECK IF NaN :

  function is_valid ( x : Complex_Number ) return boolean;

  -- DESCRIPTION :
  --   A complex number is valid if both real and imaginary parts
  --   are valid numbers.  Return true if valid, false otherwise.

-- DESTRUCTOR :

  procedure Clear ( x : in out Complex_Number );

private

  type Complex_Number is
    record
      RE,IM : double_float;
    end record;

end Standard_Complex_Numbers;
