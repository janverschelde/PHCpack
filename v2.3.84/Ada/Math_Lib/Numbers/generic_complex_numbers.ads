with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring;
with Abstract_Ring.Field;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Field is new Ring.Field(<>);

package Generic_Complex_Numbers is

-- DESCRIPTION :
--   Offers a Cartesian view on the complex numbers, where the type of the
--   reals is left for instantiation.

  use Ring; use Field;

  type Complex_Number is private;

-- CREATORS :

  function Create ( n : natural32 ) return Complex_Number;
  function Create ( i : integer32 ) return Complex_Number;
  function Create ( i : integer ) return Complex_Number;
  function Create ( f : number ) return Complex_Number;
  function Create ( re,im : number ) return Complex_Number;
  function Conjugate ( c : Complex_Number ) return Complex_Number;

-- SELECTORS :

  function REAL_PART ( x : Complex_Number ) return number;
  function IMAG_PART ( x : Complex_Number ) return number;

  function AbsVal ( x : Complex_Number ) return number;  -- x=a+bi, |x|=|a|+|b|
  function AbsVal ( x : Complex_Number ) return Complex_Number;

-- COMPARISON/COPYING :

  function Equal ( x,y : Complex_Number ) return boolean;
  procedure Copy ( x : in Complex_Number; y : in out Complex_Number );

  function "<" ( x,y : Complex_Number ) return boolean;  -- return |x|<|y|
  function ">" ( x,y : Complex_Number ) return boolean;  -- return |x|>|y|

-- ARITHMETHIC OPERATIONS AS FUNCTIONS :

  function "+" ( x : Complex_Number; y : number ) return Complex_Number;
  function "-" ( x : Complex_Number; y : number ) return Complex_Number;
  function "*" ( x : Complex_Number; y : number ) return Complex_Number;
  function "/" ( x : Complex_Number; y : number ) return Complex_Number;

  function "+" ( x : number; y : Complex_Number ) return Complex_Number;
  function "-" ( x : number; y : Complex_Number ) return Complex_Number;
  function "*" ( x : number; y : Complex_Number ) return Complex_Number;
  function "/" ( x : number; y : Complex_Number ) return Complex_Number;

  function "+" ( x,y : Complex_Number ) return Complex_Number;
  function "+" ( x : Complex_Number )   return Complex_Number;  -- copies x
  function "-" ( x,y : Complex_Number ) return Complex_Number;
  function "-" ( x : Complex_Number )   return Complex_Number;
  function "*" ( x,y : Complex_Number ) return Complex_Number;
  function "/" ( x,y : Complex_Number ) return Complex_Number;

  function "**" ( x : Complex_Number; m : integer ) return Complex_Number;

-- ARITHMETIC OPERATIONS AS PROCEDURES :

  procedure Add ( x : in out Complex_Number; y : in number );
  procedure Sub ( x : in out Complex_Number; y : in number );
  procedure Mul ( x : in out Complex_Number; y : in number );
  procedure Div ( x : in out Complex_Number; y : in number );

  procedure Add ( x : in out Complex_Number; y : in Complex_Number );
  procedure Sub ( x : in out Complex_Number; y : in Complex_Number );
  procedure Min ( x : in out Complex_Number );
  procedure Mul ( x : in out Complex_Number; y : in Complex_Number );
  procedure Div ( x : in out Complex_Number; y : in Complex_Number );

-- DESTRUCTOR :

  procedure Clear ( x : in out Complex_Number );

private

  type Complex_Number is
    record
      RE,IM : number;
    end record;

end Generic_Complex_Numbers;
