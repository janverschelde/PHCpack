with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;

package Multprec_Mathematical_Functions is

-- DESCRIPTION :
--   This package provides some special mathematical functions to ensure
--   a more portable version of the software.

-- CONSTANT :

  PI : constant :=  3.14159_26535_89793_23846_26433_83279_50288;

-- EXPONENTIAL AND LOGARITHMIC FUNCTIONS :

  function EXP ( x : Floating_Number ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the e^x accurate up to as many decimal places as the size
  --   of the number x.

  function LN ( x : Floating_Number ) return Floating_Number;
 
  -- DESCRIPTION :
  --   Returns the natural logarithm of x accurate up to as many decimal
  --   places as the size of the number x.
  --   If x is negative, then the exception numeric_error is raised. 

  function "**" ( x,y : Floating_Number ) return Floating_Number;
  function "**" ( x : double_float;
                  y : Floating_Number ) return Floating_Number;
  function "**" ( x : Floating_Number;
                  y : double_float ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns x**y.

  function LOG2 ( x : Floating_Number ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the binary logarithm.

  function LOG10 ( x : Floating_Number ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the decimal logarithm.

  function SQRT ( x : Floating_Number ) return Floating_Number;

  -- DSECRIPTION :
  --   Returns the square root of x.

-- TRIGONOMETRIC FUNCTIONS :

  function SIN ( x : Floating_Number ) return Floating_Number;
  function COS ( x : Floating_Number ) return Floating_Number;
  function TAN ( x : Floating_Number ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns sine, cosine and tangens of x.

  function ARCSIN ( x : Floating_Number ) return Floating_Number;
  function ARCCOS ( x : Floating_Number ) return Floating_Number;
  function ARCTAN ( x : Floating_Number ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns arcsin, arccos and argtan of x.

  function Radius ( x,y : Floating_Number ) return Floating_Number;
  function Angle  ( x,y : Floating_Number ) return Floating_Number;

end Multprec_Mathematical_Functions;
