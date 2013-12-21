with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Standard_Mathematical_Functions is

-- DESCRIPTION :
--   This package provides some special mathematical functions to ensure
--   a more portable version of the software.

-- CONSTANT :

  PI : constant :=  3.14159_26535_89793_23846_26433_83279_50288;

-- EXPONENTIAL AND LOGARITHMIC FUNCTIONS :

  function EXP ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns e^x, the natural exponential of x.

  function LN ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the natural logarithm of x.

  function "**" ( x,y : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns x**y.

  function LOG2 ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the binary logarithm.

  function LOG10 ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the decimal logarithm.

  function SQRT ( x : double_float ) return double_float;

  -- DSECRIPTION :
  --   Returns the square root of x.

-- TRIGONOMETRIC FUNCTIONS :

  function SIN ( x : double_float ) return double_float;
  function COS ( x : double_float ) return double_float;
  function TAN ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns sine, cosine and tangens of x.

  function ARCSIN ( x : double_float ) return double_float;
  function ARCCOS ( x : double_float ) return double_float;
  function ARCTAN ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns arcsin, arccos and argtan of x.

  function Radius ( x,y : double_float ) return double_float;
  function Angle  ( x,y : double_float ) return double_float;

end Standard_Mathematical_Functions;
