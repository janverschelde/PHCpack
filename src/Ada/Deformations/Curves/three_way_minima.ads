with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;

package Three_Way_Minima is

-- DESCRIPTION :
--   Computes the minimum of three floating-point numbers,
--   in double, double double, or quad double precision.

  function Minimum ( a,b,c : double_float ) return double_float;
  function Minimum ( a,b,c : double_double ) return double_double;
  function Minimum ( a,b,c : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Returns the minimum of a, b, and c.

  procedure Minimum ( a,b,c : in double_float; abcmin : out double_float;
                      anb,bnb,cnb : in out natural32 );
  procedure Minimum ( a,b,c : in double_double; abcmin : out double_double;
                      anb,bnb,cnb : in out natural32 );
  procedure Minimum ( a,b,c : in quad_double; abcmin : out quad_double;
                      anb,bnb,cnb : in out natural32 );

  -- DESCRIPTION :
  --   Returns in abcmin the minimum of a, b, and c.
  --   Updates the frequency for when the minimum is a, b, or c.

  -- ON ENTRY :
  --   a,b,c    numbers in double, double double, or quad double precision;
  --   anb      counts the number of times the first number was the minimum;
  --   bnb      counts the number of times the second number was the minimum;
  --   cnb      counts the number of times the third number was the minimum.

  -- ON RETURN :
  --   abcmin   the minimum of a, b, and c.
  --   anb      updated number of times the first number was the minimum;
  --   bnb      updated number of times the second number was the minimum;
  --   cnb      updated number of times the third number was the minimum.

end Three_Way_Minima;
