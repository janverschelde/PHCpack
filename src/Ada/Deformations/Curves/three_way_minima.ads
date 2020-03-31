with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;

package Three_Way_Minima is

-- DESCRIPTION :
--   Computes the minimum of three floating-point numbers,
--   in double, double double, or quad double precision.
--   Defines the update of the homotopy continuation parameter t,
--   with respect to the end value and minimum step size bounds.

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

  procedure Bounded_Update
              ( t,step : in out double_float;
                endt,minstep : in double_float );
  procedure Bounded_Update
              ( t,step : in out double_double;
                endt,minstep : in double_float );
  procedure Bounded_Update
              ( t,step : in out quad_double;
                endt,minstep : in double_float );

  -- DESCRIPTION :
  --   Updates the value for t with the step, with respect to the bounds
  --   endt and minstep as follows: t = min(t+step,endt) and if after 
  --   this update |t - endt| <= minstep, then t is set to endt.
  --   The step size is updated as the difference between the value of
  --   t on input and the value of t on output.

  -- ON ENTRY :
  --   t        current value of the homotopy continuation parameter t;
  --   step     current step size;
  --   endt     end value bound for the value of t;
  --   minstep  minimum step size.

  -- ON RETURN :
  --   t        updated value of the homotopy continuation parameter t;
  --   step     updated step size in case t + step exceeds the value of endt,
  --            within the minimum step size bound.

end Three_Way_Minima;
