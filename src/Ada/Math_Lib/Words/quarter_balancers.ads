with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package Quarter_Balancers is

-- DESCRIPTION :
--   Provides a test to check if two quarters are balanced
--   and a procedure to balance two quarters.

  function Is_Quarter_Balanced
             ( x,y : double_float; vrblvl : integer32 := 0 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the difference in exponents between
  --   two consecutive quarters x and y, is 13, or less.
  --   If vrblvl > 0, then the exponents are written.

  procedure Quarter_Balance
              ( x,y : in out double_float; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Makes x and y balanced by reduction of one bit of x
  --   and addition of one bit to y.
  --   Assumes both x and y are positive.

  procedure Quarter_Balance
              ( x0,x1,x2,x3 : in out double_float;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Checks if the consecutive quarters x0, x1, x2, x3 are balanced
  --   and for any pair that is not balanced, makes the pair balanced.

end Quarter_Balancers;
