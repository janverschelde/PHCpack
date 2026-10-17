with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_Real_Power_Series;

package Test_Double_rpSeries_Arithmetic is

-- DESCRIPTION :
--   Tests the encapsulated arithmetic on real power series
--   via a square root computation.

  function Random_Series
             ( size : integer32; vrblvl : integer32 := 0 )
             return Double_Real_Power_Series.Link_to_Series;

   -- DESCRIPTION :
   --   Returns a random series of the given size.

  procedure Write ( s : in Double_Real_Power_Series.Link_to_Series );

  -- DESCRIPTION :
  --   Writes coefficients and powers of the series s.

  procedure SQRT_Newton_Step
              ( x,A : in Double_Real_Power_Series.Link_to_Series;
                fx,dx : in out Double_Real_Power_Series.Link_to_Series );

  -- DESCRIPTION :
  --   Returns in dx the update of Newton's method to approximate
  --   the square root of A using the equation x^2 - A = 0, returning
  --   dx = fx/(2*x) at the current value for the approximation,
  --   where fx = x^2 - A.

  procedure Test_SQRT ( size : in integer32 );

  -- DESCRIPTION :
  --   Generates a random real power series of the given size and then
  --   runs Newton's method to approximate a square root of the series.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the size of a real power series,
  --   generates a random series, and then applies Newton's method
  --   to compute a square root of the series.

end Test_Double_rpSeries_Arithmetic;
