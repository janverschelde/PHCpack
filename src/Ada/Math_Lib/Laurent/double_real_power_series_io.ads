with Double_rpSeries_Vectors;
with Double_rpSeries_Matrices;

package Double_Real_Power_Series_IO is

-- DESCRIPTION :
--   Provides basic output procedure for series with real powers.

  procedure Write ( x : in Double_rpSeries_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes all components of x.

  procedure Write ( A : in Double_rpSeries_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes all elements of A.

end Double_Real_Power_Series_IO;
