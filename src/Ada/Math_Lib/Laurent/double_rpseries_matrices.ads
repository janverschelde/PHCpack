with Double_rpSeries_Ring;               use Double_rpSeries_Ring;
with Double_rpSeries_Vectors;
with Generic_Matrices;

package Double_rpSeries_Matrices is
  new Generic_Matrices(Double_rpSeries_Ring,Double_rpSeries_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of real power series.
