with Abstract_Ring.Field;
with Double_Real_Power_Series;
with Double_rpSeries_Ring;

package Double_rpSeries_Ring.FField is
  new Double_rpSeries_Ring.Field(Double_Real_Power_Series."<",
                                 Double_Real_Power_Series.">",
                                 Double_Real_Power_Series.AbsVal,
                                 Double_Real_Power_Series."/",
                                 Double_Real_Power_Series.Div);

-- DESCRIPTION :
--   Defines the extension of the ring of series with real powers 
--   and complex coefficients in double precision to a field.
