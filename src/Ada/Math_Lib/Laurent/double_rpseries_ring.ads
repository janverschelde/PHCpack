with Abstract_Ring;
with Double_Real_Power_Series;

package Double_rpSeries_Ring is
  new Abstract_Ring(Double_Real_Power_Series.Link_to_Series,
                    Double_Real_Power_Series.Make(0.0,tdx=>0),
                    Double_Real_Power_Series.Make(1.0,tdx=>0),
                    Double_Real_Power_Series.Make,
                    Double_Real_Power_Series.Is_Equal,
                    Double_Real_Power_Series.Copy,
                    Double_Real_Power_Series."+",
                    Double_Real_Power_Series."+",
                    Double_Real_Power_Series."-",
                    Double_Real_Power_Series."-",
                    Double_Real_Power_Series."*",
                    Double_Real_Power_Series.Add,
                    Double_Real_Power_Series.Sub,
                    Double_Real_Power_Series.Min,
                    Double_Real_Power_Series.Mul,
                    Double_Real_Power_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of series with real powers and complex coefficients,
--   in double precision.
