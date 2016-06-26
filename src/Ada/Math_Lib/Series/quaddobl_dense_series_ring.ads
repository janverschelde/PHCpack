with Abstract_Ring;
with QuadDobl_Dense_Series;

package QuadDobl_Dense_Series_Ring is
  new Abstract_Ring(QuadDobl_Dense_Series.Series,
                    QuadDobl_Dense_Series.Create(0.0),
                    QuadDobl_Dense_Series.Create(1.0),
                    QuadDobl_Dense_Series.Create,
                    QuadDobl_Dense_Series.Equal,
                    QuadDobl_Dense_Series.Copy,
                    QuadDobl_Dense_Series."+",
                    QuadDobl_Dense_Series."+",
                    QuadDobl_Dense_Series."-",
                    QuadDobl_Dense_Series."-",
                    QuadDobl_Dense_Series."*",
                    QuadDobl_Dense_Series.Add,
                    QuadDobl_Dense_Series.Sub,
                    QuadDobl_Dense_Series.Min,
                    QuadDobl_Dense_Series.Mul,
                    QuadDobl_Dense_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with quad double
--   complex coefficients.
