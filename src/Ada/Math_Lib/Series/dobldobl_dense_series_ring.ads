with Abstract_Ring;
with DoblDobl_Dense_Series;

package DoblDobl_Dense_Series_Ring is
  new Abstract_Ring(DoblDobl_Dense_Series.Series,
                    DoblDobl_Dense_Series.Create(0.0),
                    DoblDobl_Dense_Series.Create(1.0),
                    DoblDobl_Dense_Series.Create,
                    DoblDobl_Dense_Series.Equal,
                    DoblDobl_Dense_Series.Copy,
                    DoblDobl_Dense_Series."+",
                    DoblDobl_Dense_Series."+",
                    DoblDobl_Dense_Series."-",
                    DoblDobl_Dense_Series."-",
                    DoblDobl_Dense_Series."*",
                    DoblDobl_Dense_Series.Add,
                    DoblDobl_Dense_Series.Sub,
                    DoblDobl_Dense_Series.Min,
                    DoblDobl_Dense_Series.Mul,
                    DoblDobl_Dense_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with double double
--   complex coefficients.
