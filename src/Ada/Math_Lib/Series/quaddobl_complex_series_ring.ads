with Abstract_Ring;
with QuadDobl_Complex_Series;

package QuadDobl_Complex_Series_Ring is
  new Abstract_Ring(QuadDobl_Complex_Series.Link_to_Series,
                    QuadDobl_Complex_Series.Create(0),
                    QuadDobl_Complex_Series.Create(1),
                    QuadDobl_Complex_Series.Create,
                    QuadDobl_Complex_Series.Equal,
                    QuadDobl_Complex_Series.Copy,
                    QuadDobl_Complex_Series."+",
                    QuadDobl_Complex_Series."+",
                    QuadDobl_Complex_Series."-",
                    QuadDobl_Complex_Series."-",
                    QuadDobl_Complex_Series."*",
                    QuadDobl_Complex_Series.Add,
                    QuadDobl_Complex_Series.Sub,
                    QuadDobl_Complex_Series.Min,
                    QuadDobl_Complex_Series.Mul,
                    QuadDobl_Complex_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with complex coefficients,
--   in quad double precision.
