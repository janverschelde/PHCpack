with Abstract_Ring;
with TripDobl_Complex_Series;

package TripDobl_Complex_Series_Ring is
  new Abstract_Ring(TripDobl_Complex_Series.Link_to_Series,
                    TripDobl_Complex_Series.Create(0),
                    TripDobl_Complex_Series.Create(1),
                    TripDobl_Complex_Series.Create,
                    TripDobl_Complex_Series.Equal,
                    TripDobl_Complex_Series.Copy,
                    TripDobl_Complex_Series."+",
                    TripDobl_Complex_Series."+",
                    TripDobl_Complex_Series."-",
                    TripDobl_Complex_Series."-",
                    TripDobl_Complex_Series."*",
                    TripDobl_Complex_Series.Add,
                    TripDobl_Complex_Series.Sub,
                    TripDobl_Complex_Series.Min,
                    TripDobl_Complex_Series.Mul,
                    TripDobl_Complex_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with complex coefficients,
--   in triple double precision.
