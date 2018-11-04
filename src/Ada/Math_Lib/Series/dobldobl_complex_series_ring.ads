with Abstract_Ring;
with DoblDobl_Complex_Series;

package DoblDobl_Complex_Series_Ring is
  new Abstract_Ring(DoblDobl_Complex_Series.Link_to_Series,
                    DoblDobl_Complex_Series.Create(0),
                    DoblDobl_Complex_Series.Create(1),
                    DoblDobl_Complex_Series.Create,
                    DoblDobl_Complex_Series.Equal,
                    DoblDobl_Complex_Series.Copy,
                    DoblDobl_Complex_Series."+",
                    DoblDobl_Complex_Series."+",
                    DoblDobl_Complex_Series."-",
                    DoblDobl_Complex_Series."-",
                    DoblDobl_Complex_Series."*",
                    DoblDobl_Complex_Series.Add,
                    DoblDobl_Complex_Series.Sub,
                    DoblDobl_Complex_Series.Min,
                    DoblDobl_Complex_Series.Mul,
                    DoblDobl_Complex_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with complex coefficients,
--   in double double precision.
