with Abstract_Ring;
with PentDobl_Complex_Series;

package PentDobl_Complex_Series_Ring is
  new Abstract_Ring(PentDobl_Complex_Series.Link_to_Series,
                    PentDobl_Complex_Series.Create(0),
                    PentDobl_Complex_Series.Create(1),
                    PentDobl_Complex_Series.Create,
                    PentDobl_Complex_Series.Equal,
                    PentDobl_Complex_Series.Copy,
                    PentDobl_Complex_Series."+",
                    PentDobl_Complex_Series."+",
                    PentDobl_Complex_Series."-",
                    PentDobl_Complex_Series."-",
                    PentDobl_Complex_Series."*",
                    PentDobl_Complex_Series.Add,
                    PentDobl_Complex_Series.Sub,
                    PentDobl_Complex_Series.Min,
                    PentDobl_Complex_Series.Mul,
                    PentDobl_Complex_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with complex coefficients,
--   in penta double precision.
