with Abstract_Ring;
with DecaDobl_Complex_Series;

package DecaDobl_Complex_Series_Ring is
  new Abstract_Ring(DecaDobl_Complex_Series.Link_to_Series,
                    DecaDobl_Complex_Series.Create(0),
                    DecaDobl_Complex_Series.Create(1),
                    DecaDobl_Complex_Series.Create,
                    DecaDobl_Complex_Series.Equal,
                    DecaDobl_Complex_Series.Copy,
                    DecaDobl_Complex_Series."+",
                    DecaDobl_Complex_Series."+",
                    DecaDobl_Complex_Series."-",
                    DecaDobl_Complex_Series."-",
                    DecaDobl_Complex_Series."*",
                    DecaDobl_Complex_Series.Add,
                    DecaDobl_Complex_Series.Sub,
                    DecaDobl_Complex_Series.Min,
                    DecaDobl_Complex_Series.Mul,
                    DecaDobl_Complex_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with complex coefficients,
--   in deca double precision.
