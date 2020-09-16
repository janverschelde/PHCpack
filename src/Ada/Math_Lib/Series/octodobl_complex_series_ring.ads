with Abstract_Ring;
with OctoDobl_Complex_Series;

package OctoDobl_Complex_Series_Ring is
  new Abstract_Ring(OctoDobl_Complex_Series.Link_to_Series,
                    OctoDobl_Complex_Series.Create(0),
                    OctoDobl_Complex_Series.Create(1),
                    OctoDobl_Complex_Series.Create,
                    OctoDobl_Complex_Series.Equal,
                    OctoDobl_Complex_Series.Copy,
                    OctoDobl_Complex_Series."+",
                    OctoDobl_Complex_Series."+",
                    OctoDobl_Complex_Series."-",
                    OctoDobl_Complex_Series."-",
                    OctoDobl_Complex_Series."*",
                    OctoDobl_Complex_Series.Add,
                    OctoDobl_Complex_Series.Sub,
                    OctoDobl_Complex_Series.Min,
                    OctoDobl_Complex_Series.Mul,
                    OctoDobl_Complex_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with complex coefficients,
--   in octo double precision.
