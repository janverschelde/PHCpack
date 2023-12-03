with Abstract_Ring;
with HexaDobl_Complex_Series;

package HexaDobl_Complex_Series_Ring is
  new Abstract_Ring(HexaDobl_Complex_Series.Link_to_Series,
                    HexaDobl_Complex_Series.Create(0),
                    HexaDobl_Complex_Series.Create(1),
                    HexaDobl_Complex_Series.Create,
                    HexaDobl_Complex_Series.Equal,
                    HexaDobl_Complex_Series.Copy,
                    HexaDobl_Complex_Series."+",
                    HexaDobl_Complex_Series."+",
                    HexaDobl_Complex_Series."-",
                    HexaDobl_Complex_Series."-",
                    HexaDobl_Complex_Series."*",
                    HexaDobl_Complex_Series.Add,
                    HexaDobl_Complex_Series.Sub,
                    HexaDobl_Complex_Series.Min,
                    HexaDobl_Complex_Series.Mul,
                    HexaDobl_Complex_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with complex coefficients,
--   in hexa double precision.
