with Abstract_Ring;
with Standard_Complex_Series;

package Standard_Complex_Series_Ring is
  new Abstract_Ring(Standard_Complex_Series.Link_to_Series,
                    Standard_Complex_Series.Create(0),
                    Standard_Complex_Series.Create(1),
                    Standard_Complex_Series.Create,
                    Standard_Complex_Series.Equal,
                    Standard_Complex_Series.Copy,
                    Standard_Complex_Series."+",
                    Standard_Complex_Series."+",
                    Standard_Complex_Series."-",
                    Standard_Complex_Series."-",
                    Standard_Complex_Series."*",
                    Standard_Complex_Series.Add,
                    Standard_Complex_Series.Sub,
                    Standard_Complex_Series.Min,
                    Standard_Complex_Series.Mul,
                    Standard_Complex_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with standard
--   complex coefficients.
