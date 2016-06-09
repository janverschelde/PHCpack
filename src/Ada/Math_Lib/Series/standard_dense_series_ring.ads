with Abstract_Ring;
with Standard_Dense_Series;

package Standard_Dense_Series_Ring is
  new Abstract_Ring(Standard_Dense_Series.Series,
                    Standard_Dense_Series.Create(0.0),
                    Standard_Dense_Series.Create(1.0),
                    Standard_Dense_Series.Create,
                    Standard_Dense_Series.Equal,
                    Standard_Dense_Series.Copy,
                    Standard_Dense_Series."+",
                    Standard_Dense_Series."+",
                    Standard_Dense_Series."-",
                    Standard_Dense_Series."-",
                    Standard_Dense_Series."*",
                    Standard_Dense_Series.Add,
                    Standard_Dense_Series.Sub,
                    Standard_Dense_Series.Min,
                    Standard_Dense_Series.Mul,
                    Standard_Dense_Series.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with standard
--   complex coefficients.
