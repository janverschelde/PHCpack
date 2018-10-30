with Abstract_Ring;
with Standard_Dense_Series2;

package Standard_Dense_Series2_Ring is
  new Abstract_Ring(Standard_Dense_Series2.Link_to_Series,
                    Standard_Dense_Series2.Create(0.0),
                    Standard_Dense_Series2.Create(1.0),
                    Standard_Dense_Series2.Create,
                    Standard_Dense_Series2.Equal,
                    Standard_Dense_Series2.Copy,
                    Standard_Dense_Series2."+",
                    Standard_Dense_Series2."+",
                    Standard_Dense_Series2."-",
                    Standard_Dense_Series2."-",
                    Standard_Dense_Series2."*",
                    Standard_Dense_Series2.Add,
                    Standard_Dense_Series2.Sub,
                    Standard_Dense_Series2.Min,
                    Standard_Dense_Series2.Mul,
                    Standard_Dense_Series2.Clear);

-- DESCRIPTION :
--   Defines the ring of truncated power series with standard
--   complex coefficients.
