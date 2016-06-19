with Generic_VecVecs;
with Standard_Dense_Series_Ring;
with Standard_Dense_Series_Vectors;

package Standard_Dense_Series_VecVecs is 
  new Generic_VecVecs(Standard_Dense_Series_Ring,
                      Standard_Dense_Series_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with standard complex numbers as coefficients.
