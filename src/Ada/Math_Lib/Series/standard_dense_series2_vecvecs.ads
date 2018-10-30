with Generic_VecVecs;
with Standard_Dense_Series2_Ring;
with Standard_Dense_Series2_Vectors;

package Standard_Dense_Series2_VecVecs is 
  new Generic_VecVecs(Standard_Dense_Series2_Ring,
                      Standard_Dense_Series2_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with standard complex numbers as coefficients.
