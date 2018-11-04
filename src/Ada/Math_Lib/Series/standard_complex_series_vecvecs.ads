with Generic_VecVecs;
with Standard_Complex_Series_Ring;
with Standard_Complex_Series_Vectors;

package Standard_Complex_Series_VecVecs is 
  new Generic_VecVecs(Standard_Complex_Series_Ring,
                      Standard_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with standard complex numbers as coefficients.
