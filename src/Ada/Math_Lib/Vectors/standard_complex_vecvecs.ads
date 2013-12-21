with Generic_VecVecs;
with Standard_Complex_Ring;
with Standard_Complex_Vectors;

package Standard_Complex_VecVecs is 
  new Generic_VecVecs(Standard_Complex_Ring,Standard_Complex_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of standard complex numbers.
