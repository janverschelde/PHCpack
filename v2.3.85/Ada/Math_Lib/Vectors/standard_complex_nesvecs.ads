with Generic_NesVecs;
with Standard_Complex_Ring;
with Standard_Complex_Vectors;

package Standard_Complex_NesVecs is 
  new Generic_NesVecs(Standard_Complex_Ring,Standard_Complex_Vectors);

-- DESCRIPTION :
--   Defines nested vectors over the ring of standard complex numbers.
