with Generic_NesVecs;
with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;

package DoblDobl_Complex_NesVecs is 
  new Generic_NesVecs(DoblDobl_Complex_Ring,DoblDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines nested vectors over the ring of double double complex numbers.
