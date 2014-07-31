with Generic_NesVecs;
with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;

package QuadDobl_Complex_NesVecs is 
  new Generic_NesVecs(QuadDobl_Complex_Ring,QuadDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines nested vectors over the ring of quad double complex numbers.
