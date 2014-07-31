with QuadDobl_Complex_Ring_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_NesVecs;
with Generic_NesVecs_io;

package QuadDobl_Complex_NesVecs_io is 
  new Generic_NesVecs_io(QuadDobl_Complex_Ring_io,
                         QuadDobl_Complex_Vectors,
                         QuadDobl_Complex_Vectors_io,
                         QuadDobl_Complex_NesVecs);

-- DESCRIPTION :
--   Defines input/output of nested vectors quad double complex numbers.
