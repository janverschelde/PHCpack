with DoblDobl_Complex_Ring_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_NesVecs;
with Generic_NesVecs_io;

package DoblDobl_Complex_NesVecs_io is 
  new Generic_NesVecs_io(DoblDobl_Complex_Ring_io,
                         DoblDobl_Complex_Vectors,
                         DoblDobl_Complex_Vectors_io,
                         DoblDobl_Complex_NesVecs);

-- DESCRIPTION :
--   Defines input/output of nested vectors double double complex numbers.
