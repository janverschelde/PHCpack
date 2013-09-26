with Standard_Complex_Ring_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Standard_Complex_NesVecs;
with Generic_NesVecs_io;

package Standard_Complex_NesVecs_io is 
  new Generic_NesVecs_io(Standard_Complex_Ring_io,
                         Standard_Complex_Vectors,
                         Standard_Complex_Vectors_io,
                         Standard_Complex_NesVecs);

-- DESCRIPTION :
--   Defines input/output of nested vectors standard complex numbers.
