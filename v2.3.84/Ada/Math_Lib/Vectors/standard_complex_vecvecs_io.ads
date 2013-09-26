with Standard_Complex_Ring_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Generic_VecVecs_io;

package Standard_Complex_VecVecs_io is 
  new Generic_VecVecs_io(Standard_Complex_Ring_io,
                         Standard_Complex_Vectors,
                         Standard_Complex_Vectors_io,
                         Standard_Complex_VecVecs);

-- DESCRIPTION :
--   Defines input/output of vectors of vectors of standard complex numbers.
