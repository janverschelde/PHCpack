with Standard_Dense_Series2_Ring_io;
with Standard_Dense_Series2_Vectors;
with Generic_Vectors_io;

package Standard_Dense_Series2_Vectors_io is 
  new Generic_Vectors_io(Standard_Dense_Series2_Ring_io,
                         Standard_Dense_Series2_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of truncated power series
--   with as coefficients complex numbers in double precision.
