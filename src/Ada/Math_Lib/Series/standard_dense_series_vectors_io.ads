with Standard_Dense_Series_Ring_io;
with Standard_Dense_Series_Vectors;
with Generic_Vectors_io;

package Standard_Dense_Series_Vectors_io is 
  new Generic_Vectors_io(Standard_Dense_Series_Ring_io,
                         Standard_Dense_Series_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of truncated power series
--   with as coefficients complex numbers in double precision.
