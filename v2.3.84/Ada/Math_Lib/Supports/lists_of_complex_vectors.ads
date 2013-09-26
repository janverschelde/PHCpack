with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Generic_Lists_of_Vectors;

package Lists_of_Complex_Vectors is 
  new Generic_Lists_of_Vectors(Standard_Complex_Ring,
                               Standard_Complex_Vectors,
                               Standard_Complex_VecVecs);

-- DESCRIPTION :
--   Defines lists of links to vectors of standard complex numbers.
