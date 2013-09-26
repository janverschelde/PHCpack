with Standard_Floating_Ring;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Generic_Lists_of_Vectors;

package Lists_of_Floating_Vectors is 
  new Generic_Lists_of_Vectors(Standard_Floating_Ring,
                               Standard_Floating_Vectors,
                               Standard_Floating_VecVecs);

-- DESCRIPTION :
--   Defines lists of links to vectors of standard floating-point numbers.
