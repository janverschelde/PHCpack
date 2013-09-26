with Standard_Floating_Ring;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Lists_of_Floating_Vectors;
with Generic_Arrays_of_Vector_Lists;

package Arrays_of_Floating_Vector_Lists is 
  new Generic_Arrays_of_Vector_Lists(Standard_Floating_Ring,
                                     Standard_Floating_Vectors,
                                     Standard_Floating_VecVecs,
                                     Lists_of_Floating_Vectors);

-- DESCRIPTION :
--   Defines arrays of lists of links to vectors of standard floating numbers.
