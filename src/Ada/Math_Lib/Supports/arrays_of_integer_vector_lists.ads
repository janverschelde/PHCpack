with Standard_Integer_Ring;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Lists_of_Integer_Vectors;
with Generic_Arrays_of_Vector_Lists;

package Arrays_of_Integer_Vector_Lists is 
  new Generic_Arrays_of_Vector_Lists(Standard_Integer_Ring,
                                     Standard_Integer_Vectors,
                                     Standard_Integer_VecVecs,
                                     Lists_of_Integer_Vectors);

-- DESCRIPTION :
--   Defines arrays of lists of links to vectors of standard integer numbers.
