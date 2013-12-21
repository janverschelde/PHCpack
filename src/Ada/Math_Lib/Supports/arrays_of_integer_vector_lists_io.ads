with Standard_Integer_Ring_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;
with Generic_Arrays_of_Vector_Lists_io;

package Arrays_of_Integer_Vector_Lists_io is
  new Generic_Arrays_of_Vector_Lists_io(Standard_Integer_Ring_io,
                                        Standard_Integer_Vectors,
                                        Standard_Integer_Vectors_io,
                                        Standard_Integer_VecVecs,
                                        Lists_of_Integer_Vectors,
                                        Lists_of_Integer_Vectors_io,
                                        Arrays_of_Integer_Vector_Lists);

-- DESCRIPTION :
--   Defines input/output for arrays of lists of links to integer vectors.
