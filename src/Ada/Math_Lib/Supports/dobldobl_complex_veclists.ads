with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with Generic_Lists_of_Vectors;

package DoblDobl_Complex_VecLists is 
  new Generic_Lists_of_Vectors(DoblDobl_Complex_Ring,
                               DoblDobl_Complex_Vectors,
                               DoblDobl_Complex_VecVecs);

-- DESCRIPTION :
--   Defines lists of links to vectors of double double complex numbers.
