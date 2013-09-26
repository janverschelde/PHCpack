with Generic_VecMats;
with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;

package DoblDobl_Complex_VecMats is 
  new Generic_VecMats(DoblDobl_Complex_Ring,
                      DoblDobl_Complex_Vectors,
                      DoblDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of complex numbers
--   with real and imaginary double double floating-point numbers.
