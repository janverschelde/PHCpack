with Generic_VecMats;
with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;

package QuadDobl_Complex_VecMats is 
  new Generic_VecMats(QuadDobl_Complex_Ring,
                      QuadDobl_Complex_Vectors,
                      QuadDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of complex numbers
--   with real and imaginary quad double floating-point numbers.
