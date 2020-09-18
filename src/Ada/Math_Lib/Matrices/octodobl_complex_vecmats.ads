with Generic_VecMats;
with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Matrices;

package OctoDobl_Complex_VecMats is 
  new Generic_VecMats(OctoDobl_Complex_Ring,
                      OctoDobl_Complex_Vectors,
                      OctoDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of complex numbers
--   with real and imaginary octo double floating-point numbers.
