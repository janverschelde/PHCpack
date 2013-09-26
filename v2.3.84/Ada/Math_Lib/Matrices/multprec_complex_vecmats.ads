with Generic_VecMats;
with Multprec_Complex_Ring;
with Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;

package Multprec_Complex_VecMats is 
  new Generic_VecMats(Multprec_Complex_Ring,
                      Multprec_Complex_Vectors,
                      Multprec_Complex_Matrices);

-- DESCRIPTION :
--   Vectors of matrices over the ring of multiprecision complex numbers.
