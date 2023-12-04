with Generic_VecMats;
with HexaDobl_Complex_Ring;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_Matrices;

package HexaDobl_Complex_VecMats is 
  new Generic_VecMats(HexaDobl_Complex_Ring,
                      HexaDobl_Complex_Vectors,
                      HexaDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of complex numbers
--   with real and imaginary hexa double floating-point numbers.
