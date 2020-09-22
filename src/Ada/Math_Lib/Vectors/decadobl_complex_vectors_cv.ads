with Standard_Complex_Vectors;
with DecaDobl_Complex_Vectors;
with Multprec_Complex_Vectors;

package DecaDobl_Complex_Vectors_cv is

-- DESCRIPTION :
--   Converts between vectors of complex deca doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_DecaDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return DecaDobl_Complex_Vectors.Vector;
  function Multprec_to_DecaDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return DecaDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision vector of complex numbers
  --   into a vector of deca double complex numbers.

  function DecaDobl_Complex_to_Standard
             ( v : DecaDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function DecaDobl_Complex_to_Multprec
             ( v : DecaDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a vector of deca double complex numbers
  --   into a vector of standard or multiprecision complex numbers.

end DecaDobl_Complex_Vectors_cv;
