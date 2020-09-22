with Standard_Complex_Vectors;
with OctoDobl_Complex_Vectors;
with Multprec_Complex_Vectors;

package OctoDobl_Complex_Vectors_cv is

-- DESCRIPTION :
--   Converts between vectors of complex octo doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_OctoDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return OctoDobl_Complex_Vectors.Vector;
  function Multprec_to_OctoDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return OctoDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision vector of complex numbers
  --   into a vector of octo double complex numbers.

  function OctoDobl_Complex_to_Standard
             ( v : OctoDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function OctoDobl_Complex_to_Multprec
             ( v : OctoDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a vector of octo double complex numbers
  --   into a vector of standard or multiprecision complex numbers.

end OctoDobl_Complex_Vectors_cv;
