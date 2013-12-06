with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with Multprec_Complex_Vectors;

package DoblDobl_Complex_Vectors_cv is

-- DESCRIPTION :
--   Converts between vectors of complex double doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_DoblDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Multprec_to_DoblDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision vector of complex numbers
  --   into a vector of double double complex numbers.

  function DoblDobl_Complex_to_Standard
             ( v : DoblDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function DoblDobl_Complex_to_Multprec
             ( v : DoblDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a vector of double double complex numbers
  --   into a vector of standard or multiprecision complex numbers.

end DoblDobl_Complex_Vectors_cv;
