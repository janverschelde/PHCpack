with Standard_Complex_Vectors;
with PentDobl_Complex_Vectors;
with Multprec_Complex_Vectors;

package PentDobl_Complex_Vectors_cv is

-- DESCRIPTION :
--   Converts between vectors of complex penta doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_PentDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return PentDobl_Complex_Vectors.Vector;
  function Multprec_to_PentDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return PentDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision vector of complex numbers
  --   into a vector of penta double complex numbers.

  function PentDobl_Complex_to_Standard
             ( v : PentDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function PentDobl_Complex_to_Multprec
             ( v : PentDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a vector of penta double complex numbers
  --   into a vector of standard or multiprecision complex numbers.

end PentDobl_Complex_Vectors_cv;
