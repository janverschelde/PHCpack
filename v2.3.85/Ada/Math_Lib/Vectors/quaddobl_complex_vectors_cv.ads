with Standard_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;

package QuadDobl_Complex_Vectors_cv is

-- DESCRIPTION :
--   Converts between vectors of complex quad doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_QuadDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;
  function Multprec_to_QuadDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision vector of complex numbers
  --   into a vector of quad double complex numbers.

  function QuadDobl_Complex_to_Standard
             ( v : QuadDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function QuadDobl_Complex_to_Multprec
             ( v : QuadDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a vector of quad double complex numbers
  --   into a vector of standard or multiprecision complex numbers.

end QuadDobl_Complex_Vectors_cv;
