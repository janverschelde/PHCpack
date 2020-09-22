with Standard_Complex_Vectors;
with TripDobl_Complex_Vectors;
with Multprec_Complex_Vectors;

package TripDobl_Complex_Vectors_cv is

-- DESCRIPTION :
--   Converts between vectors of complex triple doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_TripDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return TripDobl_Complex_Vectors.Vector;
  function Multprec_to_TripDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return TripDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision vector of complex numbers
  --   into a vector of triple double complex numbers.

  function TripDobl_Complex_to_Standard
             ( v : TripDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function TripDobl_Complex_to_Multprec
             ( v : TripDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a vector of triple double complex numbers
  --   into a vector of standard or multiprecision complex numbers.

end TripDobl_Complex_Vectors_cv;
