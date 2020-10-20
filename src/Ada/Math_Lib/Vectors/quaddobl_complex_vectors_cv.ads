with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
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
  function QuadDobl_Complex_to_DoblDobl
             ( v : QuadDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function QuadDobl_Complex_to_TripDobl
             ( v : QuadDobl_Complex_Vectors.Vector )
             return TripDobl_Complex_Vectors.Vector;
  function QuadDobl_Complex_to_Multprec
             ( v : QuadDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a vector of quad double complex numbers into a vector
  --   of double, double double, triple double, 
  --   or multiprecision complex numbers.

  function to_triple_double
             ( v : QuadDobl_Complex_VecVecs.VecVec )
             return TripDobl_Complex_VecVecs.VecVec;
  function to_triple_double
             ( v : QuadDobl_Complex_VecVecs.Link_to_VecVec )
             return TripDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector v converted to triple double precision.

  -- REQUIRED : v /= null;

  function to_double_double
             ( v : QuadDobl_Complex_VecVecs.VecVec )
             return DoblDobl_Complex_VecVecs.VecVec;
  function to_double_double
             ( v : QuadDobl_Complex_VecVecs.Link_to_VecVec )
             return DoblDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector v converted to double double precision.

  -- REQUIRED : v /= null;

  function to_double
             ( v : QuadDobl_Complex_VecVecs.VecVec )
             return Standard_Complex_VecVecs.VecVec;
  function to_double
             ( v : QuadDobl_Complex_VecVecs.Link_to_VecVec )
             return Standard_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector v converted to double precision.

  -- REQUIRED : v /= null;

end QuadDobl_Complex_Vectors_cv;
