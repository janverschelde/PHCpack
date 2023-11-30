with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_VecVecs;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_VecVecs;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_VecVecs;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_VecVecs;
with Multprec_Complex_Vectors;

package HexaDobl_Complex_Vectors_cv is

-- DESCRIPTION :
--   Converts between vectors of complex hexa doubles 
--   and standard and multiprecision complex vectors.

  function Standard_to_HexaDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return HexaDobl_Complex_Vectors.Vector;
  function Multprec_to_HexaDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return HexaDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a standard or multiprecision vector of complex numbers
  --   into a vector of hexa double complex numbers.

  function HexaDobl_Complex_to_Standard
             ( v : HexaDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function HexaDobl_Complex_to_Multprec
             ( v : HexaDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a vector of hexa double complex numbers
  --   into a vector of standard or multiprecision complex numbers.

  function HexaDobl_Complex_to_DoblDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function HexaDobl_Complex_to_TripDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return TripDobl_Complex_Vectors.Vector;
  function HexaDobl_Complex_to_QuadDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;
  function HexaDobl_Complex_to_PentDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return PentDobl_Complex_Vectors.Vector;
  function HexaDobl_Complex_to_OctoDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return OctoDobl_Complex_Vectors.Vector;
  function HexaDobl_Complex_to_DecaDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return DecaDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns lower precision versions of the hexa double complex
  --   vector v to double double, triple double, quad double,
  --   penta double, octo double, and deca double precision.

  function to_deca_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return DecaDobl_Complex_VecVecs.VecVec;
  function to_deca_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return DecaDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector v converted to deca double precision.

  -- REQUIRED : v /= null;

  function to_octo_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return OctoDobl_Complex_VecVecs.VecVec;
  function to_octo_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return OctoDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector v converted to octo double precision.

  -- REQUIRED : v /= null;

  function to_penta_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return PentDobl_Complex_VecVecs.VecVec;
  function to_penta_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return PentDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector v converted to penta double precision.

  -- REQUIRED : v /= null;

  function to_quad_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return QuadDobl_Complex_VecVecs.VecVec;
  function to_quad_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return QuadDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector v converted to quad double precision.

  -- REQUIRED : v /= null;

  function to_triple_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return TripDobl_Complex_VecVecs.VecVec;
  function to_triple_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return TripDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector v converted to triple double precision.

  -- REQUIRED : v /= null;

  function to_double_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return DoblDobl_Complex_VecVecs.VecVec;
  function to_double_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return DoblDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector v converted to double double precision.

  -- REQUIRED : v /= null;

  function to_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return Standard_Complex_VecVecs.VecVec;
  function to_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return Standard_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the vector v converted to double precision.

  -- REQUIRED : v /= null;

end HexaDobl_Complex_Vectors_cv;
