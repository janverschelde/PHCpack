with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;

package Multprec_Complex_Vector_Tools is

-- DESCRIPTION :
--   Here are some tools for multi-precision vectors.

  function Create ( v : Standard_Complex_Vectors.Vector )
                  return Multprec_Complex_Vectors.Vector;

  function Create ( v : Standard_Complex_VecVecs.VecVec )
                  return Multprec_Complex_VecVecs.VecVec;

  function Create ( v : Standard_Complex_VecVecs.Array_of_VecVecs )
                  return Multprec_Complex_VecVecs.Array_of_VecVecs;

  -- DESCRIPTION :
  --   Converts standard complex vectors into multi-precision vectors.

  procedure Set_Size ( v : in out Multprec_Complex_Vectors.Vector;
                       size : in natural32 );

  procedure Set_Size ( v : in out Multprec_Complex_VecVecs.VecVec;
                       size : in natural32 );

  procedure Set_Size ( v : in out Multprec_Complex_VecVecs.Array_of_VecVecs;
                       size : in natural32 );

  -- DESCRIPTION :
  --   Adjusts the size of the numbers in v to the given size.

  function Round ( v : Multprec_Complex_Vectors.Vector )
                 return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Rounds the elements in the multi-precision vector to the nearest
  --   vector with standard complex numbers.

end Multprec_Complex_Vector_Tools;
