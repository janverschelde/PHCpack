with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Multprec_Floating_Vectors;
with Multprec_Floating_VecVecs;

package Multprec_Floating_Vector_Tools is

-- DESCRIPTION :
--   Here are some tools for multi-precision vectors.

  function Create ( v : Standard_Floating_Vectors.Vector )
                  return Multprec_Floating_Vectors.Vector;

  function Create ( v : Standard_Floating_VecVecs.VecVec )
                  return Multprec_Floating_VecVecs.VecVec;

  function Create ( v : Standard_Floating_VecVecs.Array_of_VecVecs )
                  return Multprec_Floating_VecVecs.Array_of_VecVecs;

  -- DESCRIPTION :
  --   Converts standard floating-point vectors into multi-precision vectors.

  procedure Set_Size ( v : in out Multprec_Floating_Vectors.Vector;
                       size : in natural32 );

  procedure Set_Size ( v : in out Multprec_Floating_VecVecs.VecVec;
                       size : in natural32 );

  procedure Set_Size ( v : in out Multprec_Floating_VecVecs.Array_of_VecVecs;
                       size : in natural32 );

  -- DESCRIPTION :
  --   Adjusts the size of the numbers in v to the given size.

  function Round ( v : Multprec_Floating_Vectors.Vector )
                 return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Rounds the elements in the multi-precision vector to the nearest
  --   vector with standard floating-point numbers.

end Multprec_Floating_Vector_Tools;
