with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Integer_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions;

package Floating_Integer_Convertors is

-- DESCRIPTION :
--   This package provides routines to convert lists of integer and floating-
--   point vectors into lists of floating-point and integer vectors.
--   The conversion from float to integer is done by merely rounding.

  function Convert ( v : Standard_Integer_Vectors.Vector )
                   return Standard_Floating_Vectors.Vector;
  function Convert ( v : Standard_Floating_Vectors.Vector )
                   return Standard_Integer_Vectors.Vector;

  function Convert ( L : Lists_of_Integer_Vectors.List )
                   return Lists_of_Floating_Vectors.List;
  function Convert ( L : Lists_of_Floating_Vectors.List )
                   return Lists_of_Integer_Vectors.List;

  function Convert ( L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                   return Arrays_of_Floating_Vector_Lists.Array_of_Lists;
  function Convert ( L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                   return Arrays_of_Integer_Vector_Lists.Array_of_Lists;

  function Convert ( m : Integer_Mixed_Subdivisions.Mixed_Cell )
                   return Floating_Mixed_Subdivisions.Mixed_Cell;
  function Convert ( m : Floating_Mixed_Subdivisions.Mixed_Cell )
                   return Integer_Mixed_Subdivisions.Mixed_Cell;

  function Convert ( s : Integer_Mixed_Subdivisions.Mixed_Subdivision )
                   return Floating_Mixed_Subdivisions.Mixed_Subdivision;
  function Convert ( s : Floating_Mixed_Subdivisions.Mixed_Subdivision )
                   return Integer_Mixed_Subdivisions.Mixed_Subdivision;

end Floating_Integer_Convertors;
