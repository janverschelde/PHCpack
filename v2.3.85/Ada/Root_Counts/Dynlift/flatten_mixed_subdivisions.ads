with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

package Flatten_Mixed_Subdivisions is

-- DESCRIPTION :
--   This package contains the flattening routines.

  procedure Flatten ( L : in out List );
  procedure Flatten ( L : in out Array_of_Lists );

  procedure Flatten ( mic : in out Mixed_Cell );
  procedure Flatten ( mixsub : in out Mixed_Subdivision );

end Flatten_Mixed_Subdivisions;

