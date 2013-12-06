with text_io;                      use text_io;
with Deformation_Posets;           use Deformation_Posets;

package Deformation_Posets_io is

-- DESCRIPTION :
--   Provides output routines for deformation posets.

  procedure put_size ( poset : in Array_of_Array_of_VecMats );
  procedure put_size
              ( file : in file_type; poset : in Array_of_Array_of_VecMats );

  -- DESCRIPTION :
  --   Writes the length of every vector of matrices, for every level.

end Deformation_Posets_io;
