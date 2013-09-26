with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

package Unfolding_Subdivisions is

-- DESCRIPTION :
--   This package contains routines to unfold subdivisions.

  function Different_Normals ( mixsub : Mixed_Subdivision ) return List;

  -- DESCRIPTION :
  --   Returns the list of all different normals of the cells in mixsub.

  function Extract ( normal : vector; mixsub : Mixed_Subdivision )
                   return Mixed_Subdivision;
  
  -- DESCRIPTION :
  --   Returns a list of cells with the given normal.

  function Merge_Same_Normal ( mixsub : Mixed_Subdivision ) return Mixed_Cell;
  function Merge_Same_Normal ( mixsub : Mixed_Subdivision )
                             return Mixed_Subdivision;

  -- DESCRIPTION :
  --   All cells with the same inner normal will be put in one cell,
  --   that will be contained in the mixed subdivision on return.
  --   The refinement of the cells will be discarded.

  -- REQUIRED :
  --   not Is_Null(mixsub) and all mixed cells have the same inner normal.

  function Merge ( mixsub : Mixed_Subdivision ) return Mixed_Subdivision;

  -- DESCRIPTION :
  --   All cells with the same inner normal will be put in one cell.
  --   For cells whose inner normal occurs more than once, the refinement
  --   will be discarded.

  function Relift ( mixsub : Mixed_Subdivision; point : Vector )
                  return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Returns a new mixed subdivision, derived from the original one,
  --   where all points different from the given point will be given
  --   lifting value zero.  The given point will receive lifting value 1.

  generic

    with procedure Process ( mic : in Mixed_Cell; newpts : in Array_of_Lists );
    -- DESCRIPTION :
    --   Returns the new mixed cell, with the new re-lifted points.

  procedure Unfolding ( mixsub : in out Mixed_Subdivision );

  -- DESCRIPTION :
  --   A collection of cells, all with the same normal, will be unfolded.

end Unfolding_Subdivisions;
