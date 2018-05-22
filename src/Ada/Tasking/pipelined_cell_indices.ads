with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;

package Pipelined_Cell_Indices is

-- DESCRIPTION :
--   In a pipelined processing of the cell indices,
--   one task runs DEMiCs, while the other tasks process the cells.

  procedure Produce_Cells
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Calls DEMiCs to produce the cells.

  -- ON ENTRY :
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   verbose  flag for more information.

  procedure Consume_Cells
               ( nt : in integer32;
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 process : access procedure
                   ( idtask : in integer32;
                     mix : in Standard_Integer_Vectors.Link_to_Vector;
                     idx : in Standard_Integer_Vectors.Vector ) := null;
                 verbose : in boolean := true );

  -- DESCRIPTION :
  --   The cells in DEMiCs_Output_Data are consumed by nt tasks.
  --   For each cell the procedure process is called.
  --   If verbose or if process is null, then the indices are shown on screen.
  --   Assumes DEMiCs_Output_Data contains all cells,
  --   i.e.: Produce_Cells has terminated.

  -- ON ENTRY :
  --   nt       number of tasks;
  --   mix      type of mixture;
  --   process  optional procedure which takes the task number, the type
  --            of mixture, and the integer cell indices on input;
  --   verbose  if true, then the integer cell indices are written.

end Pipelined_Cell_Indices;
