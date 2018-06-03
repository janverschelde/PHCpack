with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

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

  procedure Produce_Cells
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stable : in boolean; stlb : in double_float;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Calls DEMiCs to produce the cells, with the option
  --   to compute the stable mixed volume.

  -- ON ENTRY :
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   stable   flag to compute the stable mixed volume;
  --   stlb     the lifting bound on the artificial origins;
  --   verbose  flag for more information.

  -- ON RETURN :
  --   sup      supports with artificial origins added if stable.

  procedure Produce_Cells
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : in Standard_Floating_VecVecs.Link_to_VecVec;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Calls DEMiCs to produce the cells, for given lifting values.

  -- ON ENTRY :
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   lif      lifting values for the supports, of mix'range;
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

  procedure Pipelined_Mixed_Indices
              ( nt : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                process : access procedure
                  ( idtask : in integer32;
                    mix : in Standard_Integer_Vectors.Link_to_Vector;
                    idx : in Standard_Integer_Vectors.Vector ) := null;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Implements a 2-stage pipeline where the first task produces
  --   the cells and the other tasks consume the cells.
  --   For each cell the procedure process is called.

  -- ON ENTRY :
  --   nt       number of tasks, which must be at least 2;
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   process  optional procedure which takes the task number, the type
  --            of mixture, and the integer cell indices on input;
  --   verbose  if true, then the integer cell indices are written.

  procedure Pipelined_Mixed_Cells
              ( nt,dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : in Standard_Floating_VecVecs.Link_to_VecVec;
                lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                process : access procedure
                  ( idtask : in integer32;
                    mix : in Standard_Integer_Vectors.Link_to_Vector;
                    mic : in Mixed_Cell ) := null;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Implements a 2-stage pipeline where the first task produces
  --   the cells and the other tasks consume the cells.
  --   For each cell the procedure process is called.

  -- ON ENTRY :
  --   nt       number of tasks, which must be at least 2;
  --   dim      dimension of the points in sup, before the lifting;
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   lif      lifting values for each point in sup, of mix'range;
  --   lifsup   lifted supports of mix'range;
  --   process  optional procedure which takes the task number, the type
  --            of mixture, and the integer cell indices on input;
  --   verbose  if true, then the integer cell indices are written.

end Pipelined_Cell_Indices;
