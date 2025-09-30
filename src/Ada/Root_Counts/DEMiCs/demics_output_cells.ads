with Standard_Natural_Numbers;        use Standard_Natural_Numbers;
with Standard_Integer_Numbers;        use Standard_Integer_Numbers;
with Standard_Floating_Numbers;       use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Lists_of_Integer_Vectors;
with Floating_Mixed_Subdivisions;     use Floating_Mixed_Subdivisions;

package DEMiCs_Output_Cells is

-- DESCRIPTION :
--   This package stores the output cells computed by DEMiCs
--   and arranges the conversion to mixed cells with floating liftings.
--   Its interface is very similar to DEMiCs_Output_Data when calling
--   the C++ code directly, with the exception that the indices to the
--   cells are no longer stored as strings, but as index vectors.

  mixed_volume : integer32;    -- use -1 to indicate an error occurred
  monitor : boolean := false;  -- write cell indices with every add
  stable : boolean := false;   -- if stable mixed volume wanted
  stlb : double_float := 0.0;  -- stable lifting bound if stable
  done : boolean := false;     -- have all cell indices been computed?
  allocate : boolean := false; -- allocate memory for the mixed cells?

-- CONSTRUCTORS :

  procedure Initialize_Lifting
              ( crdsup : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Allocates memory for the lifting values,
  --   given in crdsup(i) the number of points in the i-th support.

  procedure Assign_Lifting 
              ( idxsup,idxpnt : in integer32; val : in double_float );

  -- DESCRIPTION :
  --   Assigns the point with index idxpnt in support idxsup
  --   the lifting value equal to val.

  -- REQUIRED : Initialized_Lifting() was executed and both idxsup
  --   and idxpnt are in valid ranges.

  procedure Store_Dimension_and_Mixture
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Stores the dimension dim and the type of mixture,
  --   which is necessary if allocate is set to true.

  procedure Allocate_Mixed_Cell;

  -- DESCRIPTION :
  --   Allocates memory for one new mixed cell.
  --   This procedure is called by Add_Cell_Indices if allocate.

  -- REQUIRED : Store_Dimension_and_Mixture has been executed.

  procedure Add_Cell_Indices
              ( idx : in Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Adds the cell indices in idx to the list.
  --   If the flag monitor is true,
  --   then the indices are written to screen.
  --   If allocate, then memory is allocated for the mixed cell
  --   corresponding to the given string representation.

  procedure Initialize_Cell_Indices_Pointer;

  -- DESCRIPTION :
  --   Sets the cell pointer to the beginning of the stored list
  --   of cell indices.

  procedure Initialize_Allocated_Cell_Pointer;

  -- DESCRIPTION :
  --   Initializes the internal pointer to the beginning
  --   of the allocated list of mixed cells.

-- SELECTORS :

  function Get_Labels_Size return integer32;

  -- DESCRIPTION :
  --   Returns the size of the indices vector of the labels
  --   for the mixed cells, computed when the mixture is stored.

  function Get_Mixture return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the stored mixture.

  function Retrieve_Lifting
             ( idxsup,idxpnt : integer32 ) return double_float;

  -- DESCRIPTION :
  --   Returns the lifting value of the point with index idxpnt
  --   in the support with index idxsup.

  -- REQUIRED : Initialized_Lifting() was executed and both idxsup
  --   and idxpnt are in valid ranges.

  function Lifting_Values return Standard_Floating_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns the stored lifting values.

  function Number_of_Cell_Indices return natural32;

  -- DESCRIPTION :
  --   Returns the number of cell indices stored.

  function Get_Cell_Indices
             ( index : integer32 )
             return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Given the index to indices in a list,
  --   returns the pointer to the vector at that position,
  --   or simply null if the index is out of bounds.
  --   Observe that the cost of this operation is proportional
  --   to the size of the accumulated list.
  --   The Get_Next_Cell_Indices is more efficient if all cell indices
  --   need to be retrieved.

  function Retrieve_Cell_Indices return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Returns list of all cell indices accumulated
  --   by the application of the Add_Cell_Indices procedure.

  function Get_Next_Cell_Indices
             return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns a pointer to the next cell indices.
  --   If there is no next cell, then null is returned,
  --   otherwise the pointer is moved the next cell.

  function Get_Allocated_Cells return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Returns the list of all allocated cells.

  function Get_Next_Allocated_Cell return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Returns the pointer to the next allocated cell
  --   and moves the pointer to the next cell if not empty.

-- DESTRUCTORS :

  procedure Clear_Lifting;

  -- DESCRIPTION :
  --   Deallocates the memory for the lifting values.

  procedure Clear_Cell_Indices;

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the cell indices.

  procedure Clear_Allocated_Cells;

  -- DESCRIPTION :
  --   Deallocates the memory occupied by all allocated cells.

end DEMiCs_Output_Cells;
