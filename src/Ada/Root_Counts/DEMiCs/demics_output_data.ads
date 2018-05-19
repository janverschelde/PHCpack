with Standard_Natural_Numbers;        use Standard_Natural_Numbers;
with Standard_Integer_Numbers;        use Standard_Integer_Numbers;
with Standard_Floating_Numbers;       use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with String_Splitters;
with Lists_of_Strings;

package DEMiCs_Output_Data is

-- DESCRIPTION :
--   This package stores the output data computed by DEMiCs.

  mixed_volume : integer32; -- use -1 to indicate an error occurred

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

  procedure Add_Cell_Indices ( strcell : in string );

  -- DESCRIPTION :
  --   Adds the string representation of a cell,
  --   as defined by its indices.

  procedure Initialize_Cell_Pointer;

  -- DESCRIPTION :
  --   Sets the cell pointer to the beginning of the stored list
  --   of cell indices.

-- SELECTORS :

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
             ( index : integer32 ) return String_Splitters.Link_to_String;

  -- DESCRIPTION :
  --   Given the index to a string in a list,
  --   returns the pointer to the string at that position,
  --   or simply null if the index is out of bounds.

  function Retrieve_Cell_Indices return Lists_of_Strings.List;

  -- DESCRIPTION :
  --   Returns list of the cell indices accumulated
  --   by the application of the Add_Cell_Indices procedure.

  function Get_Next_Cell return String_Splitters.Link_to_String;

  -- DESCRIPTION :
  --   Returns a pointer to the next cell indices.
  --   If there is no next cell, then null is returned,
  --   otherwise the pointer is moved the next cell.

-- DESTRUCTORS :

  procedure Clear_Lifting;

  -- DESCRIPTION :
  --   Deallocates the memory for the lifting values.

  procedure Clear_Cell_Indices;

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the cell indices.

end DEMiCs_Output_Data;
