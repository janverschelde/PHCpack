with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Complex_VecVecs;            use Standard_Complex_VecVecs;
with Standard_Complex_Poly_Systems;       use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;       use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;          use Standard_Complex_Solutions;

package Standard_Sampling_Operations is

-- DESCRIPTION :
--   This is the analogue to PHCpack_Operations to guide the path
--   tracking needed in the monodromy breakup algorithm,
--   in standard double precision.

  procedure Initialize ( p : in Poly_Sys; sols : in Solution_List;
                         k : in integer32 );
  procedure Initialize ( p : in Laur_Sys; sols : in Solution_List;
                         k : in integer32 );

  -- DESCRIPTION :
  --   Initializes the sampling machine with an embedded system p,
  --   and points in the list sols on a k-dimensional solution set.

  procedure Initialize_Slices ( n : in integer32 );

  -- DESCRIPTION :
  --   Allocates room for n more different slicing hyperplane sections.
  --   The slices for the original witness set define the 0-th slices.

  function Retrieve_First_Solutions return Solution_List;

  -- DESCRIPTION :
  --   Returns a pointer to the solutions, used with Initialize.

  function Retrieve_Start_Slices return Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns a pointer to the current set of start slices.

  function Retrieve_Target_Slices return Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns a pointer to the current set of target slices.

  function Retrieve_Slices ( i : integer32 ) return Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns a pointer to the i-th set of hyperplane sections,
  --   which is null if the i-th set does not (yet) has a value.
  --   When i = 0, the original set of slices is returned.

  procedure Add_Slices ( v : in VecVec );

  -- DESCRIPTION :
  --   Adds a copy of the given slices to the internal data structures.

  procedure Set_Target_Slices ( i : integer32 );

  -- DESCRIPTION :
  --   Sets the target hyperplanes to the i-th set of hyperplane sections,
  --   previously added with "Add_Slices".

  -- REQUIRED : i is in the range of slices that have been added.
 
  procedure Assign_Slice ( c : in Complex_Number; i,j : in integer32 );

  -- DESCRIPTION :
  --   Assigns c to the j-th coefficient of the i-th slice in the new
  --   target hyperplane sections.
  --   Notice that i starts to run from 0.

  procedure Store_Gamma ( c : in Complex_Number; i : in integer32 );

  -- DESCRIPTION :
  --   Stores c as the i-th gamma constant.
  --   Notice that i starts to run from 0.

  procedure Swap_Slices;

  -- DESCRIPTION :
  --   Swaps the start with the target set of slices,
  --   and the start solutions with the target solutions.

  procedure Sample;

  -- DESCRIPTION :
  --   Computes a new witness set on the new slices.
  --   The solutions are put in the solutions container.

  function Sample_Loop ( start,target : integer32;
                         s : Link_to_Solution ) return Link_to_Solution;

  -- DESCRIPTION :
  --   Completes one loop for one solution,
  --   starting at the slice start and ending at the slice target.

  procedure Clear;

  -- DESCRIPTION :
  --   All internal data structures managed by this package are destroyed.

end Standard_Sampling_Operations;
