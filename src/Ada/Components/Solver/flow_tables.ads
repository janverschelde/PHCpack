with text_io;                          use text_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Natural_Matrices;        use Standard_Natural_Matrices;
with Standard_Complex_Solutions;       use Standard_Complex_Solutions;

package Flow_Tables is

-- DESCRIPTION :
--   A flow table stores the transitions of the candidate witness sets
--   as computed by the equation-by-equation solver.

  type Flow_Table is private;

  function Create ( nb_eqs,nb_var : integer32 ) return Flow_Table;

  -- DESCRIPTION :
  --   Creates an empty flow table for a system with #equations = nb_eqs 
  --   and #variables = nb_var.

  procedure Store_Hypersurface_Degree 
              ( ft : in out Flow_Table; k,d : in integer32 );

  -- DESCRIPTION :
  --   Stores in the flow table the degree of the k-th hypersurface.

  procedure Update ( ft : in out Flow_Table; k : in integer32;
                     s : in Array_of_Solution_Lists );

  -- DESCRIPTION :
  --   Updates the flow table with the array of witness sets.

  procedure Write ( file : in file_type; ft : in Flow_Table );

  -- DESCRIPTION :
  --   Writes the flow table to the file.

  procedure Clear ( ft : in out Flow_Table );

  -- DESCRIPTION :
  --   Deallocation of all memory.

private

  type Flow_Table is new Link_to_Matrix;

end Flow_Tables;
