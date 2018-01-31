with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;

package Static_Columns_Queue is

-- DESCRIPTION :
--   The package gives a thread safe representation of a queue of indices.
--   Each vector of indices in the queue represents the selected columns 
--   for the first number of rows, for a parallel permanent computation.

  procedure Initialize ( stc : in Standard_Integer_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Initializes the queue with a sequence of start columns,
  --   given in the vector of vectors stc.

  function Next_Columns return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the next vector of column indices for the first rows.
  --   If all columns have been processed, then null is returned.

  function Next_Counter return integer32;

  -- DESCRIPTION :
  --   Returns the counter of the column indices currently processed.
  --   This function returns the number of times Next_Columns was called.

  procedure Clear;

  -- DESCRIPTION :
  --   Releases the pointer allocated to store the stc at initialization.

end Static_Columns_Queue;
