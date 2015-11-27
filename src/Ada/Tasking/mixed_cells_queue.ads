with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Mixed_Cells_Queue is

-- DESCRIPTION :
--   This package provides a thread safe representation of a queue:
--   getting the next element in the queue is guarded by a semaphore.

  procedure Initialize ( cells : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Initializes the internal pointer to the list cells.

  function Next return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Returns the pointer to the list of cells that has as head
  --   element the next element in the list and then moves this
  --   pointer to the next element of the list.
  --   After executing Next as many times as Length_Of(cells),
  --   where cells was the argument of Initialize, the pointer
  --   returned will be null.

  -- REQUIRED :
  --   The procedure Initialize has been called.

  function Next_Counter return integer32;

  -- DESCRIPTION :
  --   Returns the counter of the cell currently processed.
  --   This function returns the number of times Next was called.

end Mixed_Cells_Queue;
