with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;

package QuadDobl_Solutions_Queue is

-- DESCRIPTION :
--   This package provides a thread safe representation of a queue:
--   getting the next element in the queue is guarded by a semaphore.

  procedure Initialize ( sols : in Solution_List );

  -- DESCRIPTION :
  --   Initializes the internal pointer to the list sols.

  function Next return Solution_List;

  -- DESCRIPTION :
  --   Returns the pointer to the solution list that has as head
  --   element the next element in the list and then moves this
  --   pointer to the next element of the list.
  --   After executing Next as many times as Length_Of(sols),
  --   where sols was the argument of Initialize, the pointer
  --   returned will be null.

  -- REQUIRED :
  --   The procedure Initialize has been called.

  function Next_Counter return integer32;

  -- DESCRIPTION :
  --   Returns the counter of the solution currently processed.
  --   This function returns the number of times Next was called.

end QuadDobl_Solutions_Queue;
