with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;

package Mixed_Labels_Queue is

-- DESCRIPTION :
--   This package provides a thread safe representation of a queue:
--   getting the next element in the queue is guarded by a semaphore.

  procedure Start;

  -- DESCRIPTION :
  --   Initializes the start of the production/consumption cycle,
  --   which must be called before the launching of the tasks.

  procedure Append ( idx : in Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Appends the indices idx to the queue of labels.
  --   This procedure is called by the producer.

  procedure Stop;

  -- DESCRIPTION :
  --   This procedure is called by the producer to indicate that
  --   the production of labels has stopped.

  function Next return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   The consumer calls this function to get the next labels.

  function Next_Counter return integer32;

  -- DESCRIPTION :
  --   Returns the counter of the labels currently processed.

  procedure Clear;

  -- DESCRIPTION :
  --   Clears the stored labels from the queue.

end Mixed_Labels_Queue;
