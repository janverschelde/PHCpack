generic

  with procedure do_job ( idn : in integer );

  -- DESCRIPTION :
  --   Executed by worker with identification number idn.

procedure Generic_Task_Array ( p : in integer );

-- DESCRIPTION :
--   Defines an array of p workers,
--   every worker executes the do_job procedure.
