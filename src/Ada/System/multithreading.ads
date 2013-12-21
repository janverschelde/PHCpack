package Multithreading is

-- DESCRIPTION :
--   This package provides a simple interface to the GNAT.Threads.
--   Its purpose is to isolate compiler dependencies.

  generic

    with procedure Thread_Code ( id : in integer );

  procedure Start_Threads ( n : in natural );

  -- DESCRIPTION :
  --   Starts n threads.  Each thread executes the same code
  --   as defined by Thread_Code.  The parameter "id" is in the 
  --   range 1..n and is the identification number of the thread.

end Multithreading;
