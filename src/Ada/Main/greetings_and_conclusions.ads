with Ada.Calendar;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package Greetings_and_Conclusions is

-- DESCRIPTION :
--   This package defines the writing of the greetings and the conclusions,
--   as they appear at the beginning and end of phc -B,
--   called in the procedures compsolve, compsolve2, and compsolve4.

  procedure Write_Greeting ( nbtasks,precision : in natural32 );

  -- DESCRIPTION :
  --   Writes the greeting and the number of tasks in nbtasks,
  --   or writes no tasking, followed by the precision.
  --   The values for precision are 0 for standard double,
  --   1 for double double, and 2 for quad double.

  procedure Write_Conclusion
              ( start_moment : in Ada.Calendar.Time; nbtasks : in natural32 );

  -- DESCRIPTION :
  --   Writes the time stamp, number of tasks, seed, and version number
  --   to screen.  The start_moment is the start time of the computations.

  procedure Write_Conclusion
              ( file : in file_type;
                start_moment : in Ada.Calendar.Time; nbtasks : in natural32 );

  -- DESCRIPTION :
  --   Writes the time stamp, number of tasks, seed, and version number
  --   to file.  The start_moment is the start time of the computations.

end Greetings_and_Conclusions;
