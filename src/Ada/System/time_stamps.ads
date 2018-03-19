with Ada.Calendar;                    use Ada.Calendar;
with text_io;                         use text_io;

package Time_Stamps is

-- DESCRIPTION :
--   This package uses the clock time to deliver time stamps and
--   to measure the elapsed time.

  procedure Seconds_into_HMS
              ( seconds : in Duration; hour,min,sec : out integer );

  -- DESCRIPTION :
  --   Splits the seconds into hours, minutes, and seconds.

  procedure Seconds_into_HMSMS
              ( seconds : in Duration;
                hour,min,sec,millisec : out integer );

  -- DESCRIPTION :
  --   Splits the seconds into hours, minutes, seconds,
  --   and milliseconds.

  function Elapsed_Time ( before,after : Time ) return integer;

  -- DESCRIPTION :
  --   Returns the elapsed time between before and after in seconds.

  procedure Write_Time_Stamp ( file : in file_type; moment : in Time );

  -- DESCRIPTION :
  --   Puts the moment in the format like "7 October 1999, 15:40:28".

  procedure Write_Elapsed_Time
               ( file : in file_type; before,after : in Time );

  -- DESCRIPTION :
  --   Writes the total elapsed time on file, in total seconds,
  --   followed (where appropriate) by a format in hours, minutes, seconds.

end Time_Stamps;
