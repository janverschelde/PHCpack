with text_io;                            use text_io;

package Timing_Package is

-- DESCRIPTION :
--   This package allows to perform timings.

-- ACKNOWLEGMENT :
--   Originally developed by Dave Emery (emery@aries.mitre.org),
--   but modified by me.

  type Array_of_Duration is array ( integer range <> ) of duration;

  type Timing_Widget is private;

-- OPERATIONS :

  procedure tstart ( widget : out Timing_Widget );

  -- DESCRIPTION : Starts the timing.

  procedure tstop  ( widget : in out Timing_Widget );

  -- DESCRIPTION : Stops the timing.

  function Elapsed_Total_Time  ( widget : Timing_Widget ) return duration;
  function Elapsed_User_Time   ( widget : Timing_Widget ) return duration;
  function Elapsed_System_Time ( widget : Timing_Widget ) return duration;

  -- DESCRIPTION : Returns elapsed time.

  procedure print_time ( file : file_type; mach_time : duration );
  
  -- DESCRIPTION :
  --   Writes the duration in the standard way.

  procedure print_hms ( file : file_type; mach_time : duration );

  -- DESCRIPTION :
  --   Writes the duration in a hours/minutes/seconds output format.

  procedure print_times ( widget : Timing_Widget; tag : string := "" );

  -- DESCRIPTION :
  --   Prints as much information as is available on standard output.

  procedure print_times ( file : file_type;
                          widget : Timing_Widget; tag : string := "" );

  -- DESCRIPTION :
  --   Prints as much information as is available on a file.

  function times_to_string
              ( widget : Timing_Widget; delimiter : string := ":" )
              return string;

  -- DESCRIPTION :
  --   Returns a string with information.

private

  type timing_item;
  type Timing_Widget is access timing_item;

end Timing_Package;
