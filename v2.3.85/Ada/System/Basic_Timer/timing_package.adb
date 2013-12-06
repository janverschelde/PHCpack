with integer_io; use integer_io;

package body Timing_Package is

  package duration_io is new text_io.fixed_io(duration);

-- INTERFACE WITH C :

  function get_clock return integer;
  pragma Import(C,get_clock,"get_clock");

  function get_clocks_per_sec return integer;
  pragma Import(C,get_clocks_per_sec,"get_clocks_per_sec");

  clocks_per_sec : constant integer := get_clocks_per_sec;

-- DATA STRUCTURE :

  type timing_item is record
    start_time,stop_time : integer;
  end record;

-- OPERATIONS :

  procedure tstart ( widget : out Timing_Widget ) is

    item : timing_item;
    
  begin
    item.start_time := get_clock;
    widget := new timing_item'(item);
  end tstart;

  procedure tstop  ( widget : in out Timing_Widget ) is
  begin
    widget.stop_time := get_clock;
  end tstop;

  function Elapsed_Total_Time  ( widget : Timing_Widget ) return duration is

    res : duration;
    elapsed : integer := widget.stop_time - widget.start_time;

  begin
    res := duration(elapsed/(clocks_per_sec/1000));  -- also get milliseconds!
    res := res/1000.0;
    return res;
  end Elapsed_Total_Time;

  function Elapsed_User_Time   ( widget : Timing_Widget ) return duration is
  begin
    return Elapsed_Total_Time(widget);
  end Elapsed_User_Time;

  function Elapsed_System_Time ( widget : Timing_Widget ) return duration is
  begin
    return 0.0;
  end Elapsed_System_Time;

  procedure print_time ( file : file_type; mach_time : duration ) is
  begin
    duration_io.put(file,mach_time);
  end print_time;

  function truncate ( d : duration ) return integer is

    rd : integer := integer(d);

  begin
    if rd > 0
     then if duration(rd) > d
           then rd := rd-1;
          end if;
    end if;
    return rd;
  end truncate;
  
  procedure print_hms ( file : file_type; mach_time : duration ) is

    seconds : integer := truncate(mach_time);
    millsec : integer := integer((mach_time-duration(seconds))*1000);
    minutes,hours : integer;

  begin
    if millsec >= 1000                    -- could be due to rounding
     then seconds := seconds + 1;
          millsec := millsec - 1000;
    end if;
    minutes := seconds/60;
    hours := minutes/60;
    seconds := seconds - 60*minutes;
    minutes := minutes - 60*hours;
    integer_io.put(file,hours,2);   text_io.put(file,"h");
    integer_io.put(file,minutes,2); text_io.put(file,"m");
    integer_io.put(file,seconds,2); text_io.put(file,"s");
    integer_io.put(file,millsec,3); text_io.put(file,"ms");
  end print_hms;

  procedure print_times ( widget : Timing_Widget; tag : string := "" ) is
  begin
    print_times(Standard_Output,widget,tag);
  end print_times;

  procedure print_times ( file : file_type;
                          widget : Timing_Widget; tag : string := "" ) is

    elapsed : duration;

  begin
    text_io.put_line(file,"TIMING INFORMATION for " & tag);
    text_Io.put(file,"The elapsed time in seconds was ");
    elapsed := Elapsed_Total_Time(widget);
    duration_io.put(file,elapsed);
    text_io.put(file," = "); print_hms(file,elapsed);
    text_io.new_line(file);
  end print_times;

  function times_to_string ( widget : Timing_Widget; delimiter : string := ":" )
                           return string is
  begin
    return "";
  end times_to_string;

end Timing_Package;
