with Text_Io, Unix_Resource_Usage;

package body Timing_Package is

  package duration_io is new text_io.fixed_io(duration);
  package integer_io is new text_io.integer_io (integer);

  function duration_to_string (dur : duration)
      return string;
  pragma inline (duration_to_string);

  type timing_item is record
    start_time : Unix_Resource_Usage.Process_Times;
    stop_time : Unix_Resource_Usage.Process_Times;
  end record;

  type rusage_timing_stuff is
    record
      total_time : duration;
      user_time : duration;
      system_time : duration;
--      max_resident_size : natural;
--      shared_pages : natural;
--      unshared_pages : natural;
--      stack_pages : natural;
      non_io_faults : natural;
      io_faults : natural;
      swaps : natural;
--      input_blocks : natural;
--      output_blocks : natural;
--      messages_out : natural;
--      messages_in : natural;
      signals : natural;
--      vol_context_switches : natural;
--      invol_context_switches : natural;
      total_context_switches : natural;
  end record;

  function to_rusage_timing_stuff (item : Timing_Widget)
      return rusage_timing_stuff;
  pragma inline (to_rusage_timing_stuff);

  procedure tstart ( widget : out Timing_Widget ) is

    answer : timing_item;

  begin
    answer.start_time := Unix_Resource_Usage.get_process_times;
    widget := new timing_item'(answer);
  end tstart;

  procedure tstop ( widget : in out Timing_Widget ) is
  begin
    widget.all.stop_time := Unix_Resource_Usage.get_process_times;
  end tstop;

  function Elapsed_Total_Time ( widget : Timing_Widget ) return duration is
  begin
    return (Unix_Resource_Usage.Total_Time_of(widget.Stop_Time)
            - Unix_Resource_Usage.Total_Time_of(widget.Start_Time));
  end Elapsed_Total_Time;

  function Elapsed_User_Time ( widget : Timing_Widget ) return duration is
  begin
    return (Unix_Resource_Usage.User_CPU_Time_of(widget.Stop_Time)
            - Unix_Resource_Usage.User_CPU_Time_of(widget.Start_Time));
  end Elapsed_User_Time;

  function Elapsed_System_Time ( widget : Timing_Widget ) return duration is
  begin
    return (Unix_Resource_Usage.System_CPU_Time_of(widget.Stop_Time)
            - Unix_Resource_Usage.System_CPU_Time_of(widget.Start_Time));
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
                          widget : Timing_Widget; tag : string := "") is

    rusage_stuff : rusage_timing_stuff;

    printout_column : text_io.positive_count := 40;
    function "+" (l,r : text_io.positive_count)
        return text_io.positive_count   renames Text_IO."+";

  begin
    text_io.put_line (file,"TIMING INFORMATION for " & tag);
    rusage_stuff := to_rusage_timing_stuff (widget);
    -- print out total time
    Text_Io.Put (file,"The elapsed time in seconds was ");
   -- text_io.set_col (file,Text_IO.positive_count(printout_column + 3));
   -- sun4 modification :
    text_io.set_col (file,Text_IO.positive_count(printout_column));
    duration_io.put (file,rusage_stuff.total_time);
    text_io.put(file," = "); print_hms(file,rusage_stuff.total_time);
    text_io.new_line(file);
    -- print out user time
    Text_Io.Put (file,"User time in seconds was ");
   -- text_io.set_col (file,Text_IO.positive_count(printout_column + 3));
   -- sun4 modification :
    text_io.set_col (file,Text_IO.positive_count(printout_column));
    duration_io.put (file,rusage_stuff.user_time);
    text_io.put(file," = "); print_hms(file,rusage_stuff.user_time);
    text_io.new_line(file);
    -- print out system time
    Text_Io.Put (file,"System CPU time in seconds was ");
   -- text_io.set_col (file,text_io.positive_count(printout_column + 3));
   -- sun4 modification :
    text_io.set_col (file,text_io.positive_count(printout_column));
    duration_io.put (file,rusage_stuff.system_time);
    text_io.put(file," = "); print_hms(file,rusage_stuff.system_time);
    text_io.new_line(file);
    -- print out non-I/O page faults
    Text_IO.put (file,"Non-I/O page faults was ");
    text_io.set_col (file,printout_column);
    integer_io.put  (file,rusage_stuff.non_io_faults);
    Text_IO.new_line(file);
    -- print out I/O page faults
    Text_IO.put (file,"I/O page faults was ");
    text_io.set_col (file,printout_column);
    integer_io.put  (file,rusage_stuff.io_faults);
    text_io.new_line(file);
    -- print out signals
    Text_IO.put (file,"Signals delivered was ");
    text_io.set_col (file,printout_column);
    integer_io.put  (file,rusage_stuff.signals);
    text_io.new_line(file);
    -- print out swaps
    text_io.put (file,"Swaps was ");
    text_io.set_col (file,printout_column);
    integer_io.put  (file,rusage_stuff.swaps);
    text_io.new_line(file);
    -- print out total context switches
    text_io.put (file,"Total context switches was ");
    text_io.set_col (file,printout_column);
    integer_io.put  (file,rusage_stuff.total_context_switches);
    text_io.new_line(file);
 --   text_io.put_line
 --("-----------------------------------------------------------------");
  end print_times;

  function times_to_string (widget : Timing_Widget;
			    delimiter : string := ":")
      return string
  is
    rusage_stuff : rusage_timing_stuff;
  begin
    rusage_stuff := to_rusage_timing_stuff(widget);
    return   "Total Time in seconds  => "
	   & duration_to_string (rusage_stuff.total_time) & delimiter
	   & "User Time in seconds   => "
	   & duration_to_string (rusage_stuff.user_time) & delimiter
	   & "System Time in seconds => "
	   & duration_to_string (rusage_stuff.system_time) & delimiter
	   & "Non I/O Page Faults    =>       "
	   & integer'image (rusage_stuff.non_io_faults) & delimiter
	   & "I/O Page Faults        =>       "
	   & integer'image (rusage_stuff.io_faults) & delimiter
	   & "Swaps                  =>       "
	   & integer'image (rusage_stuff.swaps) & delimiter
	   & "Signals Delivered      =>       "
	   & integer'image (rusage_stuff.signals) & delimiter
	   & "Total Context Switches =>       "
	   & integer'image (rusage_stuff.total_context_switches) & delimiter;
  end times_to_string;

  function duration_to_string (dur : duration) return string is

    answer : string(1..(duration'fore + duration'aft + 1));

  begin
    duration_io.put (to => answer, item => dur);
    return answer;
  end duration_to_string;

  function to_rusage_timing_stuff (item : Timing_Widget)
                                  return rusage_timing_stuff is

    answer : rusage_timing_stuff;

  begin
    answer.total_time := Unix_Resource_Usage.total_time_of
					(item.Stop_Time)
	   	         - Unix_Resource_Usage.total_time_of
					(item.Start_Time);
    answer.user_time := Unix_Resource_Usage.user_cpu_time_of
					(Item.Stop_time)
		        - Unix_Resource_Usage.user_cpu_time_of
					(Item.Start_time);
    answer.system_time := Unix_Resource_Usage.system_cpu_time_of
					(Item.Stop_time)
		          - Unix_Resource_Usage.system_cpu_time_of
					(Item.Start_time);
    answer.non_io_faults := integer(Unix_Resource_Usage.non_io_page_faults_of
		 			(item.stop_time)
		            	    - Unix_Resource_Usage.non_io_page_faults_of
					(item.start_time));
    answer.io_faults := integer (Unix_Resource_Usage.io_page_faults_of
		 			(item.stop_time)
		                 - Unix_Resource_Usage.io_page_faults_of
					(item.start_time));
    answer.swaps := (Unix_Resource_Usage.swaps_of
		 			(item.stop_time)
		       - Unix_Resource_Usage.swaps_of
					(item.start_time));
    answer.signals := (Unix_Resource_Usage.signals_delivered_of
		 			(item.stop_time)
		       - Unix_Resource_Usage.signals_delivered_of
					(item.start_time));

    answer.total_context_switches
	:= (   Unix_Resource_Usage.voluntary_context_switches_of
    				(item.stop_time)
	       - Unix_Resource_Usage.voluntary_context_switches_of
			 	(item.start_time))
	    + (Unix_Resource_Usage.involuntary_context_switches_of
    				(item.stop_time)
	       - Unix_Resource_Usage.involuntary_context_switches_of
			 	(item.start_time));
    return answer;
  end to_rusage_timing_stuff;

end Timing_Package;
