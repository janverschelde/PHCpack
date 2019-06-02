with System;

package body Unix_Resource_Usage is

  package C_Interfaces is

    times_map : constant array (times_enum) of integer
              := (self => 0, children => -1);

    function getrusage
               ( who : integer; rusage : system.address ) return integer;
   -- pragma interface(C, getrusage); -- obsolescent feature
    pragma Import(C, getrusage);

    function timeval_to_duration ( tv : timeval ) return duration;

  end C_Interfaces;

  function Get_Process_Times
             ( who : times_enum := self ) return Process_Times is

    answer : Process_Times;
    c_result : integer;

  begin
    c_result := C_Interfaces.getrusage(who => C_Interfaces.times_map(who),
                                       rusage => answer'address);
    if (c_result = -1) then
      raise program_error;	-- something broke in Unix!
    else
      return answer;
    end if;
  end Get_Process_Times;

  function Total_Time ( t : in Process_Times ) return duration is
  begin
    return User_CPU_Time(t) + System_CPU_Time(t);
  end Total_Time;

  function User_CPU_Time ( t : in Process_Times ) return duration is
  begin
    return C_Interfaces.timeval_to_duration(t.ru_utime);
  end User_CPU_Time;

  function System_CPU_Time ( t : in Process_Times ) return duration is
  begin
    return C_Interfaces.timeval_to_duration(t.ru_stime);
  end System_CPU_Time;

  function Max_Resident_Set_Size ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_maxrss;
  end max_resident_set_size;

  function Shared_Pages_Value ( t : in Process_Times ) return page_seconds is
  begin
    return page_seconds(t.ru_ixrss);
  end Shared_Pages_Value;

  function Unshared_Data_Pages_Value
             ( t : in Process_Times ) return page_seconds is
  begin
    return page_seconds(t.ru_idrss);
  end Unshared_Data_Pages_Value;

  function Stack_Pages_Value ( t : in Process_Times ) return page_seconds is
  begin
    return page_seconds(t.ru_isrss);
  end Stack_Pages_Value;

  function Non_IO_Page_Faults ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_minflt;
  end Non_IO_Page_Faults;

  function IO_Page_Faults ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_majflt;
  end IO_Page_Faults;

  function Swaps ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_nswap;
  end Swaps;

  function Input_Blocks ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_inblock;
  end Input_Blocks;

  function Output_Blocks ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_outblock;
  end Output_Blocks;

  function Socket_Messages_Sent ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_msgsnd;
  end Socket_Messages_Sent;

  function Socket_Messages_Received
             ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_msgrcv;
  end Socket_Messages_Received;

  function Signals_Delivered ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_nsignals;
  end Signals_Delivered;

  function Voluntary_Context_Switches
             ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_nvcsw;
  end Voluntary_Context_Switches;

  function Involuntary_Context_Switches
             ( t : in Process_Times ) return long_integer is
  begin
    return t.ru_nivcsw;
  end Involuntary_Context_Switches;

  package body C_Interfaces is

    function timeval_to_duration ( tv : timeval ) return duration is

     -- package long_integer_io is new text_io.integer_io(long_integer);
     -- package integer_io is new text_io.integer_io(integer);

      answer : duration;

    begin
     -- text_io.put("seconds : "); long_integer_io.put(tv.tv_sec);
     -- text_io.new_line;
     -- text_io.put("milliseconds : "); integer_io.put(tv.tv_usec);
     -- text_io.new_line;
     -- answer := duration(tv.tv_sec); -- + duration(tv.tv_usec mod 1000)/1000;
      -- older code is below
      -- on a sun:
 answer := duration(tv.tv_sec) + duration(tv.tv_usec)/1_000_000;
      -- with the following trials, only the seconds were printed:
      -- answer := duration(tv.tv_sec + tv.tv_usec/1_000_000);
      --  answer := duration(tv.tv_sec) + duration(tv.tv_usec/1_000_000);
      return answer;
    end timeval_to_duration;

  end C_Interfaces;

end Unix_Resource_Usage;
