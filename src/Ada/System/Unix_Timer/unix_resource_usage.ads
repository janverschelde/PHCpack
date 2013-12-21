package Unix_Resource_Usage is

-- DESCRIPTION :
--   Provides a wrapper around the getrusage BSD system call.
--   The getrusage commands gives information about resource utilization.
--   Previous implementations on 32-bit machines used only integer,
--   but long_integer was required for x86_64 platforms.

  type Process_Times is private;

  type times_enum is (self, children);

  type page_seconds is new long_integer;

  --   Expressed in units of pages * clock ticks (1 tick = 1/50 second).
  --   The value is calculated by summing the number of shared memory
  --   pages in use each time the internal system clock ticks, and then
  --   averaging over 1 second intervals.

  function Get_Process_Times ( who : times_enum := self ) return Process_Times;

  function Total_Time ( t : in Process_Times ) return duration;

  function User_CPU_Time ( t : in Process_Times ) return duration;

  function System_CPU_Time ( t : in Process_Times ) return duration;

  function Max_Resident_Set_Size ( t : in Process_Times ) return long_integer;

  function Shared_Pages_Value
             ( t : in Process_Times ) return page_seconds;

  -- DESCRIPTION :
  --   Returns the amount of memory used by the text segment which was
  --   also shared among other processes.

  function Unshared_Data_Pages_Value
              ( t : in Process_Times ) return page_seconds;

  -- DESCRIPTION :
  --   Returns the amount of unshared memory residing in the data segment
  --   of the process.

  function Stack_Pages_Value ( t : in Process_Times ) return page_seconds;

  -- DESCRIPTION :
  --   Returns the amount of unshared memory residing in the stack segment
  --   of the process

  function Non_IO_Page_Faults ( t : in Process_Times ) return long_integer;

  -- DESCRIPTION :
  --   Returns the number of page faults serviced without any I/O activity;
  --   here I/O activity is avoided by "reclaiming" a page frame from the
  --   list of pages awaiting reallocation.

  function IO_Page_Faults ( t : in Process_Times ) return long_integer;

  -- DESCRIPTION :
  --   Returns the number of page faults serviced which required I/O activity.

  function Swaps ( t : in Process_Times ) return long_integer;

  -- DESCRIPTION :
  --   Returns the number of times the process was swapped out of main memory.

  function Input_Blocks ( t : in Process_Times ) return long_integer;

  -- DESCRIPTION :
  --   Returns the number of times the file system had to perform input.

  function Output_Blocks ( t : in Process_Times ) return long_integer;

  -- DESCRIPTION :
  --   Returns the number of times the file system had to perform output.

  function Socket_Messages_Sent ( t : in Process_Times ) return long_integer;

  -- DESCRIPTION :
  --   Returns the number of messages sent over sockets.

  function Socket_Messages_Received
             ( t : in Process_Times ) return long_integer;

  -- DESCRIPTION :
  --   Returns the number of messages received over sockets.

  function Signals_Delivered ( t : in Process_Times ) return long_integer;

  -- DESCRIPTION :
  --   Returns the number of signals delivered.

  function Voluntary_Context_Switches
             ( t : in Process_Times ) return long_integer;

  -- DESCRIPTION :
  --   Returns the number of times a context switch resulted due to a process
  --   voluntarily giving up the processor before its time slice was completed
  --   (usually to await availability of a resource).

  function Involuntary_Context_Switches
             ( t : in Process_Times ) return long_integer;

  -- DESCRIPTION :
  --   Returns the number of times a context switch resulted due to a
  --   higher priority process becoming runnable or because the current
  --   process exceeded its time slice.

private

  type timeval is record
    tv_sec : long_integer;      -- long_integer needed for x86_64
    tv_usec : integer;          -- but no long_integer on x86_64
  end record;

  type rusage is record
    ru_utime : timeval;
    ru_stime : timeval;
    ru_maxrss : long_integer;
    ru_ixrss : long_integer;    -- integral shared text memory size
    ru_idrss : long_integer;    -- integral unshared data size
    ru_isrss : long_integer;    -- integral unshared stack size
    ru_minflt : long_integer;   -- page reclaims
    ru_majflt : long_integer;   -- page faults
    ru_nswap : long_integer;    -- swaps
    ru_inblock : long_integer;  -- block input operations
    ru_outblock : long_integer; -- block output operations
    ru_msgsnd : long_integer;   -- messages sent
    ru_msgrcv : long_integer;   -- messages received
    ru_nsignals : long_integer; -- signals received
    ru_nvcsw : long_integer;    -- voluntary context switches
    ru_nivcsw : long_integer;   -- involuntary context switches
  end record;

  type process_times is new rusage;

  pragma inline(get_process_times,total_time,user_cpu_time,
                system_cpu_time,max_resident_set_size,
                shared_pages_value,unshared_data_pages_value,
                stack_pages_value,non_io_page_faults,io_page_faults,
                swaps,input_blocks,output_blocks,socket_messages_sent,
                socket_messages_received,signals_delivered,
                voluntary_context_switches,involuntary_context_switches);

end Unix_Resource_Usage;
