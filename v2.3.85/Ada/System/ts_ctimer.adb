with text_io,integer_io;       use text_io,integer_io;
with Timing_Package;           use Timing_Package;

procedure ts_ctimer is

  function get_clock return integer;
  pragma Import(C,get_clock,"get_clock");

  function get_clocks_per_sec return integer;
  pragma Import(C,get_clocks_per_sec,"get_clocks_per_sec");

  clocks_per_sec : constant integer := get_clocks_per_sec;

  n : constant integer := 1500000000;
  timer : Timing_Widget;
  a,b,c : float;
  save_clock,now,user_time,millisec,elapsed : integer;

begin
  put("Clocks per second : "); put(clocks_per_sec,1); new_line;
  b := 2.333; c := 3.4444;
  tstart(timer);
  save_clock := get_clock;
  for i in 1..n loop
    a := b*c;
    if i mod 10000 = 0
     then put_line("10000 multiplications done");
    end if;
  end loop;
  now := get_clock;
  tstop(timer);
  elapsed := now - save_clock;
  user_time := elapsed/clocks_per_sec;
  millisec := elapsed/(clocks_per_sec/1000);
  put("Elapsed user time : ");
  put(user_time,1); put_line(" seconds.");
  put("Elapsed user time : ");
  put(millisec,1); put_line(" milliseconds.");
  millisec := millisec - user_time*1000;
  put("Elapsed user time : ");
  put(user_time,1); put(" seconds and ");
  put(millisec,1); put_line(" milliseconds.");
  print_times(Standard_Output,timer,"Checking C timer");
end ts_ctimer;
