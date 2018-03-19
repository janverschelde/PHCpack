with Ada.Calendar;                use Ada.Calendar;
with text_io,integer_io;          use text_io,integer_io;
with Time_Stamps;                 use Time_Stamps;

procedure ts_stamps is

-- DESCRIPTION : test on Ada.Calendar.

  procedure Main is

  -- DESCRIPTION :
  --   Test on measuring the elapsed wall clock time.

    start_moment,end_moment : Time;
    ans : character;

  begin
    start_moment := Clock;
    new_line;
    put("Session starts at ");
    Write_Time_Stamp(Standard_Output,start_moment);
    put_line(".");
    new_line;
    loop
      put("Have you waited long enough ? (y/n) "); get(ans);
      exit when ans = 'y';
    end loop;
    new_line;
    end_moment := Clock;
    put("Session ends at ");
    Write_Time_Stamp(Standard_Output,end_moment);
    put_line(".");
    put("The elapsed time is ");
    put(Elapsed_Time(start_moment,end_moment),1);
    put_line(" seconds.");
    new_line;
  end Main;

begin
  Main;
end ts_stamps;
