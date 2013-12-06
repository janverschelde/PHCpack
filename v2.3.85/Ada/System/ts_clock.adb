with Ada.Calendar;                use Ada.Calendar;
with text_io,integer_io;          use text_io,integer_io;

procedure ts_clock is

-- DESCRIPTION : test on Ada.Calendar.

  this_moment : Time;
  today,this_month,this_year,this_hour,minutes,this_sec : integer;
  now_sec : Duration;

  procedure Seconds_into_HMS ( seconds : in Duration;
                               hour,min,sec : out integer ) is

    remainder : integer;

  begin
    hour := integer(seconds)/3600;
    remainder := integer(seconds)-hour*3600;
    min := remainder/60;
    remainder := remainder-min*60;
    sec := remainder/60;
  end Seconds_into_HMS;

begin
  new_line;
  put("The date today is ");
  this_moment := Clock;
  Split(this_moment,this_year,this_month,today,now_sec);
 -- today := Day(this_moment);
 -- this_month := Month(this_moment);
 -- this_year := Year(this_moment);
  Seconds_into_HMS(now_sec,this_hour,minutes,this_sec);
  put(today,2); put("/");
  put(this_month,2); put("/");
  put(this_year,4);
  put(" : "); put(this_hour,1);
  put(" : "); put(minutes,1);
  put(" : "); put(this_sec,1);
  new_line;
  new_line;
end ts_clock;
