with Ada.Calendar;                use Ada.Calendar;
with text_io,integer_io;          use text_io,integer_io;

procedure ts_clock is

  package duration_io is new text_io.fixed_io(duration);

-- DESCRIPTION : test on Ada.Calendar.

  procedure Seconds_into_HMSMS
              ( seconds : in Duration;
                hour,min,sec,millisec : out integer ) is

  -- DESCRIPTION :
  --   Splits the seconds into hours, minutes, seconds,
  --   and milliseconds.

    remainder : integer;
    wholeseconds : constant Duration := duration(integer(seconds));
    milliseconds : constant Duration := seconds - wholeseconds;

  begin
    hour := integer(seconds)/3600;
    remainder := integer(seconds)-hour*3600;
    min := remainder/60;
    remainder := remainder-min*60;
    sec := integer(remainder);
    millisec := integer(milliseconds*1000.0);
  end Seconds_into_HMSMS;

  procedure Main is

    this_moment : Time;
    today,this_month,this_year : integer;
    this_hour,minutes,this_sec,this_millisec : integer;
    now_sec : Duration;

  begin
    new_line;
    put("The date today is ");
    this_moment := Clock;
    Split(this_moment,this_year,this_month,today,now_sec);
    put("the now_sec : "); duration_io.put(now_sec); new_line;
    Seconds_into_HMSMS(now_sec,this_hour,minutes,this_sec,this_millisec);
    put(today,2); put("/");
    put(this_month,2); put("/");
    put(this_year,4);
    put(" : "); put(this_hour,1);
    put(" : "); put(minutes,1);
    put(" : "); put(this_sec,1);
    put(" : "); put(this_millisec,1);
    new_line;
    new_line;
  end Main;

begin
  Main;
end ts_clock;
