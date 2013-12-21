with text_io;                           use text_io;
with Ada.Real_Time;                     use Ada.Real_Time;

procedure ts_time is

  nowtime : constant Time := Clock;
  sc : Seconds_Count;
  ts : Time_Span;
 -- now : Duration;

  package sec_io is new text_io.integer_io(Seconds_Count);

begin
  Split(nowtime,sc,ts);
 -- now := To_Duration(ts);
  put("seconds count : "); sec_io.put(sc); new_line;
end ts_time;
