with Ada.Calendar;
with Ada.Text_Io;
with Sigint_Handler;
with Sigint_Counter;
 
procedure Sigint_Handler_Test is

  pragma Unreserve_All_Interrupts;

  task Sig_Handler;
 
  task body Sig_Handler is

    Start_Time : constant Ada.Calendar.Time
               := Ada.Calendar.Clock;
    Sig_Time : Ada.Calendar.Time;

    use Ada.Calendar;

  begin
    Sigint_Handler.Handler.Wait;
    loop
      Sig_Time := Ada.Calendar.Clock;
      Ada.text_io.put_line(" execution took"
        & Duration'Image(Sig_Time - Start_Time) & " seconds");
      Sigint_Counter.Counter.Stop;
      exit when not Sigint_Counter.continue;
    end loop;
  end Sig_Handler;
 
begin
  null;
end Sigint_Handler_Test;
