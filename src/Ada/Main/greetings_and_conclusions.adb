with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Time_Stamps;
with Write_Seed_Number;
with Write_Number_of_Tasks;
with Greeting_Banners;

package body Greetings_and_Conclusions is

  procedure Write_Greeting ( nbtasks,precision : in natural32 ) is
  begin
    put_line(Greeting_Banners.welcome & ".");
    put("Numerical irreducible decomposition");
    if nbtasks = 0
     then put(", no tasking");
     else put(", with "); put(nbtasks,1); put(" tasks");
    end if;
    if precision = 1 then
      put_line(", in double double precision.");
    elsif precision = 2 then
      put_line(", in quad double precision.");
    else
      put_line(", in double precision.");
    end if;
  end Write_Greeting;

  procedure Write_Conclusion
              ( start_moment : in Ada.Calendar.Time;
                nbtasks : in natural32 ) is

    ended_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

  begin
    new_line;
    put("PHC ran from ");
    Time_Stamps.Write_Time_Stamp(standard_output,start_moment);
    put(" till ");
    Time_Stamps.Write_Time_Stamp(standard_output,ended_moment);
    put_line(".");
    Time_Stamps.Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    Write_Number_of_Tasks(standard_output,nbtasks);
    Write_Seed_Number(standard_output);
    put_line(Greeting_Banners.Version);
  end Write_Conclusion;

  procedure Write_Conclusion
              ( file : in file_type; start_moment : in Ada.Calendar.Time;
                nbtasks : in natural32 ) is

    ended_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

  begin
    new_line(file);
    put(file,"PHC ran from ");
    Time_Stamps.Write_Time_Stamp(file,start_moment);
    put(file," till ");
    Time_Stamps.Write_Time_Stamp(file,ended_moment);
    put_line(file,".");
    Time_Stamps.Write_Elapsed_Time(file,start_moment,ended_moment);
    Write_Number_of_Tasks(file,nbtasks);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end Write_Conclusion;

end Greetings_and_Conclusions;
