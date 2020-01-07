with integer_io;                      use integer_io;
with duration_io;

package body Time_Stamps is

-- AUXILIARY FUNCTION :

  function Month_Name ( m : integer ) return string is

  -- DESCRIPTION :
  --   Returns the name of the m-th month.

  begin
    case m is
      when 1 => return "January";
      when 2 => return "February";
      when 3 => return "March";
      when 4 => return "April";
      when 5 => return "May";
      when 6 => return "June";
      when 7 => return "July";
      when 8 => return "August";
      when 9 => return "September";
      when 10 => return "October";
      when 11 => return "November";
      when 12 => return "December";
      when others => return "";
    end case;
  end Month_Name;

-- TARGET ROUTINES :

  procedure Seconds_into_HMS ( seconds : in Duration;
                               hour,min,sec : out integer ) is

    remainder : integer;

  begin
    hour := integer(seconds)/3600;
    remainder := integer(seconds)-hour*3600;
    min := remainder/60;
    remainder := remainder-min*60;
    sec := remainder;
  end Seconds_into_HMS;

  procedure Seconds_into_HMSMS
              ( seconds : in Duration;
                hour,min,sec,millisec : out integer ) is

    remainder : integer;
    wholeseconds : constant Duration := duration(integer(seconds));
    milliseconds : constant Duration := seconds - wholeseconds;

  begin
    hour := integer(seconds)/3600;
    remainder := integer(seconds)-hour*3600;
    min := remainder/60;
    remainder := remainder-min*60;
    sec := remainder;
    millisec := integer(milliseconds*1000.0);
  end Seconds_into_HMSMS;

  function Elapsed_Time ( before,after : Time ) return integer is

    res,befyr,befmo,befto,aftyr,aftmo,aftto : integer;
    befse,aftse : Duration;

  begin
    Split(before,befyr,befmo,befto,befse);
    Split(after,aftyr,aftmo,aftto,aftse);
    res := integer(Seconds(after) - Seconds(before));
    if aftto > befto
     then res := res + 86400*(aftto-befto);   -- 86400 seconds in a day
    end if;
    return res;
  end Elapsed_Time;

  procedure Write_Time_Stamp ( file : in file_type; moment : in Time ) is

    the_seconds : Duration;
    year,month,today,hours,minutes,seconds : integer;

  begin
    Split(moment,year,month,today,the_seconds);
    put(file,today,1); put(file," ");
    put(file,Month_Name(month)); put(file," ");
    put(file,year,4); put(file,", ");
    Seconds_into_HMS(the_seconds,hours,minutes,seconds);
    if hours < 10
     then put(file,"0");
    end if;
    put(file,hours,1); put(file,":");
    if minutes < 10
     then put(file,"0");
    end if;
    put(file,minutes,1); put(file,":");
    if seconds < 10
     then put(file,"0");
    end if;
    put(file,seconds,1);
  end Write_Time_Stamp;

  procedure Write_Elapsed_Time
               ( file : in file_type; before,after : in Time ) is
  
    the_seconds : constant duration := after - before;
   -- difsec : constant integer := Elapsed_Time(before,after);
    hours,minutes,remainder : integer;
    float_seconds : constant float := float(the_seconds);
    wholeseconds : constant duration
                 := duration(float'truncation(float_seconds));
    milliseconds : constant integer
                 := integer((the_seconds - wholeseconds)*1000.0);
    difsec : constant integer := integer(wholeseconds);

  begin
    put(file,"The total elapsed wall clock time is ");
    duration_io.put(file,the_seconds,1,3);
    put(file," seconds");
    if difsec < 60 then
      if milliseconds > 0 then
        new_line(file); put(file," = ");
        if difsec > 0 then
          put(file,difsec,1);
          if difsec = 1
           then put(file," second ");
           else put(file," seconds ");
          end if;
        end if;
        put(file,milliseconds,1);
        if milliseconds = 1
         then put_line(file," millisecond.");
         else put_line(file," milliseconds.");
        end if;
      else
        put_line(file,".");
      end if;
    else
      new_line(file); put(file," =");
      hours := difsec/3600;
      if hours > 0 then
        put(file," "); put(file,hours,1);
        if hours = 1 
         then put(file," hour");
         else put(file," hours");
        end if;
      end if;
      remainder := difsec-hours*3600;
      minutes := remainder/60;
      if minutes > 0 then
        put(file," "); put(file,minutes,1);
        if minutes = 1 
         then put(file," minute");
         else put(file," minutes");
        end if;
      end if;
      remainder := remainder-minutes*60;
      if remainder > 0 then
        put(file," "); put(file,remainder,1);
        if remainder = 1
         then put(file," second");
         else put(file," seconds");
        end if;
      end if;
      if milliseconds > 0 then
        put(file," "); put(file,milliseconds,1);
        if milliseconds = 1
         then put_line(file," millisecond.");
         else put_line(file," milliseconds.");
        end if;
      else
        put_line(file,".");
      end if;
    end if;
  end Write_Elapsed_Time;

end Time_Stamps;
