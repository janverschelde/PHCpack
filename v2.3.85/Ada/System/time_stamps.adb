with integer_io;                      use integer_io;

package body Time_Stamps is

-- AUXILIARY FUNCTIONS :

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

  procedure Seconds_into_HMS ( seconds : in Duration;
                               hour,min,sec : out integer ) is

  -- DESCRIPTION :
  --   Splits the seconds into hours, minutes and seconds.

    remainder : integer;

  begin
    hour := integer(seconds)/3600;
    remainder := integer(seconds)-hour*3600;
    min := remainder/60;
    remainder := remainder-min*60;
    sec := remainder;
  end Seconds_into_HMS;

-- TARGET ROUTINES :

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

    difsec : constant integer := Elapsed_Time(before,after);
    hours,minutes,remainder : integer;

  begin
    put(file,"The total elapsed time is ");
    put(file,difsec,1); put(file," seconds");
    if difsec < 60 then
      put_line(file,".");
    else
      put(file," = ");
      hours := difsec/3600;
      if hours > 0
       then put(file,hours,1); put(file," hours ");
      end if;
      remainder := difsec-hours*3600;
      minutes := remainder/60;
      if minutes > 0
       then put(file,minutes,1); put(file," minutes ");
      end if;
      remainder := remainder-minutes*60;
      if remainder > 0
       then put(file,remainder,1); put_line(file," seconds.");
      end if;
    end if;
  end Write_Elapsed_Time;

end Time_Stamps;
