with Ada.Calendar;    use Ada.Calendar;

function Bye_Bye_Message return string is

  date : constant time := Clock;

  function Convert1 ( n : integer ) return character is
  begin
    case n is
      when 1 => return '1';
      when 2 => return '2';
      when 3 => return '3';
      when 4 => return '4';
      when 5 => return '5';
      when 6 => return '6';
      when 7 => return '7';
      when 8 => return '8';
      when 9 => return '9';
      when others => return '0';
    end case;
  end Convert1;

  function Convert2 ( n : integer ) return string is
  begin
    if n < 10 then
      declare
        res : string(1..1);
      begin
        res(1) := Convert1(n);
        return res;
      end;
    else
      declare
        n10 : constant integer := n/10;
        n1 : constant integer := n mod 10;
        last : string(1..1);
        first : constant string := Convert2(n10);
      begin
        last(1) := Convert1(n1);
        return first & last;
      end;
    end if;
  end Convert2;

  function Convert ( d : Time ) return string is

    yr : Year_Number;
    mo : Month_Number;
    da : Day_Number;
    se : Day_Duration;

  begin
    Split(d,yr,mo,da,se);
    return Convert2(da) & "/" & Convert2(mo) & "/" & Convert2(yr);
  end Convert;

begin
  return "PHC worked" & " on " & Convert(date) & ".";
end Bye_Bye_Message;
