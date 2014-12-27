with String_Splitters;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;

package body Standard_Complex_Vector_Strings is

  function Write ( v : Vector ) return string is

    res : String_Splitters.Link_to_String; 

  begin
    for i in v'range loop
      declare
        nb : string(1..46);
      begin
        Standard_Complex_Numbers_io.put(nb,v(i));
        if i = v'last
         then String_Splitters.Append(res,nb);
         else String_Splitters.Append(res,nb & ASCII.LF);
        end if;
      end;
    end loop;
    return res.all;
  end Write;

  function Count_Linefeeds ( s : string ) return integer32 is

    res : integer32 := 0;

  begin
    for i in s'range loop
      if s(i) = ASCII.LF
       then res := res + 1;
      end if;
     end loop;
     return res;
  end Count_Linefeeds;

  function Next_Linefeed ( s : string ) return integer is 

    res : integer := s'first;

  begin
    while s(res) /= ASCII.LF loop
      res := res + 1;
      exit when (res > s'last);
    end loop;
    return res;
  end Next_Linefeed;

  function Parse ( s : string ) return Vector is

    cnt : constant integer32 := Count_Linefeeds(s) + 1;
    res : Vector(1..cnt);
    start : integer := s'first;
    cff : Complex_Number;
    pos,last : integer;

  begin
    for i in res'range loop
      pos := Next_Linefeed(s(start..s'last));
      if pos > s'last 
       then Standard_Complex_Numbers_io.get(s(start..s'last),cff,last);
       else Standard_Complex_Numbers_io.get(s(start..pos),cff,last);
      end if;
      res(i) := cff;
      exit when pos > s'last;
      start := pos + 1; -- skip the linefeed
      exit when start > s'last; -- careful with empty lines ...
    end loop;
    return res;
  end Parse;

end Standard_Complex_Vector_Strings; 
