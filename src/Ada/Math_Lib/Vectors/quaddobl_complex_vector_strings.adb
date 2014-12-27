with String_Splitters;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;
with Standard_Complex_Vector_Strings;

package body QuadDobl_Complex_Vector_Strings is

  function Write ( v : Vector ) return string is

    res : String_Splitters.Link_to_String; 

  begin
    for i in v'range loop
      declare
        nb : string(1..160);
      begin
        QuadDobl_Complex_Numbers_io.put(nb,v(i));
        if i = v'last
         then String_Splitters.Append(res,nb);
         else String_Splitters.Append(res,nb & ASCII.LF);
        end if;
      end;
    end loop;
    return res.all;
  end Write;

  function Parse ( s : string ) return Vector is

    cnt : constant integer32
        := Standard_Complex_Vector_Strings.Count_Linefeeds(s) + 1;
    res : Vector(1..cnt);
    start : integer := s'first;
    cff : Complex_Number;
    pos,last : integer;

  begin
    for i in res'range loop
      pos := Standard_Complex_Vector_Strings.Next_Linefeed(s(start..s'last));
      if pos > s'last 
       then QuadDobl_Complex_Numbers_io.get(s(start..s'last),cff,last);
       else QuadDobl_Complex_Numbers_io.get(s(start..pos),cff,last);
      end if;
      res(i) := cff;
      exit when pos > s'last;
      start := pos + 1; -- skip the linefeed
      exit when start > s'last; -- careful with empty lines ...
    end loop;
    return res;
  end Parse;

end QuadDobl_Complex_Vector_Strings; 
