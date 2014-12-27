with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;
with Standard_Complex_Vector_Strings;

package body Multprec_Complex_Vector_Strings is

  function Write ( v : Vector ) return string is

    res : String_Splitters.Link_to_String; 

  begin
    for i in v'range loop
      declare
        cs : constant natural32
           := Multprec_Complex_Numbers_io.Character_Size(v(i));
        nb : string(1..integer(cs));
      begin
        Multprec_Complex_Numbers_io.put(nb,v(i));
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
    pos,last : integer;

  begin
    for i in res'range loop
      pos := Standard_Complex_Vector_Strings.Next_Linefeed(s(start..s'last));
      if pos > s'last 
       then Multprec_Complex_Numbers_io.get(s(start..s'last),res(i),last);
       else Multprec_Complex_Numbers_io.get(s(start..pos),res(i),last);
      end if;
      exit when pos > s'last;
      start := pos + 1; -- skip the linefeed
      exit when start > s'last; -- careful with empty lines ...
    end loop;
    return res;
  end Parse;

end Multprec_Complex_Vector_Strings; 
