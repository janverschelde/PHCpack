with String_Splitters;
with Sets_of_Unknowns_Strings;

package body Partitions_of_Sets_Strings is

  function to_string ( z : Partition ) return string is

    res : String_Splitters.Link_to_String;

  begin
    for i in z'range loop
      String_Splitters.Append(res,Sets_of_Unknowns_Strings.to_string(z(i)));
    end loop;
    declare
      result : constant string := res.all;
    begin
      String_Splitters.Clear(res);
      return result;
    end;
  end to_string;

  function Number_of_Sets ( s : string ) return natural32 is

    res : natural32 := 0;
 
  begin
    for i in s'range loop
      if s(i) = '{'
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Sets;

  function Parse ( s : string; n : natural32 ) return Partition is

    res : Partition(1..Number_of_Sets(s));
    ind : natural32 := 0;
    current : integer := s'first;
    next : integer := current;

  begin
    loop
      while (s(next) /= '}') loop
        next := next + 1;
        exit when (next > s'last);
      end loop;
      ind := ind + 1;
      exit when (ind > n);
      res(ind) := Sets_of_Unknowns_Strings.Parse(s(current..next),n);
      current := next + 1;
      while (s(current) /= '{') loop  -- skip spaces in between sets
        current := current + 1;
        exit when (current > s'last);
      end loop;
      exit when (current > s'last);
      next := current;
    end loop;
    return res;
  end Parse;

end Partitions_of_Sets_Strings;
