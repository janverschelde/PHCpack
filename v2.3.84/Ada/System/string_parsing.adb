package body String_Parsing is

  function Scan ( s : string; banner : string ) return integer is

    ind : integer := banner'first;

  begin
    for i in s'range loop
      if s(i) = banner(ind) then
        if ind = banner'last
         then return i;
         else ind := ind + 1;
        end if;
      else
        ind := banner'first;
      end if;
    end loop;
    return -1;
  end Scan;

end String_Parsing;
