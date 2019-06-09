with Characters_and_Numbers;             use Characters_and_Numbers;
with Number_of_Cores;

package body Actions_and_Options is

  function Position ( s : string; c : character ) return integer32 is
  begin
    for i in s'range loop
      if s(i) = c
       then return integer32(i);
      end if;
    end loop;
    return integer32(s'first-1);
  end Position;

  function Number_of_Tasks ( args : Array_of_Strings ) return natural32 is

    cnt : constant integer32 := Number_of_Cores;
    res : natural32 := 0;

  begin
    for i in 1..args'last loop
      declare
        s : constant string := args(i).all;
      begin
        if s'last > 1 then
          if s(2) = 't' then
            if s(3..s'last) = "" then -- no number of cores provided
              if cnt > 0
               then res := natural32(cnt);
              end if;
            else
              res := Convert(s(3..s'last));
            end if;
            exit;
          end if;
        end if;
      end;
    end loop;
    return res;
  end Number_of_Tasks;

  function Verbose_Level ( args : Array_of_Strings ) return integer32 is

    res : natural32 := 0;

  begin
    for i in 1..args'last loop
      declare
        s : constant string := args(i).all;
      begin
        if s'last > 1 then
          if s(2) = 'V' then
            if s(3..s'last) /= ""    -- if verbose level provided
             then res := Convert(s(3..s'last));
            end if;
            exit;
          end if;
        end if;
      end;
    end loop;
    return integer32(res);
  end Verbose_Level;

  function Find_Seed ( args : Array_of_Strings ) return natural32 is

    res : natural32 := 0;

  begin
    for i in 1..args'last loop
      declare
        s : constant string := args(i).all;
      begin
        if s'last > 0 then
          if s(2) = '0'
           then res := Convert(s(3..s'last)); exit;
          end if;
        end if;
      end;
    end loop;
    return res;
  end Find_Seed;

  function Scan_Precision
             ( args : Array_of_Strings; opt : character ) return natural32 is

    res : natural32 := 1;

  begin
    for i in 1..args'last loop
      declare
        s : constant string := args(i).all;
      begin
        if s(2) = opt then
          if s'last > 2 
           then res := Convert(s(3..s'last)); exit;
          end if;
        end if;
      end;
    end loop;
    return res;
  end Scan_Precision;

end Actions_and_Options;
