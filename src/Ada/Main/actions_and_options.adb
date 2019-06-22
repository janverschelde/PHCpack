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
          if s(1) = '-' and s(2) = 't' then
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
          if s(1) = '-' and s(2) = 'V' then
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
        if s'last > 1 then
          if s(1) = '-' and s(2) = '0'
           then res := Convert(s(3..s'last)); exit;
          end if;
        end if;
      end;
    end loop;
    return res;
  end Find_Seed;

  function Scan_Precision
             ( args : Array_of_Strings; opt : character ) return natural32 is

    res : natural32 := 0; -- better default than 1

  begin
    for i in 1..args'last loop
      declare
        s : constant string := args(i).all;
      begin
        if s'last > 1 then
          if s(1) = '-' and s(2) = opt then
            if s'last > 2 
             then res := Convert(s(3..s'last)); exit;
            end if;
          end if;
        end if;
      end;
    end loop;
    return res;
  end Scan_Precision;

  function Scan_Options ( args : Array_of_Strings ) return string is

    function Scan ( k : in integer; accu : in string ) return string is

    -- DESCRIPTION :
    --   The current argument has index k in args.
    --   Accumulates the options in the string accu.
    --   Returns accu if k is out of range of args.

    begin
      if k < args'first or k > args'last then
        return accu;
      else
        declare
          s : constant string := args(k).all;
        begin
          if s'last <= 1 then
            return Scan(k+1,accu);
          else
            if s(1) = '-'
             then return Scan(k+1,accu & s(2));
             else return Scan(k+1,accu);
            end if;
          end if;
        end;
      end if;
    end Scan;

  begin
    return Scan(1,"");
  end Scan_Options;

  function Sort_Options ( opts : string ) return string is
  begin
    if opts'length <= 1 then
      return opts;
    else
      declare
        op1 : constant character := opts(opts'first);
        pos : constant integer32 := Position(actions,op1);
      begin
        if pos >= integer32(actions'first) then
          return opts;
        else
          return Sort_Options(opts(opts'first+1..opts'last)) & op1;
        end if;
      end;
    end if;
  end Sort_Options;

  function Get_Argument
             ( args : Array_of_Strings; k : integer32 ) return string is

    null_string : constant string := "";
    cnt : integer32 := 0;

  begin
    if args'last > 0 then
      for i in 1..args'last loop
        declare
          s : constant string := args(i).all;
        begin
          if s(1) /= '-' then
            cnt := cnt + 1;
            if k = cnt
             then return s;
            end if;
          end if;
        end;
      end loop;
    end if;
    return null_string;
  end Get_Argument;

end Actions_and_Options;
