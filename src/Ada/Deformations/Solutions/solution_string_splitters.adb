with String_Parsing;
with Standard_Solution_Strings;

package body Solution_String_Splitters is

  function Trim_End_to_Newline ( s : string ) return string is

    ind : integer := s'last;

  begin
    loop
      exit when (s(ind) = ASCII.CR);
      exit when (s(ind) = ASCII.LF);
      ind := ind - 1;
    end loop;
    return s(s'first..ind);
  end Trim_End_to_Newline;

  procedure Split_Coordinates
              ( s : in string; m : out integer32;
                t : out Standard_Complex_Numbers.Complex_Number;
                c : out Link_to_String; fail : out boolean ) is

    ind : integer := s'first;
    start : integer;

  begin
    Standard_Solution_Strings.Parse_Intro(s,ind,t,m,fail);
    if not fail then
      ind := String_Parsing.Scan(s(ind+1..s'last),":");
      if ind > 0 then
        ind := ind + 2;  -- scan marks the position of the ":"
                -- + 2 to skip newline after 'the solution for t :'
        loop
          start := ind+1;
          ind := String_Parsing.Scan(s(start..s'last),":");
          exit when (ind < 0);
          Append(c,Trim_End_to_Newline(s(start..ind)));
        end loop;
      end if;
    end if;
  end Split_Coordinates;

  function Coordinates ( s : string ) return Link_to_String is

    res : Link_to_String;
    ind : integer := s'first;
    t : Standard_Complex_Numbers.Complex_Number;
    m : integer32;
    fail : boolean;
    start : integer;

  begin
    Standard_Solution_Strings.Parse_Intro(s,ind,t,m,fail);
    if not fail then
      ind := String_Parsing.Scan(s(ind+1..s'last),":");
      if ind > 0 then
        ind := ind + 2;  -- scan marks the position of the ":"
                -- + 2 to skip newline after 'the solution for t :'
        loop
          start := ind+1;
          ind := String_Parsing.Scan(s(start..s'last),":");
          exit when (ind < 0);
          Append(res,Trim_End_to_Newline(s(start..ind)));
        end loop;
      end if;
    end if;
    return res;
  end Coordinates;

end Solution_String_Splitters;
