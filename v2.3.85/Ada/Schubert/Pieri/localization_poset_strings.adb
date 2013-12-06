with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Characters_and_Numbers;            use Characters_and_Numbers;

package body Localization_Poset_Strings is

  function Bracket_to_String ( b : Bracket ) return string is

    function convert ( k : integer32 ) return string is

      bk : constant string := nConvert(b(k));

    begin
      if k = b'last 
       then return bk;
       else return bk & " " & convert(k+1);
      end if;
    end convert;

  begin
    return ("[" & convert(b'first) & "]");
  end Bracket_to_String;

  function Node_to_String
            ( top,bottom : Bracket; roco : natural32 ) return string is

    t : constant string := Bracket_to_String(top);
    b : constant string := Bracket_to_String(bottom);
    r : constant string := nConvert(roco);

  begin
    return ("(" & t & "," & b & "," & r & ")");
  end Node_to_String;

  function Nodes_to_String ( nd : Node ) return string is

    s : constant string
      := Node_to_String(nd.top,nd.bottom,natural32(nd.roco));

  begin
    if nd.next_sibling = null
     then return s;
     else return (s & Nodes_to_String(nd.next_sibling.all));
    end if;
  end Nodes_to_String;

  function Poset_to_String ( p : Array_of_Nodes ) return string is

    function level_to_string ( n : in natural32 ) return string is

    -- DESCRIPTION :
    --   Returns the string representation of the current level n
    --   in the poset.

      s : constant string := "n = ";
  
    begin
      if p'last >= 10 and n < 10
       then return (s & " " & nConvert(n) & " : ");
       else return (s & nConvert(n) & " : ");
      end if;
    end level_to_string;

    function convert ( n : in integer32 ) return string is

    -- DESCRIPTION :
    --   Returns the string representation of the poset,
    --   recursively starting at level n.

      s : constant string := level_to_string(natural32(n));

    begin
      if n >= p'last then
        if p(n) = null
         then return s;
         else return (s & Nodes_to_String(p(n).all));
        end if;
      else
        if p(n) = null
         then return (s & ASCII.LF & convert(n+1));
         else return (s & Nodes_to_String(p(n).all) 
                        & ASCII.LF & convert(n+1));
        end if;
      end if;
    end convert;

  begin
    return convert(p'first);
  end Poset_to_String;

end Localization_Poset_Strings;
