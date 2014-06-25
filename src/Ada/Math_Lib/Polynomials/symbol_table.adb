with unchecked_deallocation;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Parse_Polynomial_Exceptions;        use Parse_Polynomial_Exceptions;

package body Symbol_Table is

-- INTERNAL DATA :

  type Symbol_Table ( max : integer32 ) is record
    number : natural32;     -- number of symbols that are not blank
    syms : Array_of_Symbols(1..max);
  end record;
  type Link_to_Symbol_Table is access Symbol_Table;

  st : Link_to_Symbol_Table;

-- CREATORS :

  procedure Init ( max : in natural32 ) is
  begin
    st := new Symbol_Table(integer32(max));
    st.all.number := 0;
  end Init;

  function Standard_Symbols ( n : integer32 ) return Array_of_Symbols is

    res : Array_of_Symbols(1..n);

  begin
    for i in 1..n loop
      declare
        nb : constant string := Characters_and_Numbers.Convert(i);
        sb : Symbol;
      begin
        sb := (sb'range => ' ');
        sb(sb'first) := 'x';
        for j in nb'range loop
          sb(sb'first+1+j-nb'first) := nb(j);
        end loop;
        res(i) := sb;
      end;
    end loop;
    return res;
  end Standard_Symbols;

  procedure Init ( s : in Array_of_Symbols ) is

    max : constant integer32 := s'last - s'first + 1;
    ind : integer32 := 0;

  begin
    st := new Symbol_Table(max);
    st.number := natural32(max);
    for i in s'range loop
      ind := ind+1;
      st.syms(ind) := s(i);
    end loop;
  end Init;

  procedure Enlarge ( max : in natural32 ) is
  begin
    if Empty then
      Init(max);
    else
      declare
        oldst : Array_of_Symbols(1..integer32(st.number));
        maxst : constant natural32 := max + natural32(st.max);
      begin
        for i in oldst'range loop
          oldst(i) := st.syms(i);
        end loop;
        Clear;
        Init(maxst);
        for i in oldst'range loop
          Add(oldst(i));
        end loop;
      end;
    end if;
  end Enlarge;

  procedure Downsize ( n : in natural32 ) is
  begin
    if not Empty then
      if integer32(n) >= st.max then
        Clear;
      else
        declare
          oldst : Array_of_Symbols(1..integer32(st.number));
          maxst : constant natural32 := natural32(st.max) - n;
        begin
          for i in 1..integer32(maxst) loop
            exit when (i > integer32(st.number));
            oldst(i) := st.syms(i);
          end loop;
          Clear;
          Init(maxst);
          for i in 1..maxst loop
            exit when (integer32(i) > oldst'last);
            Add(oldst(integer32(i)));
          end loop;
        end;
      end if;
    end if;
  end Downsize;

  procedure Replace ( i : in natural32; sb : in Symbol ) is

    tab : Symbol_Table renames st.all;

  begin
    if i <= tab.number then
      for j in sb'range loop
        tab.syms(integer32(i))(j) := sb(j);
      end loop;
    end if;
  end Replace;

  function Create ( s : string ) return Symbol is

    res : Symbol;
    k : integer;

  begin
    for i in res'range loop
      res(i) := ' ';
    end loop;
    k := s'first;
    for i in res'range loop
      exit when (k > s'last);
      res(i) := s(k);
      k := k + 1;
    end loop;
    return res;
  end Create;

-- CONSTRUCTORS :

  procedure Add_String ( s : in string ) is

    sb : constant Symbol := Create(s);

  begin
    Add(sb);
  end Add_String;

  procedure Add_String ( s : in string; pos : out natural32 ) is

    sb : constant Symbol := Create(s);

  begin
    Add(sb,pos);
  end Add_String;

  procedure Add ( sb : in Symbol; pos : out natural32 ) is

    tab : Symbol_Table renames st.all;

  begin
    tab.number := tab.number + 1;
    for i in sb'range loop
      tab.syms(integer32(tab.number))(i) := sb(i);
    end loop;
    pos := tab.number;
  exception
    when others => -- put("adding symbol "); Symbol_Table_io.put(sb); new_line;
      raise OVERFLOW_IN_THE_SYMBOL_TABLE;
  end Add;

  procedure Add ( sb : in Symbol ) is

    pos : natural32;

  begin
    Add(sb,pos);
  end Add;

  procedure Remove ( sb : in Symbol ) is

    pos : constant natural32 := Get(sb);

  begin
    Remove(pos);
  end Remove;

  procedure Remove ( i : in natural32 ) is

    tab : Symbol_Table renames st.all;

  begin
    if ((i /= 0) and then (tab.number >= i)) then
      tab.number := tab.number - 1;                -- reduce #symbols
      for j in i..tab.number loop                  -- shift symbol table
        for k in tab.syms(integer32(j))'range loop
          tab.syms(integer32(j))(k) := tab.syms(integer32(j)+1)(k);
        end loop;
      end loop;
    end if;
  end Remove;

-- SELECTORS :

  function Length ( s : Symbol ) return natural32 is

    res : natural32 := 0;

  begin
    for i in s'range loop
      exit when s(i) = ' ';
      res := res + 1;
    end loop;
    return res;
  end Length;

  function Equal ( s1,s2 : Symbol ) return boolean is
  begin
    for i in s1'range loop
      if s1(i) /= s2(i)
       then return false;
      end if;
      exit when ((s1(i) = ' ') and (s2(i) = ' '));
    end loop;
    return true;
  end Equal;

  function Is_Valid ( sb : Symbol ) return boolean is
  begin
    if Convert(sb(1)) < 10 then
      return false;
    else
      for j in 2..sb'last loop
        case sb(j) is
          when '*' | '+' | '-' | '^' | '/' => return false;
          when delimiter | '(' | ')' => return false;
          when others => null;
        end case;
      end loop;
    end if;
    return true;
  end Is_Valid;

  function Check_Symbol ( n : natural32; sb : Symbol ) return natural32 is

    res : natural32;

  begin
    if not Is_Valid(sb)
     then raise INVALID_SYMBOL;
    end if;
    res := Get(sb);         -- search for the number of the symbol
   -- put("in Check_Symbol, res = "); put(res,1); new_line;
    if res = 0 then
      declare
      begin
        Add(sb,res);
      exception
        when OVERFLOW_IN_THE_SYMBOL_TABLE 
          => -- put("looking up symbol "); Symbol_Table_io.put(sb); new_line;
             raise OVERFLOW_OF_UNKNOWNS;
      end;
    end if;
    if res > n then
      -- put("checking symbol "); Symbol_Table_io.put(sb); new_line;
      raise OVERFLOW_OF_UNKNOWNS;
    end if; 
    return res;
  end Check_Symbol;

  function "<" ( s1,s2 : Symbol ) return boolean is
  begin
    for i in s1'range loop
      if s1(i) < s2(i) then
        return true;
      elsif s1(i) > s2(i) then
        return false;
      end if;
    end loop;
    return false;
  end "<";

  function ">" ( s1,s2 : Symbol ) return boolean is
  begin
    for i in s1'range loop
      if s1(i) > s2(i) then
        return true;
      elsif s1(i) < s2(i) then
        return false;
      end if;
    end loop;
    return false;
  end ">";

  function Maximal_Size return natural32 is
  begin
    if st = null
     then return 0;
     else return natural32(st.max);
    end if;
  end Maximal_Size;

  function Number return natural32 is
  begin
    if st = null
     then return 0;
     else return st.all.number;
    end if;
  end Number;
  
  function Empty return boolean is
  begin
    return (st = null);
  end Empty;

  function Get ( sb : Symbol ) return natural32 is

    tab : Symbol_Table renames st.all;

  begin
    for i in 1..tab.number loop
      --if tab.syms(i) = sb
      if Equal(tab.syms(integer32(i)),sb)
       then return i;
      end if;
    end loop;
    return 0;
  end Get;

  function Get ( i : natural32 ) return Symbol is

    tab : Symbol_Table renames st.all;

  begin
    if i > tab.number
     then raise INDEX_OUT_OF_RANGE;
     else return tab.syms(integer32(i));
    end if;
  end Get;

  function Content return Array_of_Symbols is

    res : Array_of_Symbols(st.syms'range);

  begin
    for i in st.syms'range loop
      res(i) := st.syms(i);
    end loop;
    return res;
  end Content;

-- DESTRUCTOR :

  procedure Clear ( ls : in out Link_to_Array_of_Symbols ) is

    procedure free is
      new unchecked_deallocation(Array_of_Symbols,Link_to_Array_of_Symbols);

  begin
    free(ls);
  end Clear;

  procedure Clear is

    procedure free is 
      new unchecked_deallocation(Symbol_Table,Link_to_Symbol_Table);

  begin
    free(st);
  end Clear;

end Symbol_Table;
