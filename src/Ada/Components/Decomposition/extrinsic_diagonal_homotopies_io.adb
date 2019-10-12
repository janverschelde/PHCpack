with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Symbol_Table_io;
with Characters_and_Numbers;
with Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;       use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;       use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;       use QuadDobl_Complex_Solutions_io;

package body Extrinsic_Diagonal_Homotopies_io is

  procedure Write ( sa : in Array_of_Symbols ) is
  begin
    Write(Standard_Output,sa);
  end Write;

  procedure Write ( file : in file_type; sa : in Array_of_Symbols ) is
  begin
    for i in sa'range loop
      put(file," "); Symbol_Table_io.put(file,sa(i));
    end loop;
    new_line(file);
  end Write;

  function Add_Suffix ( sb : Symbol; c : character ) return Symbol is

    res : Symbol := sb;
    ind : integer := sb'first;

  begin
    while sb(ind) /= ' ' loop
      ind := ind + 1;
    end loop;
    res(ind) := c;
    return res;
  end Add_Suffix;

  function Add_Suffix ( sa : Array_of_Symbols; c : character )
                      return Array_of_Symbols is

    res : Array_of_Symbols(sa'range);
 
  begin
    for i in sa'range loop
      res(i) := Add_Suffix(sa(i),c);
    end loop;
    return res;
  end Add_Suffix;

  function Suffix ( sb : Symbol ) return character is

    ind : integer := sb'first;

  begin
    while (sb(ind) /= ' ') and (ind < sb'last) loop
      ind := ind + 1;
    end loop;
    if ind = sb'last then
      return sb(ind);
    else
      return sb(ind-1);
    end if;
  end Suffix;

  function Remove_Suffix ( sb : Symbol ) return Symbol is

    res : Symbol := sb;
    ind : integer := sb'first;

  begin
    while (sb(ind) /= ' ') and (ind < sb'last) loop
      ind := ind + 1;
    end loop;
    res(ind-1) := ' ';
    return res;
  end Remove_Suffix;

  function Get_Symbols return Array_of_Symbols is

    res : Array_of_Symbols(1..integer32(Symbol_Table.Number));

  begin
    for i in 1..Symbol_Table.Number loop
      res(integer32(i)) := Symbol_Table.Get(i);
    end loop;
    return res;
  end Get_Symbols;

  function Get_Link_to_Symbols return Link_to_Array_of_Symbols is

    res : Link_to_Array_of_Symbols;
    res_rep : constant Array_of_Symbols := Get_Symbols;

  begin
    res := new Array_of_Symbols'(res_rep);
    return res;
  end Get_Link_to_Symbols;

  function Suffixed_Symbols ( c : character ) return Array_of_Symbols is

    res : Array_of_Symbols(1..integer32(Symbol_Table.Number));

  begin
    for i in 1..Symbol_Table.Number loop 
      declare
        sb : constant Symbol := Symbol_Table.get(i);
      begin
        res(integer32(i)) := Add_Suffix(sb,c);
      end;
    end loop;
    return res;
  end Suffixed_Symbols;

  function Is_Embed_Symbol ( sb : Symbol ) return boolean is
  begin
    if sb(1) /= 'z' then 
      return false;
    elsif sb(2) /= 'z' then
      return false;
    else
      return true;
    end if;
  end Is_Embed_Symbol;

  function Retrieve_Suffixed_Symbols
              ( n : in integer32; c : character ) return Array_of_Symbols is

    res : Array_of_Symbols(1..n);
    ind : integer32 := 0;

  begin
    for i in 1..Symbol_Table.Number loop
      declare
        sb : constant Symbol := Symbol_Table.get(i);
      begin
        if not Is_Embed_Symbol(sb) then
          if Suffix(sb) = c then
            ind := ind + 1;
            res(ind) := sb;         
          end if;
        end if;
      end;
      exit when (ind = n);
    end loop;
    return res;
  end Retrieve_Suffixed_Symbols;

  procedure Retrieve_Suffixed_Symbols
              ( n : in integer32; c : character;
                s : out Array_of_Symbols; p : out Permutation ) is

    ind : integer32 := 0;

  begin
    for i in 1..Symbol_Table.Number loop
      declare
        sb : constant Symbol := Symbol_Table.get(i);
      begin
        if not Is_Embed_Symbol(sb) then
          if Suffix(sb) = c then
            ind := ind + 1;
            s(ind) := sb;         
            p(ind) := integer32(i);
          end if;
        end if;
      end;
      exit when (ind = n);
    end loop;
  end Retrieve_Suffixed_Symbols;

  function Remove_Embed_Symbols 
              ( sa : Array_of_Symbols ) return Array_of_Symbols is

    res : Array_of_Symbols(sa'range);
    ind : integer32 := res'first-1;

  begin
    for i in sa'range loop
      if not Is_Embed_Symbol(sa(i))
       then ind := ind + 1; res(ind) := sa(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Remove_Embed_Symbols;

  function Look_for_Position 
             ( sa : Array_of_Symbols; sb : Symbol ) return integer32 is
  begin
    for i in sa'range loop
      if Equal(sa(i),sb)
       then return i;
      end if;
    end loop;
    return 0;
  end Look_for_Position;

  function Match_Symbols
              ( s1,s2 : in Array_of_Symbols ) return Permutation is

    prm : Permutation(s2'range);

  begin
    for i in s2'range loop
      prm(i) := Look_for_Position(s1,s2(i));
    end loop;
    return prm;
  end Match_Symbols;

  function Equals_mod_Suffix ( sb1,sb2 : Symbol ) return boolean is

    ind : natural := sb1'first;

  begin
    while sb1(ind) /= ' ' loop
      ind := ind + 1;
    end loop;
    if sb2(ind) /= ' ' then
      return false;
    else
      ind := ind - 1;
      for i in sb1'first..(ind-1) loop
        if sb1(i) /= sb2(i)
         then return false;
        end if;
      end loop;
    end if;
    return true;
  end Equals_mod_Suffix;

  function Search_Position 
              ( sa : Array_of_Symbols; sb : Symbol ) return integer32 is
  begin
    for i in sa'range loop
      if Equals_mod_Suffix(sa(i),sb)
       then return i;
      end if;
    end loop;
    return 0;
  end Search_Position;

  function Matching_Permutation
              ( s1,s2 : Array_of_Symbols ) return Permutation is

    prm : Permutation(s2'range);

  begin
    for i in s2'range loop
     -- put("Matching the "); put(i,1); put("-th symbol ");
     -- Symbol_Table_io.put(s2(i)); new_line;
      prm(i) := Search_Position(s1,s2(i));
    end loop;
    return prm;
  end Matching_Permutation;

  function Combine_Permutations
              ( n,d : integer32; p1,p2,prm : Permutation )
              return Permutation is

    res : Permutation(1..2*n+d);

  begin
    for i in 1..n loop
      res(i) := p1(i);
      res(prm(i)+n) := p2(i);
    end loop;
    for i in 2*n+1..2*n+d loop
      res(i) := i;
    end loop;
    return res;
  end Combine_Permutations;

  procedure Permute ( prm : in Permutation; s : in out Array_of_Symbols ) is

    res : Array_of_Symbols(s'range) := s;

  begin
    for i in prm'range loop
      res(prm(i)) := s(i);
    end loop;
    s := res;
  end Permute;

  procedure Write_Symbol_Table is
  begin
    put_line("The symbol table :");
    for i in 1..Symbol_Table.Number loop
      put(" ");
      Symbol_Table_io.put(Symbol_Table.get(i));
    end loop;
    new_line;
  end Write_Symbol_Table;

  procedure Assign_Symbol_Table ( sa : in Array_of_Symbols ) is
  begin
    Symbol_Table.Clear;
    Symbol_Table.Init(sa'length);
    for i in sa'range loop
      Symbol_Table.Add(sa(i));
    end loop;
  end Assign_Symbol_Table;

  procedure Assign_Symbol_Table ( s1,s2 : in Array_of_Symbols ) is
  begin
    Symbol_Table.Clear;
    Symbol_Table.Init(s1'length+s2'length);
    for i in s1'range loop
      Symbol_Table.Add(s1(i));
    end loop;
    for i in s2'range loop
      Symbol_Table.Add(s2(i));
    end loop;
  end Assign_Symbol_Table;

  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    strd : constant string
         := Characters_and_Numbers.convert(integer32(d));
    filename : constant string := name & "_sw" & strd;
    witfile : file_type;

    use Standard_Complex_Solutions;

  begin
    put(file,"writing to file "); put(file,filename); new_line(file);
    create(witfile,out_file,filename);
    put(witfile,p'last,1); new_line(witfile);
    Standard_Complex_Poly_Systems_io.put(witfile,p);
    new_line(witfile);
    put_line(witfile,"THE SOLUTIONS :");
    put(witfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    close(witfile);
  end Write_Witness_Set;

  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    strd : constant string
         := Characters_and_Numbers.convert(integer32(d));
    filename : constant string := name & "_sw" & strd;
    witfile : file_type;

    use DoblDobl_Complex_Solutions;

  begin
    put(file,"writing to file "); put(file,filename); new_line(file);
    create(witfile,out_file,filename);
    put(witfile,p'last,1); new_line(witfile);
    DoblDobl_Complex_Poly_Systems_io.put(witfile,p);
    new_line(witfile);
    put_line(witfile,"THE SOLUTIONS :");
    put(witfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    close(witfile);
  end Write_Witness_Set;

  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    strd : constant string
         := Characters_and_Numbers.convert(integer32(d));
    filename : constant string := name & "_sw" & strd;
    witfile : file_type;

    use QuadDobl_Complex_Solutions;

  begin
    put(file,"writing to file "); put(file,filename); new_line(file);
    create(witfile,out_file,filename);
    put(witfile,p'last,1); new_line(witfile);
    QuadDobl_Complex_Poly_Systems_io.put(witfile,p);
    new_line(witfile);
    put_line(witfile,"THE SOLUTIONS :");
    put(witfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    close(witfile);
  end Write_Witness_Set;

  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    strd : constant string
         := Characters_and_Numbers.convert(integer32(d));
    filename : constant string := name & "_sw" & strd;
    witfile : file_type;

    use Standard_Complex_Solutions;

  begin
    put(file,"writing to file "); put(file,filename); new_line(file);
    create(witfile,out_file,filename);
    put(witfile,p'last,1); new_line(witfile);
    Standard_Complex_Laur_Systems_io.put(witfile,p);
    new_line(witfile);
    put_line(witfile,"THE SOLUTIONS :");
    put(witfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    close(witfile);
  end Write_Witness_Set;

  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    strd : constant string
         := Characters_and_Numbers.convert(integer32(d));
    filename : constant string := name & "_sw" & strd;
    witfile : file_type;

    use DoblDobl_Complex_Solutions;

  begin
    put(file,"writing to file "); put(file,filename); new_line(file);
    create(witfile,out_file,filename);
    put(witfile,p'last,1); new_line(witfile);
    DoblDobl_Complex_Laur_Systems_io.put(witfile,p);
    new_line(witfile);
    put_line(witfile,"THE SOLUTIONS :");
    put(witfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    close(witfile);
  end Write_Witness_Set;

  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    strd : constant string
         := Characters_and_Numbers.convert(integer32(d));
    filename : constant string := name & "_sw" & strd;
    witfile : file_type;

    use QuadDobl_Complex_Solutions;

  begin
    put(file,"writing to file "); put(file,filename); new_line(file);
    create(witfile,out_file,filename);
    put(witfile,p'last,1); new_line(witfile);
    QuadDobl_Complex_Laur_Systems_io.put(witfile,p);
    new_line(witfile);
    put_line(witfile,"THE SOLUTIONS :");
    put(witfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    close(witfile);
  end Write_Witness_Set;

end Extrinsic_Diagonal_Homotopies_io;
