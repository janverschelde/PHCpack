with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with String_Splitters;
with Standard_Integer_Vectors;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Poly_Systems;
with Parse_Dimensions;
with Witness_Sets_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;

package body Symbol_Table_Interface is

  function Symbol_Table_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in symbol_table_interface.Symbol_Table_Size ...");
    end if;
    Assign(integer32(Symbol_Table.Number),a); 
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in symbol_table_interface.");
        put_line("Symbol_Table_Size.");
      end if;
      return 293;
  end Symbol_Table_Size;

  function Symbol_Table_Write
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in symbol_table_interface.Symbol_Table_Write ...");
    end if;
    for i in 1..Symbol_Table.Number loop
      put(" "); Symbol_Table_io.put(Symbol_Table.Get(i));
    end loop;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in symbol_table_interface.");
        put_line("Symbol_Table_Write.");
      end if;
      return 294;
  end Symbol_Table_Write;

  function Characters_of_Symbol ( i : natural32 ) return string is

  -- DESCRIPTION :
  --   Returns the non-blank characters of the i-th symbol.

    sb : constant Symbol_Table.Symbol := Symbol_Table.Get(i);
    res : string(sb'range);
    cnt : natural32 := 0;

  begin
    for i in sb'range loop
      exit when (sb(i) = ' ');
      cnt := cnt + 1;
      res(i) := sb(i);
    end loop;
    return res(1..integer(cnt));
  end Characters_of_Symbol;

  function String_of_Symbols ( i : natural32 ) return string is

  -- DESCRIPTION :
  --   Returns a string of all symbols in the symbol table,
  --   separated by one space.

  begin
    if i > Symbol_Table.Number then
      return "";
    elsif i = Symbol_Table.Number then
      return Characters_of_Symbol(i);
    else
      return Characters_of_Symbol(i) & " " & String_of_Symbols(i+1);
    end if;  
  end String_of_Symbols;

  function Symbol_Table_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    s : constant string := String_of_Symbols(1);
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(s);

  begin
    if vrblvl > 0
     then put_line("-> in symbol_table_interface.Symbol_Table_String ...");
    end if;
    Assign(integer32(s'last),a);
    Assign(sv,b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in symbol_table_interface.");
        put_line("Symbol_Table_String.");
      end if;
      return 295;
  end Symbol_Table_String;

  function Symbol_Table_Sort_Embedded 
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    dim : natural32;

  begin
    if vrblvl > 0 then
      put("-> in symbol_table_interface.");
      put_line("Symbol_Table_Sort_Embedded ...");
    end if;
    if lp /= null then
      dim := Witness_Sets_io.Count_Embed_Symbols(natural32(lp'last),"zz");
      if dim > 0 then
        Witness_Sets_io.Swap_Symbols_to_End(natural32(lp'last),dim,"zz",lp.all);
        if dim > 1 then
          Witness_Sets_io.Sort_Embed_Symbols
            (natural32(lp'last),natural32(lp'last)-dim,dim,lp.all);
        end if;
      end if;
      Assign(integer32(dim),a);
    else
      Assign(-1,a);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in symbol_table_interface.");
        put_line("Symbol_Table_Sort_Embedded.");
      end if;
      return 292;
  end Symbol_Table_Sort_Embedded;

  function Symbol_Table_Add
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is 

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant integer := integer(v_a(v_a'first));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    smb : constant String(1..nc)
        := C_Integer_Array_to_String(natural32(nc),v_b);

  begin
    if vrblvl > 0
     then put_line("-> in symbol_table_interface.Symbol_Table_Add ...");
    end if;
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add_String(smb);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in symbol_table_interface.");
        put_line("Symbol_Table_Add.");
      end if;
      return 897;
  end Symbol_Table_Add;

  function Symbol_Table_Remove_by_Index
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));

  begin
    if vrblvl > 0 then
      put("-> in symbol_table_interface.");
      put_line("Symbol_Table_Remove_by_Index ...");
    end if;
    Symbol_Table.Remove(n);
    Symbol_Table.Downsize(1);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in symbol_table_interface.");
        put_line("Symbol_Table_Remove_by_Index.");
      end if;
      return 291;
  end Symbol_Table_Remove_by_Index;

  function Symbol_Table_Remove_by_Name
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant natural32 := natural32(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..integer(nc)) := C_Integer_Array_to_String(nc,vb);   
    sb : Symbol_Table.Symbol;
    ind : natural32;

  begin
    if vrblvl > 0 then
      put("-> in symbol_table_interface.");
      put_line("Symbol_Table_Remove_by_Name ...");
    end if;
    for i in 1..integer(nc) loop
      sb(i) := sv(i);
    end loop;
    for i in integer(nc)+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    Symbol_Table.Remove(ind);
    Symbol_Table.Downsize(1);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in symbol_table_interface.");
        put_line("Symbol_Table_Remove_by_Name.");
      end if;
      return 296;
  end Symbol_Table_Remove_by_Name;

  function Symbol_Table_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in symbol_table_interface.Symbol_Table_Clear ...");
    end if;
    Symbol_Table.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in symbol_table_interface.");
        put_line("Symbol_Table_Clear.");
      end if;
      return 29;
  end Symbol_Table_Clear;

  function Symbol_Table_Scan
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbc : constant natural32 := natural32(v_a(v_a'first));
    nbc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nbc-1);
    v_b : constant C_Integer_Array(0..nbc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc));
    s : constant String(1..integer(nbc)) := C_Integer_Array_to_String(nbc,v_b);
    cnt : constant natural := String_Splitters.Count_Delimiters(s,';');
    ls : String_Splitters.Array_of_Strings(1..integer(cnt))
       := String_Splitters.Split(cnt,s,';');
    dim : natural32;
    maxvar : constant natural32 := 2*natural32(cnt); -- better than 1024

  begin
    if vrblvl > 0
     then put_line("-> in symbol_table_interface.Symbol_Table_Scan ...");
    end if;
    dim := Parse_Dimensions.Dim(maxvar,ls);
    Assign(integer32(dim),a);
    String_Splitters.Clear(ls);
    Parse_Dimensions.Clear;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in symbol_table_interface.");
        put_line("Symbol_Table_Scan.");
      end if;
      return 439;
  end Symbol_Table_Scan;

end Symbol_Table_Interface;
