with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Standard_Solution_Strings;
with Standard_System_and_Solutions_io;
with Solution_Drops;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;
with PHCpack_Operations;
with Standard_PolySys_Container;
with Standard_Solutions_Container;

package body Standard_Solutions_Interface is

  function Standard_Solutions_Read 
             ( vrblvl : integer32 := 0 ) return integer32 is

    sols : Standard_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Read ...");
    end if;
    Standard_Complex_Solutions_io.Read(sols);
    Standard_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Read.");
      end if;
      return 30;
  end Standard_Solutions_Read;

  function Standard_Solutions_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Read_from_File ...");
    end if;
   -- new_line;
   -- put_line("Opening the file with name " & sv & " ...");
    declare
      file : file_type;
      sols : Standard_Complex_Solutions.Solution_List;
    begin
      Open(file,in_file,sv);
      Standard_Complex_Solutions_io.get(file,sols);
      Standard_Solutions_Container.Clear;
      Standard_Solutions_Container.Initialize(sols);
      exception 
        when NAME_ERROR =>
          put_line("File with name " & sv & " could not be found!");
          return 916;
        when USE_ERROR =>
          put_line("File with name " & sv & " could not be read!");
          return 916;
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Read_from_file.");
      end if;
      return 916;
  end Standard_Solutions_Read_from_File;

  function Standard_System_Solutions_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_System_Solutions_Read_from_File ...");
    end if;
   -- new_line;
   -- put_line("Opening the file with name " & sv & " ...");
    declare
      file : file_type;
      p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
      sols : Standard_Complex_Solutions.Solution_List;
    begin
      Open(file,in_file,sv);
      Standard_System_and_Solutions_io.get(file,p,sols);
      Standard_PolySys_Container.Clear;
      Standard_PolySys_Container.Initialize(p.all);
      Standard_Solutions_Container.Clear;
      Standard_Solutions_Container.Initialize(sols);
      exception 
        when NAME_ERROR =>
          put_line("File with name " & sv & " could not be found!");
          return 544;
        when USE_ERROR =>
          put_line("File with name " & sv & " could not be read!");
          return 544;
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_System_Solutions_Read_from_file.");
      end if;
      return 544;
  end Standard_System_Solutions_Read_from_File;

  function Standard_Solutions_Write
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;
    use Standard_Complex_Solutions_io;

    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Write ...");
    end if;
    if not Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        put(PHCpack_Operations.output_file,
            Length_Of(sols),natural32(Head_Of(sols).n),sols);
      else
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Write.");
      end if;
      return 31;
  end Standard_Solutions_Write;

  function Standard_Solutions_Size 
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Size ...");
    end if;
    Assign(integer32(Standard_Solutions_Container.Length),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Size.");
      end if;
      return 32;
  end Standard_Solutions_Size;

  function Standard_Solutions_Dimension
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Dimension ...");
    end if;
    Assign(integer32(Standard_Solutions_Container.Dimension),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Dimension.");
      end if;
      return 33;
  end Standard_Solutions_Dimension;

  function Standard_Solutions_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
           
    use Standard_Complex_Solutions;

    ls : Link_to_Solution;
    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Get ...");
    end if;
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail
     then return 34;
     else Assign_Solution(ls,b,c); return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Get.");
      end if;
      return 34;
  end Standard_Solutions_Get;

  function Standard_Solutions_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
           
    use Standard_Complex_Solutions;

    ls : Link_to_Solution := Convert_to_Solution(b,c);
    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Read ...");
    end if;
    Standard_Solutions_Container.Replace(k,ls.all,fail);
    Clear(ls);
    if fail
     then return 35;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Read.");
      end if;
      return 35;
  end Standard_Solutions_Update;

  function Standard_Solutions_Add
             ( b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
           
    use Standard_Complex_Solutions;

    ls : constant Link_to_Solution := Convert_to_Solution(b,c);

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Add ...");
    end if;
    Standard_Solutions_Container.Append(ls);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Add.");
      end if;
      return 36;
  end Standard_Solutions_Add;

  function Standard_Solutions_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Link_to_Solution;
    fail : boolean;
    n : natural32;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_String_Size ...");
    end if;
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      Assign(0,b);
      return 200;
    else
      n := Standard_Solution_Strings.Length(ls.all);
      Assign(integer32(n),b);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_String_Size.");
      end if;
      return 200;
  end Standard_Solutions_String_Size;

  function Standard_Solutions_Get_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is


    use Interfaces.C;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Get_String ...");
    end if;
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 201;
    else
      declare
        s  : constant string := Standard_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      end;
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Get_String.");
      end if;
      return 201;
  end Standard_Solutions_Get_String;

  function Standard_Solutions_Drop_by_Index
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant natural32 := natural32(v_a(v_a'first));
    use Standard_Complex_Solutions;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    dropped : constant Solution_List := Solution_Drops.Drop(sols,ind);

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Drop_by_Index ...");
    end if;
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(dropped);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Drop_by_Index.");
      end if;
      return 38;
  end Standard_Solutions_Drop_by_Index;

  function Standard_Solutions_Drop_by_Name
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant integer := integer(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);
    sb : Symbol_Table.Symbol;
    ind : natural32;
    use Standard_Complex_Solutions;
    sols,dropped : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Drop_by_Name ...");
    end if;
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    sols := Standard_Solutions_Container.Retrieve;
    dropped := Solution_Drops.Drop(sols,ind);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(dropped);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Drop_by_Name.");
      end if;
      return 146;
  end Standard_Solutions_Drop_by_Name;

  function Standard_Solutions_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Clear ...");
    end if;
    Standard_Solutions_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Clear.");
      end if;
      return 37;
  end Standard_Solutions_Clear;

end Standard_Solutions_Interface;
