with text_io;                           use text_io;
with Interfaces.C;
with File_Scanning;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Standard_Solution_Strings;
with Standard_System_and_Solutions_io;
with Solution_Drops;
with Projective_Transformations;
with Multi_Projective_Transformations;
with Total_Degree_Start_Systems;        use Total_Degree_Start_Systems;
with Lexicographic_Root_Enumeration;    use Lexicographic_Root_Enumeration;
with Drivers_to_Track_Standard_Paths;   use Drivers_to_Track_Standard_Paths;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;
with File_Management;
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
      close(file);
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
      close(file);
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

  function Standard_Solutions_Add_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use Standard_Complex_Solutions;

    va : constant C_Integer_Array 
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(va(va'first));
    nc : constant integer := integer(va(va'first+1));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : String(1..nc);   
    ind : integer := 1;
    sol : Solution(integer32(nv));
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Add_String ...");
    end if;
   -- put_line("Inside the Job 38 ...");
    sv := C_Integer_Array_to_String(natural32(nc),vb);
   -- put_line("The string received in Job 38 : "); put_line(sv);
    Standard_Solution_Strings.Parse(sv,ind,nv,sol,fail);
   -- put_line("The parsed solution : ");
   -- Standard_Complex_Solutions_io.put(sol);
    if fail then
     -- put_line("Failure occurred !");
      return 208;
    else
     -- put_line("Appending the solution to the container...");
      Standard_Solutions_Container.Append(sol);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Add_String.");
      end if;
      return 208;
  end Standard_Solutions_Add_String;

  function Standard_Solutions_Replace_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use Standard_Complex_Solutions;

    va : constant C_Integer_Array 
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    k : constant natural32 := natural32(va(va'first));
    nv : constant natural32 := natural32(va(va'first+1));
    nc : constant integer := integer(va(va'first+2));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : String(1..nc);   
    ind : integer := 1;
    sol : Solution(integer32(nv));
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Replace_String ...");
    end if;
    sv := C_Integer_Array_to_String(natural32(nc),vb);
    Standard_Solution_Strings.Parse(sv,ind,nv,sol,fail);
    if fail then
      return 209;
    else
      Standard_Solutions_Container.Replace(k,sol,fail);
      if fail
       then return 209;
       else return 0;
      end if;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Replace_String.");
      end if;
      return 209;
  end Standard_Solutions_Replace_String;

  function Standard_Solutions_Intro_String_Size
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
      put_line("Standard_Solutions_Intro_String_Size ...");
    end if;
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      Assign(0,b);
      return 202;
    else
      n := Standard_Solution_Strings.Length_Intro(ls.all);
      Assign(integer32(n),b);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Intro_String_Size.");
      end if;
      return 202;
  end Standard_Solutions_Intro_String_Size;

  function Standard_Solutions_Intro_String
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
    s : string(1..integer(n));
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Intro_String ...");
    end if;
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 203;
    else
      s := Standard_Solution_Strings.Write_Intro(ls.all);
      sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Intro_String.");
      end if;
      return 203;
  end Standard_Solutions_Intro_String;

  function Standard_Solutions_Vector_String_Size
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
      put_line("Standard_Solutions_Vector_String_Size ...");
    end if;
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      Assign(0,b);
      return 204;
    else
      n := Standard_Solution_Strings.Length_Vector(ls.all);
      Assign(integer32(n),b);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Vector_String_Size.");
      end if;
      return 204;
  end Standard_Solutions_Vector_String_Size;

  function Standard_Solutions_Vector_String
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
    s : string(1..integer(n));
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Vector_String ...");
    end if;
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 205;
    else
      s := Standard_Solution_Strings.Write_Vector(ls.all);
      sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Vector_String.");
      end if;
      return 205;
  end Standard_Solutions_Vector_String;

  function Standard_Solutions_Diagnostics_String_Size
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : natural32;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Diagnostics_String_Size ...");
    end if;
    n := Standard_Solution_Strings.Length_Diagnostics;
    Assign(integer32(n),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Diagnostics_String_Size.");
      end if;
      return 206;
  end Standard_Solutions_Diagnostics_String_Size;

  function Standard_Solutions_Diagnostics_String
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
    s : string(1..integer(n));
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Diagnostics_String ...");
    end if;
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 207;
    else
      s := Standard_Solution_Strings.Write_Diagnostics(ls.all);
      sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Diagnostics_String.");
      end if;
      return 207;
  end Standard_Solutions_Diagnostics_String;

  function Standard_Solutions_Retrieve_Next
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;

    ls : Link_to_Solution;
    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v(v'first));
    idx : natural32;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Retrieve_Next ...");
    end if;
    if k = 0 then
      Standard_Solutions_Container.Retrieve_Next_Initialize;
    else
      Standard_Solutions_Container.Retrieve_Next(ls,idx);
      Assign(integer32(idx),a);
      if idx = 0
       then return 276;
       else Assign_Solution(ls,b,c);
      end if;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Retrieve_Next.");
      end if;
      return 276;
  end Standard_Solutions_Retrieve_Next;

  function Standard_Solutions_Move_Pointer
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    ind : natural32;

  begin
    if vrblvl > 0 then
      put("-> in Standard_solutions_interface.");
      put_line("Standard_Solutions_Move_Pointer ...");
    end if;
    Standard_Solutions_Container.Move_Current(ind);
    Assign(integer32(ind),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Move_Pointer.");
      end if;
      return 454;
  end Standard_Solutions_Move_Pointer;

  function Standard_Solutions_Current_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;

    ind,len : natural32;
    ls : Link_to_Solution;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Current_Size ...");
    end if;
    Standard_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      len := Standard_Solution_Strings.Length(ls.all);
      Assign(integer32(len),b);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Current_Size.");
      end if;
      return 525;
  end Standard_Solutions_Current_Size;

  function Standard_Solutions_Current_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    len : constant natural32 := natural32(v_a(v_a'first));
    ind : natural32;
    ls : Link_to_Solution;
    sv : Standard_Integer_Vectors.Vector(1..integer32(len));

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Current_String ...");
    end if;
    Standard_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      declare
        s : constant string := Standard_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(len,s));
      end;
      Assign(sv,b);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Current_String.");
      end if;
      return 533;
  end Standard_Solutions_Current_String;

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

  function Standard_Solutions_Make_Homogeneous
             ( vrblvl : integer32 := 0 ) return integer32 is

    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Make_Homogeneous ...");
    end if;
    if not Standard_Complex_Solutions.Is_Null(sols)
     then Projective_Transformations.Projective_Transformation(sols);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Make_Homogeneous.");
      end if;
      return 894;
  end Standard_Solutions_Make_Homogeneous;

  function Standard_Solutions_Multi_Homogeneous
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    m : constant natural32 := natural32(v_a(v_a'first));
    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Multi_Homogeneous ...");
    end if;
    if not Standard_Complex_Solutions.Is_Null(sols)
     then Multi_Projective_Transformations.Add_Ones(sols,m);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Multi_Homogeneous.");
      end if;
      return 910;
  end Standard_Solutions_Multi_Homogeneous;

  function Standard_Solutions_1Hom2Affine
             ( vrblvl : integer32 := 0 ) return integer32 is

    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_1Hom2Affine ...");
    end if;
    if not Standard_Complex_Solutions.Is_Null(sols)
     then Projective_Transformations.Affine_Transformation(sols);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_1Hom2Affine.");
      end if;
      return 898;
  end Standard_Solutions_1Hom2Affine;

  function Standard_Solutions_mHom2Affine
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    mhom : constant natural32 := natural32(v_a(v_a'first+1));
    idz : Standard_Natural_Vectors.Vector(1..integer32(nvr));
    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_mHom2Affine ...");
    end if;
    Assign(nvr,b,idz);
    if not Standard_Complex_Solutions.Is_Null(sols)
     then Multi_Projective_Transformations.Make_Affine(sols,mhom,idz);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_mHom2Affine.");
      end if;
      return 913;
  end Standard_Solutions_mHom2Affine;

  function Standard_Solutions_Tzero
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;

    sols : Solution_List := Standard_Solutions_Container.Retrieve;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

  begin
    if vrblvl > 0 then
      put("-> in standard_solution_interface.");
      put_line("Standard_Solutions_Tzero ...");
    end if;
    if not Is_Null(sols)
     then Set_Continuation_Parameter(sols,zero);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Tzero.");
      end if;
      return 875;
  end Standard_Solutions_Tzero;

  function Standard_Solutions_Read_Next
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;
    use Standard_Complex_Solutions_io;

    dim : natural32;
    ls : Link_to_Solution;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Read_Next ...");
    end if;
    Assign(a,integer32(dim));
   -- put("Dimension : "); put(dim,1); put_line(", calling Read_Next ...");
    Read_Next(File_Management.Link_to_Input.all,dim,ls);
   -- put_line("The solution read : "); put(ls.all); new_line;
    Assign_Solution(ls,b,c);
    Clear(ls);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Read_Next.");
      end if;
      return 135;
  end Standard_Solutions_Read_Next;

  function Standard_Solutions_Write_Next
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;
    use Standard_Complex_Solutions_io;

    cnt : natural32;
    ls : Link_to_Solution := Convert_to_Solution(b,c);

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Write_Next ...");
    end if;
    Assign(a,integer32(cnt));
    Write_Next(File_Management.Link_to_Output.all,cnt,ls);
   -- put_line("Written solution : "); put(ls.all); new_line;
    Assign(integer32(cnt),a);
    Clear(ls);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Write_Next.");
      end if;
      return 136;
  end Standard_Solutions_Write_Next;

  function Standard_Solutions_Next_to_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;
    use Standard_Complex_Solutions_io;

    cnt : natural32;
    ls : Link_to_Solution := Convert_to_Solution(b,c);

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Next_to_File ...");
    end if;
    Assign(a,integer32(cnt));
    if PHCpack_Operations.Is_File_Defined
     then Write_Next(PHCpack_Operations.output_file,cnt,ls);
     else Write_Next(standard_output,cnt,ls);
    end if;
   -- put_line("Written solution : "); put(ls.all); new_line;
    Assign(integer32(cnt),a);
    Clear(ls);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Next_to_File.");
      end if;
      return 141;
  end Standard_Solutions_Next_to_File;

  function Standard_Solutions_Total_Degree
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant natural32 := natural32(v_a(v_a'first));
    i : constant natural32 := natural32(v_a(v_a'first+1));
    lq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    d,s : Standard_Natural_Vectors.Vector(1..integer32(n));
    cff,sol : Standard_Complex_Vectors.Vector(1..integer32(n));
    ls : Link_to_Solution;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Total_Degree ...");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lq);
    d := Degrees(lq.all);
    cff := Coefficients(lq.all);
    s := Root_Map(n,i,d);
    sol := Root(d,s,cff);
    ls := new Solution'(Create(sol));
    Assign_Solution(ls,b,c);
    Clear(ls);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Total_Degree.");
      end if;
      return 142;
  end Standard_Solutions_Total_Degree;

  function Standard_Solutions_Next_Product
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant natural32 := natural32(v_a(v_a'first));
    i : constant natural32 := natural32(v_a(v_a'first+1));
    lq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    d : Standard_Natural_Vectors.Vector(1..integer32(n));
   -- cp : Standard_Natural_Vectors.Vector(1..n-1);
    ls : Link_to_Solution;
    cnt,len : natural32;
    tol : constant double_float := 1.0E-10;
    fail : boolean;
    new_a : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Next_Product ...");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lq);
    d := Degrees(lq.all);
    len := Product(d);
   -- cp := Consecutive_Products(d);
    cnt := i;
   -- Next_Lex_Linear_Product_Start_Solution
   --   (false,n,d,cp,cnt,len,5,tol,ls,fail);
    Get_Next_Linear_Product_Start_Solution(false,n,cnt,len,5,tol,ls,fail);
    new_a(1) := n;   
    new_a(2) := cnt;
    Assign(new_a,a);
    if fail then
      return 143;
    else
      Assign_Solution(ls,b,c);
      Clear(ls);
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Next_Product.");
      end if;
      return 143;
  end Standard_Solutions_Next_Product;

  function Standard_Solutions_Lex_Product
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant natural32 := natural32(v_a(v_a'first));
    i : constant natural32 := natural32(v_a(v_a'first+1));
    lq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    d : Standard_Natural_Vectors.Vector(1..integer32(n));
    cp : Standard_Natural_Vectors.Vector(1..integer32(n)-1);
    ls : Link_to_Solution;
    tol : constant double_float := 1.0E-10;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Lex_Product ...");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lq);
    d := Degrees(lq.all);
    cp := Consecutive_Products(d);
    Get_Lex_Linear_Product_Start_Solution(false,n,d,cp,i,5,tol,ls,fail);
    if fail then
      return 144;
    else
      Assign_Solution(ls,b,c);
      Clear(ls);
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Lex_Product.");
      end if;
      return 144;
  end Standard_Solutions_Lex_Product;

  function Standard_Solutions_Next_Witness
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use Standard_Complex_Solutions;
    use Standard_Complex_Solutions_io;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;

  begin
   -- put("reading next witness point from set "); put(k,1); new_line;
   -- put("  solution vector has length "); put(n,1); new_line;
    Read_Next(File_Management.Link_to_Input(k).all,n,ls,
              Standard_Solutions_Container.Retrieve_Symbol_Table(0).all);
   -- was the following:
   --           Standard_Solutions_Container.Retrieve_Symbol_Table(k).all);
   -- put_line("The solution vector read : "); put(ls.all); new_line;
    Assign_Solution(ls,b,c);
    Clear(ls);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Next_Witness.");
      end if;
      return 145;
  end Standard_Solutions_Next_Witness;

  function Standard_Solutions_Scan_Banner
             ( vrblvl : integer32 := 0 ) return integer32 is

    use File_Scanning;

    found : boolean;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Scan_Banner.");
    end if;
    Scan_and_Skip(File_Management.Link_to_Input.all,"SOLUTIONS",found);
    if found
     then return 0;
     else return 132;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Scan_Banner.");
      end if;
      return 132;
  end Standard_Solutions_Scan_Banner;

  function Standard_Solutions_Read_Dimensions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions_io;

    len,dim : natural32;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Read_Dimensions.");
    end if;
    Read_First(File_Management.Link_to_Input.all,len,dim);
    Assign(integer32(len),a);
    Assign(integer32(dim),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Read_Dimensions.");
      end if;
      Assign(0,a); Assign(0,b);
      return 133;
  end Standard_Solutions_Read_Dimensions;

  function Standard_Solutions_Write_Dimensions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions_io;

    len,dim : natural32;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Write_Dimensions.");
    end if;
    Assign(a,integer32(len));
    Assign(b,integer32(dim));
    Write_First(File_Management.Link_to_Output.all,len,dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Write_Dimensions.");
      end if;
      return 134;
  end Standard_Solutions_Write_Dimensions;

  function Standard_Solutions_Close_Input_File
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v_a(v_a'first));

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Close_Input_File.");
    end if;
    if k = 0
     then File_Management.Close_Input_File;
     else File_Management.Close_Input_File(k);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Close_Input_File.");
      end if;
      return 137;
  end Standard_Solutions_Close_Input_File;

  function Standard_Solutions_Close_Output_File
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Close_Output_File.");
    end if;
    File_Management.Close_Output_File;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Close_Output_File.");
      end if;
      return 138;
  end Standard_Solutions_Close_Output_File;

  function Standard_Solutions_Banner_to_Output
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Banner_to_Output.");
    end if;
    if PHCpack_Operations.Is_File_Defined then
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE SOLUTIONS");
    else
      new_line(standard_output);
      put_line(standard_output,"THE SOLUTIONS");
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Banner_to_Output.");
      end if;
      return 139;
  end Standard_Solutions_Banner_to_Output;

  function Standard_Solutions_Dimensions_to_Output
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions_io;

    len,dim : natural32;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Dimensions_to_Output.");
    end if;
    Assign(a,integer32(len));
    Assign(b,integer32(dim));
    if PHCpack_Operations.Is_File_Defined
     then Write_First(PHCpack_Operations.output_file,len,dim);
     else Write_First(standard_output,len,dim);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solutions_interface.");
        put_line("Standard_Solutions_Dimensions_to_Output.");
      end if;
      return 140;
  end Standard_Solutions_Dimensions_to_Output;

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
