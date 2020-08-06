with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with DoblDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Symbol_Table;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;
with DoblDobl_Solution_Strings;
with DoblDobl_System_and_Solutions_io;
with Solution_Drops;
with Projective_Transformations;
with Multi_Projective_Transformations;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;
with PHCpack_Operations;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;

package body DoblDobl_Solutions_Interface is

  function DoblDobl_Solutions_Read 
             ( vrblvl : integer32 := 0 ) return integer32 is

    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_Solutions_interface.");
      put_line("DoblDobl_Solutions_Read ...");
    end if;
    DoblDobl_Complex_Solutions_io.Read(sols);
    DoblDobl_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_Solutions_interface.");
        put_line("DoblDobl_Solutions_Read.");
      end if;
      return 340;
  end DoblDobl_Solutions_Read;

  function DoblDobl_Solutions_Read_from_File
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
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_Read_from_File ...");
    end if;
   -- new_line;
   -- put_line("Opening the file with name " & sv & " ...");
    declare
      file : file_type;
      sols : DoblDobl_Complex_Solutions.Solution_List;
    begin
      Open(file,in_file,sv);
      DoblDobl_Complex_Solutions_io.get(file,sols);
      DoblDobl_Solutions_Container.Clear;
      DoblDobl_Solutions_Container.Initialize(sols);
      close(file);
      exception 
        when NAME_ERROR =>
          put_line("File with name " & sv & " could not be found!");
          return 917;
        when USE_ERROR =>
          put_line("File with name " & sv & " could not be read!");
          return 917;
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Read_from_file.");
      end if;
      return 917;
  end DoblDobl_Solutions_Read_from_File;

  function DoblDobl_System_Solutions_Read_from_File
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
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_System_Solutions_Read_from_File ...");
    end if;
   -- new_line;
   -- put_line("Opening the file with name " & sv & " ...");
    declare
      file : file_type;
      p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
      sols : DoblDobl_Complex_Solutions.Solution_List;
    begin
      Open(file,in_file,sv);
      DoblDobl_System_and_Solutions_io.get(file,p,sols);
      DoblDobl_PolySys_Container.Clear;
      DoblDobl_PolySys_Container.Initialize(p.all);
      DoblDobl_Solutions_Container.Clear;
      DoblDobl_Solutions_Container.Initialize(sols);
      close(file);
      exception 
        when NAME_ERROR =>
          put_line("File with name " & sv & " could not be found!");
          return 545;
        when USE_ERROR =>
          put_line("File with name " & sv & " could not be read!");
          return 545;
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_System_Solutions_Read_from_file.");
      end if;
      return 545;
  end DoblDobl_System_Solutions_Read_from_File;

  function DoblDobl_Solutions_Write
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Solutions_io;

    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_Solutions_interface.");
      put_line("DoblDobl_Solutions_Write ...");
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
        put("Exception raised in DoblDobl_Solutions_interface.");
        put_line("DoblDobl_Solutions_Write.");
      end if;
      return 341;
  end DoblDobl_Solutions_Write;

  function DoblDobl_Solutions_Size 
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_Solutions_interface.");
      put_line("DoblDobl_Solutions_Size ...");
    end if;
    Assign(integer32(DoblDobl_Solutions_Container.Length),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_Solutions_interface.");
        put_line("DoblDobl_Solutions_Size.");
      end if;
      return 342;
  end DoblDobl_Solutions_Size;

  function DoblDobl_Solutions_Dimension
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_Solutions_interface.");
      put_line("DoblDobl_Solutions_Dimension ...");
    end if;
    Assign(integer32(DoblDobl_Solutions_Container.Dimension),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_Solutions_interface.");
        put_line("DoblDobl_Solutions_Dimension.");
      end if;
      return 343;
  end DoblDobl_Solutions_Dimension;

  function DoblDobl_Solutions_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
           
    use DoblDobl_Complex_Solutions;

    ls : Link_to_Solution;
    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_Solutions_interface.");
      put_line("DoblDobl_Solutions_Get ...");
    end if;
    DoblDobl_Solutions_Container.Retrieve(k,ls,fail);
    if fail
     then return 344;
     else Assign_Solution(ls,b,c); return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_Solutions_interface.");
        put_line("DoblDobl_Solutions_Get.");
      end if;
      return 344;
  end DoblDobl_Solutions_Get;

  function DoblDobl_Solutions_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
           
    use DoblDobl_Complex_Solutions;

    ls : Link_to_Solution := Convert_to_Solution(b,c);
    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_Solutions_interface.");
      put_line("DoblDobl_Solutions_Read ...");
    end if;
    DoblDobl_Solutions_Container.Replace(k,ls.all,fail);
    Clear(ls);
    if fail
     then return 345;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_Solutions_interface.");
        put_line("DoblDobl_Solutions_Read.");
      end if;
      return 345;
  end DoblDobl_Solutions_Update;

  function DoblDobl_Solutions_Add
             ( b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
           
    use DoblDobl_Complex_Solutions;

    ls : constant Link_to_Solution := Convert_to_Solution(b,c);

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_Solutions_interface.");
      put_line("DoblDobl_Solutions_Add ...");
    end if;
    DoblDobl_Solutions_Container.Append(ls);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_Solutions_interface.");
        put_line("DoblDobl_Solutions_Add.");
      end if;
      return 346;
  end DoblDobl_Solutions_Add;

  function DoblDobl_Solutions_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Link_to_Solution;
    fail : boolean;
    n : natural32;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_String_Size ...");
    end if;
    DoblDobl_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      Assign(0,b);
      return 240;
    else
      n := DoblDobl_Solution_Strings.Length(ls.all);
      Assign(integer32(n),b);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_String_Size.");
      end if;
      return 240;
  end DoblDobl_Solutions_String_Size;

  function DoblDobl_Solutions_Get_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;
    sv : Standard_Integer_Vectors.Vector(1..integer32(n));

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_Get_String ...");
    end if;
    DoblDobl_Solutions_Container.Retrieve(k,ls,fail);
    if fail then
      return 241;
    else
      declare
        s : constant string := DoblDobl_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(n,s));
      end;
      Assign(sv,b);
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Get_String.");
      end if;
      return 241;
  end DoblDobl_Solutions_Get_String;

  function DoblDobl_Solutions_Add_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use DoblDobl_Complex_Solutions;

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
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_Add_String ...");
    end if;
    sv := C_Integer_Array_to_String(natural32(nc),vb);
    DoblDobl_Solution_Strings.Parse(sv,ind,nv,sol,fail);
   -- put_line("the parsed solution : ");
   -- DoblDobl_Complex_Solutions_io.put(sol);
    if fail then
      return 378;
    else
      DoblDobl_Solutions_Container.Append(sol);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Add_String.");
      end if;
      return 378;
  end DoblDobl_Solutions_Add_String;

  function DoblDobl_Solutions_Retrieve_Next
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;

    ls : Link_to_Solution;
    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v(v'first));
    idx : natural32;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_Retrieve_Next ...");
    end if;
    if k = 0 then
      DoblDobl_Solutions_Container.Retrieve_Next_Initialize;
    else
      DoblDobl_Solutions_Container.Retrieve_Next(ls,idx);
      Assign(integer32(idx),a);
      if idx = 0
       then return 277;
       else Assign_Solution(ls,b,c);
      end if;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Retrieve_Next.");
      end if;
      return 277;
  end DoblDobl_Solutions_Retrieve_Next;

  function DoblDobl_Solutions_Move_Pointer
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    ind : natural32;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_Move_Pointer ...");
    end if;
    DoblDobl_Solutions_Container.Move_Current(ind);
    Assign(integer32(ind),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Move_Pointer.");
      end if;
      return 455;
  end DoblDobl_Solutions_Move_Pointer;

  function DoblDobl_Solutions_Current_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;

    ind,len : natural32;
    ls : Link_to_Solution;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_Current_Size ...");
    end if;
    DoblDobl_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      len := DoblDobl_Solution_Strings.Length(ls.all);
      Assign(integer32(len),b);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Current_Size.");
      end if;
      return 526;
  end DoblDobl_Solutions_Current_Size;

  function DoblDobl_Solutions_Current_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    len : constant natural32 := natural32(v_a(v_a'first));
    ind : natural32;
    ls : Link_to_Solution;
    sv : Standard_Integer_Vectors.Vector(1..integer32(len));

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_Current_String ...");
    end if;
    DoblDobl_Solutions_Container.Retrieve_Current(ls,ind);
    Assign(integer32(ind),a);
    if ind /= 0 then
      declare
        s : constant string := DoblDobl_Solution_Strings.Write(ls.all);
      begin
        sv := String_to_Integer_Vector(Pad_with_Spaces(len,s));
      end;
      Assign(sv,b);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Current_String.");
      end if;
      return 534;
  end DoblDobl_Solutions_Current_String;

  function DoblDobl_Solutions_Drop_by_Index
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant natural32 := natural32(v_a(v_a'first));
    use DoblDobl_Complex_Solutions;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    dropped : constant Solution_List := Solution_Drops.Drop(sols,ind);

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_Drop_by_Index ...");
    end if;
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(dropped);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Drop_by_Index.");
      end if;
      return 348;
  end DoblDobl_Solutions_Drop_by_Index;

  function DoblDobl_Solutions_Drop_by_Name
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
    use DoblDobl_Complex_Solutions;
    sols,dropped : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_Drop_by_Name ...");
    end if;
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    sols := DoblDobl_Solutions_Container.Retrieve;
    dropped := Solution_Drops.Drop(sols,ind);
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(dropped);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Drop_by_Name.");
      end if;
      return 349;
  end DoblDobl_Solutions_Drop_by_Name;

  function DoblDobl_Solutions_Make_Homogeneous
             ( vrblvl : integer32 := 0 ) return integer32 is

    sols : DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in standard_solutions_interface.");
      put_line("Standard_Solutions_Make_Homogeneous ...");
    end if;
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then Projective_Transformations.Projective_Transformation(sols);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Make_Homogeneous.");
      end if;
      return 895;
  end DoblDobl_Solutions_Make_Homogeneous;

  function DoblDobl_Solutions_Multi_Homogeneous
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    m : constant natural32 := natural32(v_a(v_a'first));
    sols : DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_Multi_Homogeneous ...");
    end if;
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then Multi_Projective_Transformations.Add_Ones(sols,m);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Multi_Homogeneous.");
      end if;
      return 911;
  end DoblDobl_Solutions_Multi_Homogeneous;

  function DoblDobl_Solutions_1Hom2Affine
             ( vrblvl : integer32 := 0 ) return integer32 is

    sols : DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_1Hom2Affine ...");
    end if;
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then Projective_Transformations.Affine_Transformation(sols);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_1Hom2Affine.");
      end if;
      return 899;
  end DoblDobl_Solutions_1Hom2Affine;

  function DoblDobl_Solutions_mHom2Affine
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    mhom : constant natural32 := natural32(v_a(v_a'first+1));
    idz : Standard_Natural_Vectors.Vector(1..integer32(nvr));
    sols : DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solutions_interface.");
      put_line("DoblDobl_Solutions_mHom2Affine ...");
    end if;
    Assign(nvr,b,idz);
    if not DoblDobl_Complex_Solutions.Is_Null(sols)
     then Multi_Projective_Transformations.Make_Affine(sols,mhom,idz);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_mHom2Affine.");
      end if;
      return 914;
  end DoblDobl_Solutions_mHom2Affine;

  function DoblDobl_Solutions_Tzero
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;

    sols : Solution_List := DoblDobl_Solutions_Container.Retrieve;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer32(0));

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_solution_interface.");
      put_line("DoblDobl_Solutions_Tzero ...");
    end if;
    if not Is_Null(sols)
     then Set_Continuation_Parameter(sols,zero);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_solutions_interface.");
        put_line("DoblDobl_Solutions_Tzero.");
      end if;
      return 876;
  end DoblDobl_Solutions_Tzero;

  function DoblDobl_Solutions_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_Solutions_interface.");
      put_line("DoblDobl_Solutions_Clear ...");
    end if;
    DoblDobl_Solutions_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_Solutions_interface.");
        put_line("DoblDobl_Solutions_Clear.");
      end if;
      return 347;
  end DoblDobl_Solutions_Clear;

end DoblDobl_Solutions_Interface;
