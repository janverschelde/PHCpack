with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;
with PHCpack_Operations;
with QuadDobl_Solutions_Container;

package body QuadDobl_Solutions_Interface is

  function QuadDobl_Solutions_Read 
             ( vrblvl : integer32 := 0 ) return integer32 is

    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_Solutions_interface.");
      put_line("QuadDobl_Solutions_Read ...");
    end if;
    QuadDobl_Complex_Solutions_io.Read(sols);
    QuadDobl_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_Solutions_interface.");
        put_line("QuadDobl_Solutions_Read.");
      end if;
      return 390;
  end QuadDobl_Solutions_Read;

  function QuadDobl_Solutions_Write
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Solutions_io;

    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_Solutions_interface.");
      put_line("QuadDobl_Solutions_Write ...");
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
        put("Exception raised in QuadDobl_Solutions_interface.");
        put_line("QuadDobl_Solutions_Write.");
      end if;
      return 391;
  end QuadDobl_Solutions_Write;

  function QuadDobl_Solutions_Size 
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_Solutions_interface.");
      put_line("QuadDobl_Solutions_Size ...");
    end if;
    Assign(integer32(QuadDobl_Solutions_Container.Length),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_Solutions_interface.");
        put_line("QuadDobl_Solutions_Size.");
      end if;
      return 392;
  end QuadDobl_Solutions_Size;

  function QuadDobl_Solutions_Dimension
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_Solutions_interface.");
      put_line("QuadDobl_Solutions_Dimension ...");
    end if;
    Assign(integer32(QuadDobl_Solutions_Container.Dimension),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_Solutions_interface.");
        put_line("QuadDobl_Solutions_Dimension.");
      end if;
      return 393;
  end QuadDobl_Solutions_Dimension;

  function QuadDobl_Solutions_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
           
    use QuadDobl_Complex_Solutions;

    ls : Link_to_Solution;
    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_Solutions_interface.");
      put_line("QuadDobl_Solutions_Get ...");
    end if;
    QuadDobl_Solutions_Container.Retrieve(k,ls,fail);
    if fail
     then return 394;
     else Assign_Solution(ls,b,c); return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_Solutions_interface.");
        put_line("QuadDobl_Solutions_Get.");
      end if;
      return 394;
  end QuadDobl_Solutions_Get;

  function QuadDobl_Solutions_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
           
    use QuadDobl_Complex_Solutions;

    ls : Link_to_Solution := Convert_to_Solution(b,c);
    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v(v'first));
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_Solutions_interface.");
      put_line("QuadDobl_Solutions_Read ...");
    end if;
    QuadDobl_Solutions_Container.Replace(k,ls.all,fail);
    Clear(ls);
    if fail
     then return 395;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_Solutions_interface.");
        put_line("QuadDobl_Solutions_Read.");
      end if;
      return 395;
  end QuadDobl_Solutions_Update;

  function QuadDobl_Solutions_Add
             ( b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
           
    use QuadDobl_Complex_Solutions;

    ls : constant Link_to_Solution := Convert_to_Solution(b,c);

  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_Solutions_interface.");
      put_line("QuadDobl_Solutions_Add ...");
    end if;
    QuadDobl_Solutions_Container.Append(ls);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_Solutions_interface.");
        put_line("QuadDobl_Solutions_Add.");
      end if;
      return 396;
  end QuadDobl_Solutions_Add;

  function QuadDobl_Solutions_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_Solutions_interface.");
      put_line("QuadDobl_Solutions_Clear ...");
    end if;
    QuadDobl_Solutions_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_Solutions_interface.");
        put_line("QuadDobl_Solutions_Clear.");
      end if;
      return 397;
  end QuadDobl_Solutions_Clear;

end QuadDobl_Solutions_Interface;
