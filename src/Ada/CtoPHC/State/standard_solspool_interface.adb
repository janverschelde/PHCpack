with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Solutions_Pool;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;

package body Standard_SolsPool_Interface is

  function Standard_SolsPool_Initialize
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(v(v'first));

  begin
    if vrblvl > 0 then
      put("-> in standard_solspool_interface.");
      put_line("Standard_SolsPool_Initialize ...");
    end if;
    Solutions_Pool.Initialize(n);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solspool_interface.");
        put_line("Standard_SolsPool_Initialize.");
      end if;
      return 320;
  end Standard_SolsPool_Initialize;

  function Standard_SolsPool_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : constant natural32 := Solutions_Pool.Size;

  begin
    if vrblvl > 0 then
      put("-> in standard_solspool_interface.");
      put_line("Standard_SolsPool_Size ...");
    end if;
    Assign(integer32(n),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solspool_interface.");
        put_line("Standard_SolsPool_Initialize.");
      end if;
      return 321;
  end Standard_SolsPool_Size;

  function Standard_SolsPool_Length
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    res : constant natural32 := Solutions_Pool.Length(k);

  begin
    if vrblvl > 0 then
      put("-> in standard_solspool_interface.");
      put_line("Standard_SolsPool_Length ...");
    end if;
    Assign(integer32(res),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solspool_interface.");
        put_line("Standard_SolsPool_Length.");
      end if;
      return 322;
  end Standard_SolsPool_Length;

  function Standard_SolsPool_Dimension
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    res : constant natural32 := Solutions_Pool.Dimension(k);

  begin
    if vrblvl > 0 then
      put("-> in standard_solspool_interface.");
      put_line("Standard_SolsPool_Dimension ...");
    end if;
    Assign(integer32(res),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solspool_interface.");
        put_line("Standard_SolsPool_Dimension.");
      end if;
      return 323;
  end Standard_SolsPool_Dimension;

  function Standard_SolsPool_Add
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    ls : constant Link_to_Solution := Convert_to_Solution(b,c);

  begin
    if vrblvl > 0 then
      put("-> in standard_solspool_interface.");
      put_line("Standard_SolsPool_Add ...");
    end if;
   -- put("Thread "); put(k,1); put_line(" appends a solution ...");
    Solutions_Pool.Append(k,ls);
   -- delay(1.0);
   -- put("Number of solutions in pool "); put(k,1);
   -- put(" is "); put(Solutions_Pool.Length(k),1); new_line;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solspool_interface.");
        put_line("Standard_SolsPool_Add.");
      end if;
      return 324;
  end Standard_SolsPool_Add;

  function Standard_SolsPool_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant integer32 := integer32(v_a(v_a'first));
    i : constant natural32 := natural32(v_a(v_a'first+1));
    ls : Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in standard_solspool_interface.");
      put_line("Standard_SolsPool_Get ...");
    end if;
    Solutions_Pool.Retrieve(k,i,ls,fail);
    if fail then
      return 325;
    else
      Assign_Solution(ls,b,c);
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_solspool_interface.");
        put_line("Standard_SolsPool_Get.");
      end if;
      return 325;
  end Standard_SolsPool_Get;

end Standard_SolsPool_Interface;
