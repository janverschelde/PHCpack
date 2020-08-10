with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with PHCpack_Operations;                use PHCpack_Operations;
with Standard_PolySys_Container;
with Standard_Systems_Pool;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

package body Standard_SysPool_Interface is

  function Standard_SysPool_Initialize
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(v(v'first));

  begin
    if vrblvl > 0 then
      put("-> in standard_syspool_interface.");
      put_line("Standard_SysPool_Initialize ...");
    end if;
    Standard_Systems_Pool.Initialize(n);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_syspool_interface.");
        put_line("Standard_SysPool_Initialize.");
      end if;
      return 300;
  end Standard_SysPool_Initialize;

  function Standard_SysPool_Read
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    
  begin
    if vrblvl > 0 then
      put("-> in standard_syspool_interface.");
      put_line("Standard_SysPool_Read ...");
    end if;
    new_line;
    put_line("Reading a polynomial system ...");
    get(p);
    Standard_Systems_Pool.Create(k,p.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_syspool_interface.");
        put_line("Standard_SysPool_Read.");
      end if;
      return 302;
  end Standard_SysPool_Read;

  function Standard_SysPool_Write
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    p : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Standard_Systems_Pool.Retrieve(k);

    use Standard_Complex_Poly_Systems;

  begin
    if vrblvl > 0 then
      put("-> in standard_syspool_interface.");
      put_line("Standard_SysPool_Write ...");
    end if;
    if p /= null then
      if PHCpack_Operations.Is_File_Defined
       then put(PHCpack_Operations.output_file,natural32(p'last),p.all);
       else put(standard_output,natural32(p'last),p.all);
      end if;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_syspool_interface.");
        put_line("Standard_SysPool_Write.");
      end if;
      return 303;
  end Standard_SysPool_Write;

  function Standard_SysPool_from_Container
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    p : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Standard_PolySys_Container.Retrieve;

    use Standard_Complex_Poly_Systems;

  begin
    if vrblvl > 0 then
      put("-> in standard_syspool_interface.");
      put_line("Standard_SysPool_from_Container ...");
    end if;
    if p /= null
     then Standard_Systems_Pool.Create(k,p.all);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_syspool_interface.");
        put_line("Standard_SysPool_from_Container.");
      end if;
      return 304;
  end Standard_SysPool_from_Container;

  function Standard_SysPool_into_Container
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    ind : constant integer32 := integer32(v_a(v_a'first));
    sys : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
        := Standard_Systems_Pool.Retrieve(ind);

  begin
    if vrblvl > 0 then
      put("-> in standard_syspool_interface.");
      put_line("Standard_SysPool_into_Container ...");
    end if;
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(sys.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_syspool_interface.");
        put_line("Standard_SysPool_into_Container.");
      end if;
      return 313;
  end Standard_SysPool_into_Container;

  function Standard_SysPool_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : constant natural32 := Standard_Systems_Pool.Size;

  begin
    if vrblvl > 0 then
      put("-> in standard_syspool_interface.");
      put_line("Standard_SysPool_Size ...");
    end if;
    Assign(integer32(n),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_syspool_interface.");
        put_line("Standard_SysPool_Size.");
      end if;
      return 301;
  end Standard_SysPool_Size;

  function Standard_SysPool_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in standard_syspool_interface.");
      put_line("Standard_SysPool_Clear ...");
    end if;
    Standard_Systems_Pool.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_syspool_interface.");
        put_line("Standard_SysPool_Clear.");
      end if;
      return 697;
  end Standard_SysPool_Clear;

end Standard_SysPool_Interface;
