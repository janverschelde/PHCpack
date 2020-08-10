with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_PolySys_Container;
with QuadDobl_Systems_Pool;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

package body QuadDobl_SysPool_Interface is

  function QuadDobl_SysPool_Initialize
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(v(v'first));

  begin
    if vrblvl > 0 then
      put("-> in quaddobl_syspool_interface.");
      put_line("QuadDobl_SysPool_Initialize ...");
    end if;
    QuadDobl_Systems_Pool.Initialize(n);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in quaddobl_syspool_interface.");
        put_line("QuadDobl_SysPool_Initialize.");
      end if;
      return 319;
  end QuadDobl_SysPool_Initialize;

  function QuadDobl_SysPool_from_Container
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    p : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := QuadDobl_PolySys_Container.Retrieve;

    use QuadDobl_Complex_Poly_Systems;

  begin
    if vrblvl > 0 then
      put("-> in quaddobl_syspool_interface.");
      put_line("QuadDobl_SysPool_from_Container ...");
    end if;
    if p /= null
     then QuadDobl_Systems_Pool.Create(k,p.all);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in quaddobl_syspool_interface.");
        put_line("QuadDobl_SysPool_from_Container.");
      end if;
      return 609;
  end QuadDobl_SysPool_from_Container;

  function QuadDobl_SysPool_into_Container
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    ind : constant integer32 := integer32(v_a(v_a'first));
    sys : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
        := QuadDobl_Systems_Pool.Retrieve(ind);

  begin
    if vrblvl > 0 then
      put("-> in quaddobl_syspool_interface.");
      put_line("QuadDobl_SysPool_into_Container ...");
    end if;
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(sys.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in quaddobl_syspool_interface.");
        put_line("QuadDobl_SysPool_into_Container.");
      end if;
      return 315;
  end QuadDobl_SysPool_into_Container;

  function QuadDobl_SysPool_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : constant natural32 := QuadDobl_Systems_Pool.Size;

  begin
    if vrblvl > 0 then
      put("-> in quaddobl_syspool_interface.");
      put_line("QuadDobl_SysPool_Size ...");
    end if;
    Assign(integer32(n),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in quaddobl_syspool_interface.");
        put_line("QuadDobl_SysPool_Size.");
      end if;
      return 317;
  end QuadDobl_SysPool_Size;

  function QuadDobl_SysPool_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in quaddobl_syspool_interface.");
      put_line("QuadDobl_SysPool_Clear ...");
    end if;
    QuadDobl_Systems_Pool.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in quaddobl_syspool_interface.");
        put_line("QuadDobl_SysPool_Clear.");
      end if;
      return 699;
  end QuadDobl_SysPool_Clear;

end QuadDobl_SysPool_Interface;
