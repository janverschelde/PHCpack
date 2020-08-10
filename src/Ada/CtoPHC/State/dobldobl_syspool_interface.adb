with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_PolySys_Container;
with DoblDobl_Systems_Pool;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

package body DoblDobl_SysPool_Interface is

  function DoblDobl_SysPool_Initialize
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(v(v'first));

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_syspool_interface.");
      put_line("DoblDobl_SysPool_Initialize ...");
    end if;
    DoblDobl_Systems_Pool.Initialize(n);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_syspool_interface.");
        put_line("DoblDobl_SysPool_Initialize.");
      end if;
      return 318;
  end DoblDobl_SysPool_Initialize;

  function DoblDobl_SysPool_from_Container
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    p : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := DoblDobl_PolySys_Container.Retrieve;

    use DoblDobl_Complex_Poly_Systems;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_syspool_interface.");
      put_line("DoblDobl_SysPool_from_Container ...");
    end if;
    if p /= null
     then DoblDobl_Systems_Pool.Create(k,p.all);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_syspool_interface.");
        put_line("DoblDobl_SysPool_from_Container.");
      end if;
      return 608;
  end DoblDobl_SysPool_from_Container;

  function DoblDobl_SysPool_into_Container
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    ind : constant integer32 := integer32(v_a(v_a'first));
    sys : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
        := DoblDobl_Systems_Pool.Retrieve(ind);

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_syspool_interface.");
      put_line("DoblDobl_SysPool_into_Container ...");
    end if;
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(sys.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_syspool_interface.");
        put_line("DoblDobl_SysPool_into_Container.");
      end if;
      return 314;
  end DoblDobl_SysPool_into_Container;

  function DoblDobl_SysPool_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : constant natural32 := DoblDobl_Systems_Pool.Size;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_syspool_interface.");
      put_line("DoblDobl_SysPool_Size ...");
    end if;
    Assign(integer32(n),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_syspool_interface.");
        put_line("DoblDobl_SysPool_Size.");
      end if;
      return 316;
  end DoblDobl_SysPool_Size;

  function DoblDobl_SysPool_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_syspool_interface.");
      put_line("DoblDobl_SysPool_Clear ...");
    end if;
    DoblDobl_Systems_Pool.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_syspool_interface.");
        put_line("DoblDobl_SysPool_Clear.");
      end if;
      return 698;
  end DoblDobl_SysPool_Clear;

end DoblDobl_SysPool_Interface;
