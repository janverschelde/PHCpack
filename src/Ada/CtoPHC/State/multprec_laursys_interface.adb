with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Symbol_Table;
with Multprec_Complex_Laur_Systems;
with Multprec_Complex_Laur_Systems_io;  use Multprec_Complex_Laur_Systems_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with PHCpack_Operations;
with Multprec_LaurSys_Container;

package body Multprec_LaurSys_Interface is

  function Multprec_LaurSys_Read
             ( vrblvl : integer32 := 0 ) return integer32 is

    lp : Multprec_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put_line("-> in multprec_laursys_interface.Multprec_LaurSys_Read ...");
    end if;
    new_line;
    put_line("Reading a Laurnomial system ...");
    get(lp);
    Multprec_LaurSys_Container.Initialize(lp.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in multprec_laursys_interface.");
        put_line("Multprec_LaurSys_Read.");
      end if;
      return 130;
  end Multprec_LaurSys_Read;

  function Multprec_LaurSys_Write
             ( vrblvl : in integer32 := 0 ) return integer32 is

   -- use Multprec_Complex_Laurnomials;
    use Multprec_Complex_Laur_Systems;

    lp : constant Link_to_Laur_Sys := Multprec_LaurSys_Container.Retrieve;
   -- nvr : natural32;

  begin
    if vrblvl > 0 then
      put_line("-> in multprec_laursys_interface.multprec_LaurSys_Write ...");
    end if;
    if lp /= null then
      if PHCpack_Operations.file_okay
      -- then put(PHCpack_Operations.output_file,lp'last,lp.all);
       then put(PHCpack_Operations.output_file,lp.all);
      -- else put(standard_output,lp'last,lp.all);
       else put(standard_output,lp.all);
      end if;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in multprec_laursys_interface.");
        put_line("Multprec_LaurSys_Write.");
      end if;
      return 131;
  end Multprec_LaurSys_Write;

  function Multprec_LaurSys_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in multprec_laursys_interface.");
      put_line("-> Multprec_LaurSys_Get_Dimension ...");
    end if;
    Assign(integer32(Multprec_LaurSys_Container.Dimension),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in multprec_laursys_interface.");
        put_line("Multprec_LaurSys_Get_Dimension.");
      end if;
      return 132;
  end Multprec_LaurSys_Get_Dimension;

  function Multprec_LaurSys_Set_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant integer32 := integer32(v(v'first));

  begin
    if vrblvl > 0 then
      put("-> in multprec_laursys_interface.");
      put_line("Multprec_LaurSys_Set_Dimension ...");
    end if;
    Multprec_LaurSys_Container.Initialize(n);
    Symbol_Table.Init(natural32(n));
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in multprec_laursys_interface.");
        put_line("Multprec_LaurSys_Set_Dimension.");
      end if;
      return 133;
  end Multprec_LaurSys_Set_Dimension;

  function Multprec_LaurSys_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    i : constant integer32 := integer32(v(v'first+1));

  begin
    if vrblvl > 0 then
      put_line("-> in multprec_laursys_interface.Multprec_LaurSys_Size ...");
    end if;
    Assign(integer32(Multprec_LaurSys_Container.Number_of_Terms(i)),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in multprec_laursys_interface.");
        put_line("Multprec_LaurSys_Size");
      end if;
      return 134;
  end Multprec_LaurSys_Size;

  function Multprec_LaurSys_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put_line("-> in multprec_laursys_interface.Multprec_LaurSys_Clear ...");
    end if;
    Multprec_LaurSys_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in multprec_laursys_interface.");
        put_line("Multprec_LaurSys_Clear.");
      end if;
      return 137;
  end Multprec_LaurSys_Clear;

end Multprec_LaurSys_Interface;
