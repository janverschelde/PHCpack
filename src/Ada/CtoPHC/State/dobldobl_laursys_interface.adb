with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Symbol_Table;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;  use DoblDobl_Complex_Laur_Systems_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with PHCpack_Operations;
with DoblDobl_LaurSys_Container;

package body DoblDobl_LaurSys_Interface is

  function DoblDobl_LaurSys_Read
             ( vrblvl : integer32 := 0 ) return integer32 is

    lp : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put_line("-> in DoblDobl_LaurSys_interface.DoblDobl_LaurSys_Read ...");
    end if;
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    DoblDobl_LaurSys_Container.Initialize(lp.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_LaurSys_interface.");
        put_line("DoblDobl_LaurSys_Read.");
      end if;
      return 110;
  end DoblDobl_LaurSys_Read;

  function DoblDobl_LaurSys_Write
             ( vrblvl : in integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Laur_Systems;

    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    nvr : natural32;

  begin
    if vrblvl > 0 then
      put_line("-> in DoblDobl_LaurSys_interface.DoblDobl_LaurSys_Write ...");
    end if;
    if lp /= null then
      nvr := Number_of_Unknowns(lp(lp'first));
      if PHCpack_Operations.file_okay then
        if integer32(nvr) = lp'last then
          put(PHCpack_Operations.output_file,natural32(lp'last),lp.all);
        else
          put(PHCpack_Operations.output_file,natural32(lp'last),nvr,lp.all);
        end if;
      elsif integer32(nvr) = lp'last then
        put(standard_output,natural32(lp'last),lp.all);
      else
        put(standard_output,natural32(lp'last),nvr,lp.all);
      end if;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_LaurSys_interface.");
        put_line("DoblDobl_LaurSys_Write.");
      end if;
      return 111;
  end DoblDobl_LaurSys_Write;

  function DoblDobl_LaurSys_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_LaurSys_interface.");
      put_line("-> DoblDobl_LaurSys_Get_Dimension ...");
    end if;
    Assign(integer32(DoblDobl_LaurSys_Container.Dimension),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_LaurSys_interface.");
        put_line("DoblDobl_LaurSys_Get_Dimension.");
      end if;
      return 112;
  end DoblDobl_LaurSys_Get_Dimension;

  function DoblDobl_LaurSys_Set_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant integer32 := integer32(v(v'first));

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_LaurSys_interface.");
      put_line("DoblDobl_LaurSys_Set_Dimension ...");
    end if;
    DoblDobl_LaurSys_Container.Initialize(n);
    Symbol_Table.Init(natural32(n));
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_LaurSys_interface.");
        put_line("DoblDobl_LaurSys_Set_Dimension.");
      end if;
      return 113;
  end DoblDobl_LaurSys_Set_Dimension;

  function DoblDobl_LaurSys_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    i : constant integer32 := integer32(v(v'first+1));

  begin
    if vrblvl > 0 then
      put_line("-> in DoblDobl_LaurSys_interface.DoblDobl_LaurSys_Size ...");
    end if;
    Assign(integer32(DoblDobl_LaurSys_Container.Number_of_Terms(i)),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_LaurSys_interface.");
        put_line("DoblDobl_LaurSys_Size");
      end if;
      return 114;
  end DoblDobl_LaurSys_Size;

  function DoblDobl_LaurSys_Get_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array(0..2)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    i : constant integer32 := integer32(v(1));
    j : constant natural32 := natural32(v(2));
    t : constant DoblDobl_Complex_Laurentials.Term
      := DoblDobl_LaurSys_Container.Retrieve_Term(i,j);

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_LaurSys_interface.");
      put_line("DoblDobl_LaurSys_Get_Term ...");
    end if;
    Assign(t.cf,c);
    Assign(t.dg.all,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_LaurSys_interface.");
        put_line("DoblDobl_LaurSys_Get_Term");
      end if;
      return 115;
  end DoblDobl_LaurSys_Get_Term;

  function DoblDobl_LaurSys_Add_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(v(0));
    i : constant integer32 := integer32(v(1));
    e : Standard_Integer_Vectors.Vector(1..n);
    t : DoblDobl_Complex_Laurentials.Term;

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_LaurSys_interface.");
      put_line("DoblDobl_LaurSys_Add_Term ...");
    end if;
    Assign(c,t.cf);
    Assign(natural32(n),b,e);
    t.dg := new Standard_Integer_Vectors.Vector'(e);
    DoblDobl_LaurSys_Container.Add_Term(i,t);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_LaurSys_interface.");
        put_line("DoblDobl_LaurSys_Add_Term");
      end if;
      return 116;
  end DoblDobl_LaurSys_Add_Term;

  function DoblDobl_LaurSys_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put_line("-> in DoblDobl_LaurSys_interface.DoblDobl_LaurSys_Clear ...");
    end if;
    DoblDobl_LaurSys_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_LaurSys_interface.");
        put_line("DoblDobl_LaurSys_Clear.");
      end if;
      return 117;
  end DoblDobl_LaurSys_Clear;

end DoblDobl_LaurSys_Interface;
