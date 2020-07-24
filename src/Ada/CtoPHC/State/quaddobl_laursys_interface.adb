with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Symbol_Table;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;  use QuadDobl_Complex_Laur_Systems_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with PHCpack_Operations;
with QuadDobl_LaurSys_Container;

package body QuadDobl_LaurSys_Interface is

  function QuadDobl_LaurSys_Read
             ( vrblvl : integer32 := 0 ) return integer32 is

    lp : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put_line("-> in QuadDobl_LaurSys_interface.QuadDobl_LaurSys_Read ...");
    end if;
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    QuadDobl_LaurSys_Container.Initialize(lp.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_LaurSys_interface.");
        put_line("QuadDobl_LaurSys_Read.");
      end if;
      return 120;
  end QuadDobl_LaurSys_Read;

  function QuadDobl_LaurSys_Write
             ( vrblvl : in integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Systems;

    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    nvr : natural32;

  begin
    if vrblvl > 0 then
      put_line("-> in QuadDobl_LaurSys_interface.QuadDobl_LaurSys_Write ...");
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
        put("Exception raised in QuadDobl_LaurSys_interface.");
        put_line("QuadDobl_LaurSys_Write.");
      end if;
      return 121;
  end QuadDobl_LaurSys_Write;

  function QuadDobl_LaurSys_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_LaurSys_interface.");
      put_line("-> QuadDobl_LaurSys_Get_Dimension ...");
    end if;
    Assign(integer32(QuadDobl_LaurSys_Container.Dimension),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_LaurSys_interface.");
        put_line("QuadDobl_LaurSys_Get_Dimension.");
      end if;
      return 122;
  end QuadDobl_LaurSys_Get_Dimension;

  function QuadDobl_LaurSys_Set_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant integer32 := integer32(v(v'first));

  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_LaurSys_interface.");
      put_line("QuadDobl_LaurSys_Set_Dimension ...");
    end if;
    QuadDobl_LaurSys_Container.Initialize(n);
    Symbol_Table.Init(natural32(n));
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_LaurSys_interface.");
        put_line("QuadDobl_LaurSys_Set_Dimension.");
      end if;
      return 123;
  end QuadDobl_LaurSys_Set_Dimension;

  function QuadDobl_LaurSys_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    i : constant integer32 := integer32(v(v'first+1));

  begin
    if vrblvl > 0 then
      put_line("-> in QuadDobl_LaurSys_interface.QuadDobl_LaurSys_Size ...");
    end if;
    Assign(integer32(QuadDobl_LaurSys_Container.Number_of_Terms(i)),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_LaurSys_interface.");
        put_line("QuadDobl_LaurSys_Size");
      end if;
      return 124;
  end QuadDobl_LaurSys_Size;

  function QuadDobl_LaurSys_Get_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array(0..2)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    i : constant integer32 := integer32(v(1));
    j : constant natural32 := natural32(v(2));
    t : constant QuadDobl_Complex_Laurentials.Term
      := QuadDobl_LaurSys_Container.Retrieve_Term(i,j);

  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_LaurSys_interface.");
      put_line("QuadDobl_LaurSys_Get_Term ...");
    end if;
    Assign(t.cf,c);
    Assign(t.dg.all,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_LaurSys_interface.");
        put_line("QuadDobl_LaurSys_Get_Term");
      end if;
      return 125;
  end QuadDobl_LaurSys_Get_Term;

  function QuadDobl_LaurSys_Add_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(v(0));
    i : constant integer32 := integer32(v(1));
    e : Standard_Integer_Vectors.Vector(1..n);
    t : QuadDobl_Complex_Laurentials.Term;

  begin
    if vrblvl > 0 then
      put("-> in QuadDobl_LaurSys_interface.");
      put_line("QuadDobl_LaurSys_Add_Term ...");
    end if;
    Assign(c,t.cf);
    Assign(natural32(n),b,e);
    t.dg := new Standard_Integer_Vectors.Vector'(e);
    QuadDobl_LaurSys_Container.Add_Term(i,t);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_LaurSys_interface.");
        put_line("QuadDobl_LaurSys_Add_Term");
      end if;
      return 126;
  end QuadDobl_LaurSys_Add_Term;

  function QuadDobl_LaurSys_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put_line("-> in QuadDobl_LaurSys_interface.QuadDobl_LaurSys_Clear ...");
    end if;
    QuadDobl_LaurSys_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in QuadDobl_LaurSys_interface.");
        put_line("QuadDobl_LaurSys_Clear.");
      end if;
      return 127;
  end QuadDobl_LaurSys_Clear;

end QuadDobl_LaurSys_Interface;
