with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Symbol_Table;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Strings;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;  use QuadDobl_Complex_Laur_Systems_io;
with Polynomial_Drops;
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

  function QuadDobl_LaurSys_String_Save
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    nc : constant integer := integer(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    k : constant integer32 := integer32(v_a(v_a'first+2));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    p : QuadDobl_Complex_Laurentials.Poly;

  begin
    if vrblvl > 0 then
      put("-> in quaddobl_laursys_interface.");
      put_line("QuadDobl_LaurSys_String_Save ...");
    end if;
   -- put("Polynomial "); put(k,1);
   -- put(" given as string of "); put(nc,1);
   -- put_line(" characters.");
   -- put("The string : "); put_line(s);
    if Symbol_Table.Empty then
      Symbol_Table.Init(n);
    elsif Symbol_Table.Maximal_Size < n then
      Symbol_Table.Clear;
      Symbol_Table.Init(n);
    end if;
    p := QuadDobl_Complex_Laur_Strings.Parse(n,s);
    QuadDobl_LaurSys_Container.Add_Poly(k,p);
    QuadDobl_Complex_Laurentials.Clear(p);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in quaddobl_Laursys_interface.");
        put_line("QuadDobl_LaurSys_String_Save.");
      end if;
      return 568;
  end QuadDobl_LaurSys_String_Save;

  function QuadDobl_LaurSys_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    equ : constant integer32 := integer32(v_a(v_a'first));
    p : constant QuadDobl_Complex_Laurentials.Poly
      := QuadDobl_LaurSys_Container.Retrieve_Poly(equ);
    sz : constant integer32
       := integer32(QuadDobl_Complex_Laur_Strings.Size_Limit(p));

  begin
    if vrblvl > 0 then
      put("-> in quaddobl_laursys_interface.");
      put_line("QuadDobl_LaurSys_String_Size ...");
    end if;
    Assign(sz,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in quaddobl_laursys_interface.");
        put_line("QuadDobl_LaurSys_String_Size.");
      end if;
      return 606;
  end QuadDobl_LaurSys_String_Size;

  function QuadDobl_LaurSys_String_Load 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    equ : constant integer32 := integer32(v_a(v_a'first));
    p : constant QuadDobl_Complex_Laurentials.Poly
      := QuadDobl_LaurSys_Container.Retrieve_Poly(equ);
    s : constant string := QuadDobl_Complex_Laur_Strings.Write(p);
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(s);
    slast : constant integer32 := integer32(s'last);

  begin
    if vrblvl > 0 then
      put("-> in quadDobl_laursys_interface.");
      put_line("QuadDobl_LaurSys_String_Load.");
    end if;
   -- put("Polynomial "); put(equ,1); put(" : "); put_line(s);
   -- put("#characters : "); put(s'last,1); new_line;
    Assign(slast,a);
    Assign(sv,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in quaddobl_laursys_interface.");
        put_line("QuadDobl_LaurSys_String_Load.");
      end if;
      return 569;
  end QuadDobl_LaurSys_String_Load;

  function QuadDobl_LaurSys_Drop_by_Index
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant integer32 := integer32(v_a(v_a'first));
    use QuadDobl_Complex_Laur_Systems;
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    dropped : constant Laur_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,ind);

  begin
    if vrblvl > 0 then
      put("-> in quaddobl_laursys_interface.");
      put_line("QuadDobl_LaurSys_Drop_by_Index ...");
    end if;
    QuadDobl_LaurSys_Container.Clear;
    QuadDobl_LaurSys_Container.Initialize(dropped);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in quaddobl_laursys_interface.");
        put_line("QuadDobl_LaurSys_Drop_by_Index.");
      end if;
      return 830;
  end QuadDobl_LaurSys_Drop_by_Index;

  function QuadDobl_LaurSys_Drop_by_Name
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);
    sb : Symbol_Table.Symbol;
    ind : natural32;
    use Quaddobl_Complex_Laur_Systems;
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    dropped : Laur_Sys(lp'range);

  begin
    if vrblvl > 0 then
      put("-> in quaddobl_laursys_interface.");
      put_line("QuadDobl_LaurSys_Drop_by_Name ...");
    end if;
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    dropped := Polynomial_Drops.Drop(lp.all,integer32(ind));
    QuadDobl_LaurSys_Container.Clear;
    QuadDobl_LaurSys_Container.Initialize(dropped);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in quaddobl_laursys_interface.");
        put_line("QuadDobl_LaurSys_Drop_by_Name.");
      end if;
      return 833;
  end QuadDobl_LaurSys_Drop_by_Name;

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
