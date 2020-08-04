with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Symbol_Table;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Laur_Strings;
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
      return 570;
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
      return 571;
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
      return 572;
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
      return 573;
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
      return 574;
  end Multprec_LaurSys_Size;

  function Multprec_LaurSys_String_Save
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));
    nc : constant integer := integer(v_a(v_a'first));
    n : constant natural32 := natural32(v_a(v_a'first+1));
    k : constant integer32 := integer32(v_a(v_a'first+2));
    size : constant natural32 := natural32(v_a(v_a'first+3));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    p : Multprec_Complex_Laurentials.Poly;

  begin
    if vrblvl > 0 then
      put("-> in multprec_laursys_interface.");
      put_line("Multprec_LaurSys_String_Save ...");
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
    p := Multprec_Complex_Laur_Strings.Parse(n,size,s);
    Multprec_LaurSys_Container.Add_Poly(k,p);
    Multprec_Complex_Laurentials.Clear(p);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in multprec_Laursys_interface.");
        put_line("Multprec_LaurSys_String_Save.");
      end if;
      return 578;
  end Multprec_LaurSys_String_Save;

  function Multprec_LaurSys_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    equ : constant integer32 := integer32(v_a(v_a'first));
    p : constant Multprec_Complex_Laurentials.Poly
      := Multprec_LaurSys_Container.Retrieve_Poly(equ);
    sz : constant integer32
       := integer32(Multprec_Complex_Laur_Strings.Size_Limit(p));

  begin
    if vrblvl > 0 then
      put("-> in multprec_laursys_interface.");
      put_line("Multprec_LaurSys_String_Size ...");
    end if;
    Assign(sz,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in multprec_laursys_interface.");
        put_line("Multprec_LaurSys_String_Size.");
      end if;
      return 604;
  end Multprec_LaurSys_String_Size;

  function Multprec_LaurSys_String_Load 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    equ : constant integer32 := integer32(v_a(v_a'first));
    p : constant Multprec_Complex_Laurentials.Poly
      := Multprec_LaurSys_Container.Retrieve_Poly(equ);
    s : constant string := Multprec_Complex_Laur_Strings.Write(p);
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(s);
    slast : constant integer32 := integer32(s'last);

  begin
    if vrblvl > 0 then
      put("-> in multprec_laursys_interface.");
      put_line("Multprec_LaurSys_String_Load.");
    end if;
   -- put("Polynomial "); put(equ,1); put(" : "); put_line(s);
   -- put("#characters : "); put(s'last,1); new_line;
    Assign(slast,a);
    Assign(sv,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in multprec_laursys_interface.");
        put_line("Multprec_LaurSys_String_Load.");
      end if;
      return 579;
  end Multprec_LaurSys_String_Load;

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
      return 577;
  end Multprec_LaurSys_Clear;

end Multprec_LaurSys_Interface;
