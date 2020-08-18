with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Symbol_Table;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Strings;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Random_Polynomials;
with Polynomial_Drops;
with Homogenization;
with Projective_Transformations;
with Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;
with Multi_Projective_Transformations;
with Affine_Transformations;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with PHCpack_Operations;
with PHCpack_Operations_io;
with DoblDobl_PolySys_Container;

package body DoblDobl_PolySys_Interface is

  function DoblDobl_PolySys_Read
             ( vrblvl : integer32 := 0 ) return integer32 is

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put_line("-> in DoblDobl_PolySys_Interface.DoblDobl_PolySys_Read ...");
    end if;
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    DoblDobl_PolySys_Container.Initialize(lp.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_PolySys_Interface.");
        put_line("DoblDobl_PolySys_Read.");
      end if;
      return 200;
  end DoblDobl_PolySys_Read;

  function DoblDobl_PolySys_Read_from_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_PolySys_Interface.");
      put_line("DoblDobl_PolySys_Read_from_File ...");
    end if;
   -- new_line;
   -- put_line("Opening the file with name " & sv & " ...");
    declare
      file : file_type;
      p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    begin
      Open(file,in_file,sv);
      get(file,p);
      DoblDobl_PolySys_Container.Clear;
      DoblDobl_PolySys_Container.Initialize(p.all);
      exception 
        when NAME_ERROR =>
          put_line("File with name " & sv & " could not be found!");
          return 541;
        when USE_ERROR =>
          put_line("File with name " & sv & " could not be read!");
          return 541;
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_PolySys_Interface.");
        put_line("DoblDobl_PolySys_Read_from_File.");
      end if;
      return 541;
  end DoblDobl_PolySys_Read_from_File;

  function DoblDobl_PolySys_Write
             ( vrblvl : in integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    nvr : natural32;

  begin
    if vrblvl > 0 then
      put_line("-> in DoblDobl_PolySys_Interface.DoblDobl_PolySys_Write ...");
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
        put("Exception raised in DoblDobl_PolySys_Interface.");
        put_line("DoblDobl_PolySys_Write.");
      end if;
      return 201;
  end DoblDobl_PolySys_Write;

  function DoblDobl_PolySys_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_PolySys_Interface.");
      put_line("DoblDobl_PolySys_Get_Dimension ...");
    end if;
    Assign(integer32(DoblDobl_PolySys_Container.Dimension),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_PolySys_Interface.");
        put_line("DoblDobl_PolySys_Get_Dimension.");
      end if;
      return 202;
  end DoblDobl_PolySys_Get_Dimension;

  function DoblDobl_PolySys_Set_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant integer32 := integer32(v(v'first));

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_PolySys_Interface.");
      put_line("DoblDobl_PolySys_Set_Dimension ...");
    end if;
    DoblDobl_PolySys_Container.Initialize(n);
    Symbol_Table.Init(natural32(n));
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_PolySys_Interface.");
        put_line("DoblDobl_PolySys_Set_Dimension.");
      end if;
      return 203;
  end DoblDobl_PolySys_Set_Dimension;

  function DoblDobl_PolySys_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    i : constant integer32 := integer32(v(v'first+1));

  begin
    if vrblvl > 0 then
      put_line("-> in DoblDobl_PolySys_Interface.DoblDobl_PolySys_Size ...");
    end if;
    Assign(integer32(DoblDobl_PolySys_Container.Number_of_Terms(i)),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_PolySys_Interface.");
        put_line("DoblDobl_PolySys_Size");
      end if;
      return 204;
  end DoblDobl_PolySys_Size;

  function DoblDobl_PolySys_Degree
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    equ : constant integer32 := integer32(v_a(v_a'first));
    deg : constant integer32 := DoblDobl_PolySys_Container.Degree(equ);

  begin
    if vrblvl > 0 then
      put_line("-> in dobldobl_polysys_interface.DoblDobl_PolySys_Degree ...");
    end if;
    Assign(deg,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_Degree");
      end if;
      return 209;
  end DoblDobl_PolySys_Degree;

  function DoblDobl_PolySys_Get_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array(0..2)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    i : constant integer32 := integer32(v(1));
    j : constant natural32 := natural32(v(2));
    t : constant DoblDobl_Complex_Polynomials.Term
      := DoblDobl_PolySys_Container.Retrieve_Term(i,j);

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_PolySys_Interface.");
      put_line("DoblDobl_PolySys_Get_Term ...");
    end if;
    Assign(t.cf,c);
    Assign(t.dg.all,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_PolySys_Interface.");
        put_line("DoblDobl_PolySys_Get_Term");
      end if;
      return 205;
  end DoblDobl_PolySys_Get_Term;

  function DoblDobl_PolySys_Add_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(v(0));
    i : constant integer32 := integer32(v(1));
    e : Standard_Natural_Vectors.Vector(1..n);
    t : DoblDobl_Complex_Polynomials.Term;

  begin
    if vrblvl > 0 then
      put("-> in DoblDobl_PolySys_Interface.");
      put_line("DoblDobl_PolySys_Add_Term ...");
    end if;
    Assign(c,t.cf);
    Assign(natural32(n),b,e);
    t.dg := new Standard_Natural_Vectors.Vector'(e);
    DoblDobl_PolySys_Container.Add_Term(i,t);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_PolySys_Interface.");
        put_line("DoblDobl_PolySys_Add_Term");
      end if;
      return 206;
  end DoblDobl_PolySys_Add_Term;

  function DoblDobl_PolySys_String_Save
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
    p : DoblDobl_Complex_Polynomials.Poly;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_String_Save.");
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
    p := DoblDobl_Complex_Poly_Strings.Parse(n,s);
    DoblDobl_PolySys_Container.Add_Poly(k,p);
    DoblDobl_Complex_Polynomials.Clear(p);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_String_Save.");
      end if;
      return 338;
  end DoblDobl_PolySys_String_Save;

  function DoblDobl_PolySys_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    equ : constant integer32 := integer32(v_a(v_a'first));
    p : constant DoblDobl_Complex_Polynomials.Poly
      := DoblDobl_PolySys_Container.Retrieve_Poly(equ);
    sz : constant integer32
       := integer32(DoblDobl_Complex_Poly_Strings.Size_Limit(p));

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_String_Size ...");
    end if;
    Assign(sz,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_String_Size.");
      end if;
      return 601;
  end DoblDobl_PolySys_String_Size;

  function DoblDobl_PolySys_String_Load 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    equ : constant integer32 := integer32(v_a(v_a'first));
    p : constant DoblDobl_Complex_Polynomials.Poly
      := DoblDobl_PolySys_Container.Retrieve_Poly(equ);
    s : constant string := DoblDobl_Complex_Poly_Strings.Write(p);
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(s);
    slast : constant integer32 := integer32(s'last);

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_String_Load.");
    end if;
   -- put("Polynomial "); put(equ,1); put(" : "); put_line(s);
   -- put("#characters : "); put(s'last,1); new_line;
    Assign(slast,a);
    Assign(sv,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_String_Load.");
      end if;
      return 106;
  end DoblDobl_PolySys_String_Load;

  function DoblDobl_PolySys_Drop_by_Index
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant integer32 := integer32(v_a(v_a'first));
    use DoblDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    dropped : constant Poly_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,ind);

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_Drop_by_Index ...");
    end if;
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(dropped);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_Drop_by_Index.");
      end if;
      return 307;
  end DoblDobl_PolySys_Drop_by_Index;

  function DoblDobl_PolySys_Drop_by_Name
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
    use DoblDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    dropped : Poly_Sys(lp'range);

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_Drop_by_Name ...");
    end if;
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    dropped := Polynomial_Drops.Drop(lp.all,integer32(ind));
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(dropped);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_Drop_by_Name.");
      end if;
      return 310;
  end DoblDobl_PolySys_Drop_by_Name;

  function DoblDobl_PolySys_Random_System
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    neq : constant integer32 := integer32(v_a(v_a'first+1));
    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..neq);
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(3));
    m : constant natural32 := natural32(v_b(v_b'first));
    d : constant natural32 := natural32(v_b(v_b'first+1));
    c : constant natural32 := natural32(v_b(v_b'first+2));

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_Random_System ...");
    end if;
    for i in p'range loop
      if m = 0 then
        p(i) := DoblDobl_Random_Polynomials.Random_Dense_Poly(nvr,d,c);
      else
        p(i) := DoblDobl_Random_Polynomials.Random_Sparse_Poly(nvr,d,m,c);
      end if;
    end loop;
    DoblDobl_PolySys_Container.Clear; 
    DoblDobl_PolySys_Container.Initialize(p); 
   -- must initialize the symbol table with actual symbols for printing
    Symbol_Table.Init(Symbol_Table.Standard_Symbols(integer32(nvr)));
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_Random_System.");
      end if;
      return 548;
  end Dobldobl_PolySys_Random_System;

  function DoblDobl_PolySys_Make_Homogeneous
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    opt : constant natural32 := natural32(v_a(v_a'first));
    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last+1);

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_Make_Homogeneous ...");
    end if;
    Projective_Transformations.Projective_Transformation(lp.all);
    if opt = 0
     then res := Homogenization.Add_Random_Hyperplanes(lp.all,1,false);
     else res := Homogenization.Add_Standard_Hyperplanes(lp.all,1);
    end if;
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(res);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_Make_Homogeneous.");
      end if;
      return 892;
  end DoblDobl_PolySys_Make_Homogeneous;

  function DoblDobl_PolySys_Multi_Homogeneous
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    mhom : constant natural32 := natural32(v_a(v_a'first+1));
    opt : constant natural32 := natural32(v_a(v_a'first+2));
    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    md : constant integer32 := integer32(mhom);
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last+md);
    idz : Standard_Natural_Vectors.Vector(1..integer32(nvr));
    z : Partitions_of_Sets_of_Unknowns.Partition(1..mhom);

    use Multi_Projective_Transformations;

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_Multi_Homogeneous ...");
    end if;
    Assign(nvr,b,idz);
    z := Partitions_of_Sets_of_Unknowns_io.Make_Partition(nvr,mhom,idz);
    if opt = 0
     then res := Multi_Projective_Transformation(lp.all,mhom,z,false);
     else res := Multi_Projective_Transformation(lp.all,mhom,z,true);
    end if;
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(res);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_Multi_Homogeneous.");
      end if;
      return 905;
  end DoblDobl_PolySys_Multi_Homogeneous;

  function DoblDobl_PolySys_1Hom2Affine
             ( vrblvl : integer32 := 0 ) return integer32 is

    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    res : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last-1)
        := Affine_Transformations.Make_Affine(lp.all);

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_1Hom2Affine ...");
    end if;
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(res);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_1Hom2Affine.");
      end if;
      return 902;
  end DoblDobl_PolySys_1Hom2Affine;

  function DoblDobl_PolySys_mHom2Affine
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    mhom : constant natural32 := natural32(v_a(v_a'first));
    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    res : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
        := Affine_Transformations.Make_Affine(lp.all,mhom);

  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_mHom2Affine ...");
    end if;
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(res);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_mHom2Affine.");
      end if;
      return 907;
  end DoblDobl_PolySys_mHom2Affine;

  function DoblDobl_PolySys_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put_line("-> in DoblDobl_PolySys_Interface.DoblDobl_PolySys_Clear ...");
    end if;
    DoblDobl_PolySys_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in DoblDobl_PolySys_Interface.");
        put_line("DoblDobl_PolySys_Clear.");
      end if;
      return 207;
  end DoblDobl_PolySys_Clear;

  function DoblDobl_PolySys_Prompt_for_Target
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in dobldobl_polysys_interface.");
      put_line("DoblDobl_PolySys_Prompt_for_Target ...");
    end if;
    PHCpack_Operations_io.Read_DoblDobl_Target_System_without_Solutions;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in dobldobl_polysys_interface.");
        put_line("DoblDobl_PolySys_Prompt_for_Target.");
      end if;
      return 872;
  end DoblDobl_PolySys_Prompt_for_Target;

end DoblDobl_PolySys_Interface;
