with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Strings;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Random_Polynomials;
with Homogenization;
with Projective_Transformations;
with Affine_Transformations;
with Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;
with Multi_Projective_Transformations;
with Polynomial_Drops;
with Total_Degree_Start_Systems;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with PHCpack_Operations;
with PHCpack_Operations_io;
with Standard_PolySys_Container;

package body Standard_PolySys_Interface is

  function Standard_PolySys_Read
             ( vrblvl : integer32 := 0 ) return integer32 is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_polysys_interface.Standard_PolySys_Read ...");
    end if;
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    Standard_PolySys_Container.Initialize(lp.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Read.");
      end if;
      return 20;
  end Standard_PolySys_Read;

  function Standard_PolySys_Read_from_File
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
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Read_from_File ...");
    end if;
    declare
      file : file_type;
      p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    begin
      Open(file,in_file,sv);
      get(file,p);
      Standard_PolySys_Container.Clear;
      Standard_PolySys_Container.Initialize(p.all);
      exception 
        when NAME_ERROR =>
          put_line("File with name " & sv & " could not be found!");
          return 540;
        when USE_ERROR =>
          put_line("File with name " & sv & " could not be read!");
          return 540;
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Read_from_File.");
      end if;
      return 540;
  end Standard_PolySys_Read_from_File;

  function Standard_PolySys_Write
             ( vrblvl : in integer32 := 0 ) return integer32 is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    nvr : natural32;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_polysys_interface.Standard_PolySys_Write ...");
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
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Write.");
      end if;
      return 21;
  end Standard_PolySys_Write;

  function Standard_PolySys_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Get_Dimension ...");
    end if;
    Assign(integer32(Standard_PolySys_Container.Dimension),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Get_Dimension.");
      end if;
      return 22;
  end Standard_PolySys_Get_Dimension;

  function Standard_PolySys_Set_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant integer32 := integer32(v(v'first));

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Set_Dimension ...");
    end if;
    Standard_PolySys_Container.Initialize(n);
    Symbol_Table.Init(natural32(n));
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Set_Dimension.");
      end if;
      return 23;
  end Standard_PolySys_Set_Dimension;

  function Standard_PolySys_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    i : constant integer32 := integer32(v(v'first+1));

  begin
    if vrblvl > 0 then
      put_line("-> in standard_polysys_interface.Standard_PolySys_Size ...");
    end if;
    Assign(integer32(Standard_PolySys_Container.Number_of_Terms(i)),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Size");
      end if;
      return 24;
  end Standard_PolySys_Size;

  function Standard_PolySys_Degree
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    equ : constant integer32 := integer32(v_a(v_a'first));
    deg : constant integer32 := Standard_PolySys_Container.Degree(equ);

  begin
    if vrblvl > 0 then
      put_line("-> in standard_polysys_interface.Standard_PolySys_Degree ...");
    end if;
    Assign(deg,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Degree");
      end if;
      return 20;
  end Standard_PolySys_Degree;

  function Standard_PolySys_Total_Degree
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    td : natural32;

    use Total_Degree_Start_Systems;

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Total_Degree ...");
    end if;
    if lp = null then
      return 1;
    else
      td := Product(Total_Degree_Start_Systems.Degrees(lp.all));
      Assign(integer32(td),a);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Total_Degree");
      end if;
      return 28;
  end Standard_PolySys_Total_Degree;

  function Standard_PolySys_Get_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array(0..2)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    i : constant integer32 := integer32(v(1));
    j : constant natural32 := natural32(v(2));
    t : constant Standard_Complex_Polynomials.Term
      := Standard_PolySys_Container.Retrieve_Term(i,j);

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Get_Term ...");
    end if;
    Assign(t.cf,c);
    Assign(t.dg.all,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Get_Term");
      end if;
      return 25;
  end Standard_PolySys_Get_Term;

  function Standard_PolySys_Add_Term
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(v(0));
    i : constant integer32 := integer32(v(1));
    e : Standard_Natural_Vectors.Vector(1..n);
    t : Standard_Complex_Polynomials.Term;

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Add_Term ...");
    end if;
    Assign(c,t.cf);
    Assign(natural32(n),b,e);
    t.dg := new Standard_Natural_Vectors.Vector'(e);
    Standard_PolySys_Container.Add_Term(i,t);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Add_Term");
      end if;
      return 26;
  end Standard_PolySys_Add_Term;

  function Standard_PolySys_String_Save
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
    p : Standard_Complex_Polynomials.Poly;

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_String_Save.");
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
    p := Standard_Complex_Poly_Strings.Parse(n,s);
    Standard_PolySys_Container.Add_Poly(k,p);
    Standard_Complex_Polynomials.Clear(p);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_String_Save.");
      end if;
      return 76;
  end Standard_PolySys_String_Save;

  function Standard_PolySys_String_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    equ : constant integer32 := integer32(v_a(v_a'first));
    p : constant Standard_Complex_Polynomials.Poly
      := Standard_PolySys_Container.Retrieve_Poly(equ);
    sz : constant integer32
       := integer32(Standard_Complex_Poly_Strings.Size_Limit(p));

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_String_Size ...");
    end if;
    Assign(sz,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_String_Size.");
      end if;
      return 600;
  end Standard_PolySys_String_Size;

  function Standard_PolySys_String_Load 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    equ : constant integer32 := integer32(v_a(v_a'first));
    p : constant Standard_Complex_Polynomials.Poly
      := Standard_PolySys_Container.Retrieve_Poly(equ);
    s : constant string := Standard_Complex_Poly_Strings.Write(p);
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(s);
    slast : constant integer32 := integer32(s'last);

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_String_Load.");
    end if;
   -- put("Polynomial "); put(equ,1); put(" : "); put_line(s);
   -- put("#characters : "); put(s'last,1); new_line;
    Assign(slast,a);
    Assign(sv,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_String_Load.");
      end if;
      return 67;
  end Standard_PolySys_String_Load;

  function Standard_PolySys_Drop_by_Index
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant integer32 := integer32(v_a(v_a'first));
    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    dropped : constant Poly_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,ind);

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Drop_by_Index ...");
    end if;
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(dropped);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Drop_by_Index.");
      end if;
      return 306;
  end Standard_PolySys_Drop_by_Index; 

  function Standard_PolySys_Drop_by_Name
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc)
       := C_Integer_Array_to_String(natural32(nc),vb);   
    sb : Symbol_Table.Symbol;
    ind : natural32;
    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    dropped : Poly_Sys(lp'range);

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Drop_by_Name ...");
    end if;
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    dropped := Polynomial_Drops.Drop(lp.all,integer32(ind));
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(dropped);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Drop_by_Name.");
      end if;
      return 309;
  end Standard_PolySys_Drop_by_Name;

  function Standard_PolySys_Random_System
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    neq : constant integer32 := integer32(v_a(v_a'first+1));
    p : Standard_Complex_Poly_Systems.Poly_Sys(1..neq);
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(3));
    m : constant natural32 := natural32(v_b(v_b'first));
    d : constant natural32 := natural32(v_b(v_b'first+1));
    c : constant natural32 := natural32(v_b(v_b'first+2));

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Random_System ...");
    end if;
    for i in p'range loop
      if m = 0 then
        p(i) := Standard_Random_Polynomials.Random_Dense_Poly(nvr,d,c);
      else
        p(i) := Standard_Random_Polynomials.Random_Sparse_Poly(nvr,d,m,c);
      end if;
    end loop;
    Standard_PolySys_Container.Clear; 
    Standard_PolySys_Container.Initialize(p); 
   -- must initialize the symbol table with actual symbols for printing
    Symbol_Table.Init(Symbol_Table.Standard_Symbols(integer32(nvr)));
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Random_System.");
      end if;
      return 109;
  end Standard_PolySys_Random_System;

  function Standard_PolySys_Make_Homogeneous
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    opt : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    res : Standard_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last+1);

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Make_Homogeneous ...");
    end if;
    Projective_Transformations.Projective_Transformation(lp.all);
    if opt = 0
     then res := Homogenization.Add_Random_Hyperplanes(lp.all,1,false);
     else res := Homogenization.Add_Standard_Hyperplanes(lp.all,1);
    end if;
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(res);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Make_Homogeneous.");
      end if;
      return 891;
  end Standard_PolySys_Make_Homogeneous;

  function Standard_PolySys_Multi_Homogeneous
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    mhom : constant natural32 := natural32(v_a(v_a'first+1));
    opt : constant natural32 := natural32(v_a(v_a'first+2));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    md : constant integer32 := integer32(mhom);
    res : Standard_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last+md);
    idz : Standard_Natural_Vectors.Vector(1..integer32(nvr));
    z : Partitions_of_Sets_of_Unknowns.Partition(1..mhom);

    use Multi_Projective_Transformations;

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Multi_Homogeneous ...");
    end if;
    Assign(nvr,b,idz);
    z := Partitions_of_Sets_of_Unknowns_io.Make_Partition(nvr,mhom,idz);
    if opt = 0
     then res := Multi_Projective_Transformation(lp.all,mhom,z,false);
     else res := Multi_Projective_Transformation(lp.all,mhom,z,true);
    end if;
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(res);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Multi_Homogeneous.");
      end if;
      return 904;
  end Standard_PolySys_Multi_Homogeneous;

  function Standard_PolySys_1Hom2Affine
             ( vrblvl : integer32 := 0 ) return integer32 is

    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    res : constant Standard_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last-1)
        := Affine_Transformations.Make_Affine(lp.all);

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_1Hom2Affine ...");
    end if;
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(res);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_1Hom2Affine.");
      end if;
      return 906;
  end Standard_PolySys_1Hom2Affine;

  function Standard_PolySys_mHom2Affine
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    mhom : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    res : constant Standard_Complex_Poly_Systems.Poly_Sys
        := Affine_Transformations.Make_Affine(lp.all,mhom);

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_mHom2Affine ...");
    end if;
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(res);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_mHom2Affine.");
      end if;
      return 907;
  end Standard_PolySys_mHom2Affine;

  function Standard_PolySys_Make_Function
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Make_Function ...");
    end if;
    if Standard_PolySys_Container.Retrieve = null then
      return 147;
    else
      Standard_PolySys_Container.Create_Evaluator;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Make_Function.");
      end if;
      return 147;
  end Standard_PolySys_Make_Function;

  function Standard_PolySys_Jacobian_Function
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Jacobian_Function ...");
    end if;
    if Standard_PolySys_Container.Retrieve = null then
      return 11;
    else
      Standard_PolySys_Container.Create_Jacobian_Evaluator;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Jacobian_Function.");
      end if;
      return 148;
  end Standard_PolySys_Jacobian_Function;

  function Standard_PolySys_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_polysys_interface.Standard_PolySys_Clear ...");
    end if;
    Standard_PolySys_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Clear.");
      end if;
      return 27;
  end Standard_PolySys_Clear;

  function Standard_PolySys_Prompt_for_Target
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Prompt_for_Target ...");
    end if;
    PHCpack_Operations_io.Read_Target_System_without_Solutions;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Prompt_for_Target.");
      end if;
      return 150;
  end Standard_PolySys_Prompt_for_Target;

  function Standard_PolySys_Prompt_for_Start
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Prompt_for_Start ...");
    end if;
    PHCpack_Operations_io.Read_Start_System_without_Solutions;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Prompt_for_Start.");
      end if;
      return 151;
  end Standard_PolySys_Prompt_for_Start;

  function Standard_PolySys_Read_Target_on_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant integer := integer(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    s : constant String(1..n) := C_Integer_Array_to_String(natural32(n),v_b);

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Read_Target_on_File ...");
    end if;
   -- put_line("opening the file " & s & " for the target system ...");
    PHCpack_Operations_io.Read_Target_System_without_Solutions(s);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Read_Target_on_File.");
      end if;
      return 161;
  end Standard_PolySys_Read_Target_on_File;

  function Standard_PolySys_Read_Start_on_File
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant integer := integer(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    s : constant String(1..n) := C_Integer_Array_to_String(natural32(n),v_b);

  begin
    if vrblvl > 0 then
      put("-> in standard_polysys_interface.");
      put_line("Standard_PolySys_Read_Start_on_File ...");
    end if;
   -- put_line("opening the file " & s & " for the start system ...");
    PHCpack_Operations_io.Read_Start_System_without_Solutions(s);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in standard_polysys_interface.");
        put_line("Standard_PolySys_Read_Start_on_File.");
      end if;
      return 162;
  end Standard_PolySys_Read_Start_on_File;

end Standard_PolySys_Interface;
