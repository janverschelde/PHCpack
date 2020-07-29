with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
-- with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
--   with Standard_Integer_Numbers_io;
--    use Standard_Integer_Numbers_io;
--   with Standard_Complex_Polynomials_io;
--    use Standard_Complex_Polynomials_io;
with Standard_Random_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Symbol_Table;
with DoblDobl_Random_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Random_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Multprec_Complex_Laurentials;
with Polynomial_Drops;
with Projective_Transformations;
with Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;
with Multi_Projective_Transformations;
with Affine_Transformations;
with Homogenization;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Standard_LaurSys_Container;
with DoblDobl_LaurSys_Container;
with QuadDobl_LaurSys_Container;
with Multprec_LaurSys_Container;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Interface;
with Standard_LaurSys_Interface;
with DoblDobl_PolySys_Interface;
with DoblDobl_LaurSys_Interface;
with QuadDobl_PolySys_Interface;
with QuadDobl_LaurSys_Interface;
with Multprec_PolySys_Interface;
with Multprec_LaurSys_Interface;

function use_syscon ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function Job135 return integer32 is -- returns a term of a Laurential

    v : constant C_Integer_Array(0..2)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    i : constant integer32 := integer32(v(1));
    j : constant natural32 := natural32(v(2));
    t : constant Multprec_Complex_Laurentials.Term 
      := Multprec_LaurSys_Container.Retrieve_Term(i,j);

  begin
   -- Assign(t.cf,c);
    Assign(t.dg.all,b);
    return 0;
  end Job135;

  function Job136 return integer32 is -- add a term to a Laurential

    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(v(0));
    i : constant integer32 := integer32(v(1));
    e : Standard_Integer_Vectors.Vector(1..n);
    t : Multprec_Complex_Laurentials.Term;

  begin
   -- Assign(c,t.cf);
    Assign(natural32(n),b,e);
    t.dg := new Standard_Integer_Vectors.Vector'(e);
    Multprec_LaurSys_Container.Add_Term(i,t);
    return 0;
  end Job136;

  function Job10 return integer32 is -- creates an evaluator

    use Standard_Complex_Poly_Systems;

  begin
    if Standard_PolySys_Container.Retrieve = null then
      return 10;
    else
      Standard_PolySys_Container.Create_Evaluator;
      return 0;
    end if;
  end Job10;

  function Job11 return integer32 is -- creates a Jacobian matrix evaluator

    use Standard_Complex_Poly_Systems;

  begin
    if Standard_PolySys_Container.Retrieve = null then
      return 11;
    else
      Standard_PolySys_Container.Create_Jacobian_Evaluator;
      return 0;
    end if;
  end Job11;

  function Job71 return integer32 is -- puts random system in the container

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
    when others => return 71;
  end Job71;

  function Job78 return integer32 is -- random dobldobl system

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
    when others => return 78;
  end Job78;

  function Job79 return integer32 is -- random quaddobl system

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    neq : constant integer32 := integer32(v_a(v_a'first+1));
    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..neq);
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(3));
    m : constant natural32 := natural32(v_b(v_b'first));
    d : constant natural32 := natural32(v_b(v_b'first+1));
    c : constant natural32 := natural32(v_b(v_b'first+2));

  begin
    for i in p'range loop
      if m = 0 then
        p(i) := QuadDobl_Random_Polynomials.Random_Dense_Poly(nvr,d,c);
      else
        p(i) := QuadDobl_Random_Polynomials.Random_Sparse_Poly(nvr,d,m,c);
      end if;
    end loop;
    QuadDobl_PolySys_Container.Clear; 
    QuadDobl_PolySys_Container.Initialize(p); 
   -- must initialize the symbol table with actual symbols for printing
    Symbol_Table.Init(Symbol_Table.Standard_Symbols(integer32(nvr)));
    return 0;
  exception
    when others => return 79;
  end Job79;

-- The jobs to drop a coordinate of a system come in two flavors:
-- (1) by index: given the index of the variable in a[0];
-- (2) by name: given the number of characters of the symbol in a[0]
--     and the characters for the symbol name in b.

  function Job12 return integer32 is -- drop variable by index

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant integer32 := integer32(v_a(v_a'first));
    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    dropped : constant Poly_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,ind);

  begin
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(dropped);
    return 0;
  end Job12;

  function Job13 return integer32 is -- drop variable by index

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant integer32 := integer32(v_a(v_a'first));
    use DoblDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    dropped : constant Poly_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,ind);

  begin
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(dropped);
    return 0;
  end Job13;

  function Job14 return integer32 is -- drop variable by index

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant integer32 := integer32(v_a(v_a'first));
    use QuadDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    dropped : constant Poly_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,ind);

  begin
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(dropped);
    return 0;
  end Job14;

  function Job15 return integer32 is -- standard drop variable by name

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
  end Job15;

  function Job16 return integer32 is -- dobldobl drop variable by name

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
  end Job16;

  function Job17 return integer32 is -- quaddobl drop variable by name

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);
    sb : Symbol_Table.Symbol;
    ind : natural32;
    use Quaddobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    dropped : Poly_Sys(lp'range);

  begin
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    dropped := Polynomial_Drops.Drop(lp.all,integer32(ind));
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(dropped);
    return 0;
  end Job17;

  function Job22 return integer32 is -- drop standard Laurent variable by idx

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant integer32 := integer32(v_a(v_a'first));
    use Standard_Complex_Laur_Systems;
    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    dropped : constant Laur_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,ind);

  begin
    Standard_LaurSys_Container.Clear;
    Standard_LaurSys_Container.Initialize(dropped);
    return 0;
  end Job22;

  function Job23 return integer32 is -- drop dobldobl Laurent variable by idx

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant integer32 := integer32(v_a(v_a'first));
    use DoblDobl_Complex_Laur_Systems;
    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    dropped : constant Laur_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,ind);

  begin
    DoblDobl_LaurSys_Container.Clear;
    DoblDobl_LaurSys_Container.Initialize(dropped);
    return 0;
  end Job23;

  function Job24 return integer32 is -- drop quaddobl Laurent variable by idx

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    ind : constant integer32 := integer32(v_a(v_a'first));
    use QuadDobl_Complex_Laur_Systems;
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    dropped : constant Laur_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,ind);

  begin
    QuadDobl_LaurSys_Container.Clear;
    QuadDobl_LaurSys_Container.Initialize(dropped);
    return 0;
  end Job24;

  function Job25 return integer32 is -- standard Laurent drop variable by name

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc)
       := C_Integer_Array_to_String(natural32(nc),vb);   
    sb : Symbol_Table.Symbol;
    ind : natural32;
    use Standard_Complex_Laur_Systems;
    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    dropped : Laur_Sys(lp'range);

  begin
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    dropped := Polynomial_Drops.Drop(lp.all,integer32(ind));
    Standard_LaurSys_Container.Clear;
    Standard_LaurSys_Container.Initialize(dropped);
    return 0;
  end Job25;

  function Job26 return integer32 is -- dobldobl Laurent drop variable by name

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant natural := natural(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),vb);
    sb : Symbol_Table.Symbol;
    ind : natural32;
    use DoblDobl_Complex_Laur_Systems;
    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    dropped : Laur_Sys(lp'range);

  begin
    for i in 1..nc loop
      sb(i) := sv(i);
    end loop;
    for i in nc+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    dropped := Polynomial_Drops.Drop(lp.all,integer32(ind));
    DoblDobl_LaurSys_Container.Clear;
    DoblDobl_LaurSys_Container.Initialize(dropped);
    return 0;
  end Job26;

  function Job27 return integer32 is -- quaddobl Laurent drop variable by name

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
  end Job27;

  function Job891 return integer32 is -- 1-homogeneous standard system

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    opt : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    res : Standard_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last+1);

  begin
    Projective_Transformations.Projective_Transformation(lp.all);
    if opt = 0
     then res := Homogenization.Add_Random_Hyperplanes(lp.all,1,false);
     else res := Homogenization.Add_Standard_Hyperplanes(lp.all,1);
    end if;
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(res);
    return 0;
  end Job891;

  function Job892 return integer32 is -- 1-homogeneous dobldobl system

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    opt : constant natural32 := natural32(v_a(v_a'first));
    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last+1);

  begin
    Projective_Transformations.Projective_Transformation(lp.all);
    if opt = 0
     then res := Homogenization.Add_Random_Hyperplanes(lp.all,1,false);
     else res := Homogenization.Add_Standard_Hyperplanes(lp.all,1);
    end if;
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(res);
    return 0;
  end Job892;

  function Job893 return integer32 is -- 1-homogeneous quaddobl system

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    opt : constant natural32 := natural32(v_a(v_a'first));
    lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last+1);

  begin
    Projective_Transformations.Projective_Transformation(lp.all);
    if opt = 0
     then res := Homogenization.Add_Random_Hyperplanes(lp.all,1,false);
     else res := Homogenization.Add_Standard_Hyperplanes(lp.all,1);
    end if;
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(res);
    return 0;
  end Job893;

  function Job904 return integer32 is -- m-homogeneous standard system

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
    Assign(nvr,b,idz);
    z := Partitions_of_Sets_of_Unknowns_io.Make_Partition(nvr,mhom,idz);
    if opt = 0
     then res := Multi_Projective_Transformation(lp.all,mhom,z,false);
     else res := Multi_Projective_Transformation(lp.all,mhom,z,true);
    end if;
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(res);
    return 0;
  end Job904;

  function Job905 return integer32 is -- m-homogeneous dobldobl system

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
    Assign(nvr,b,idz);
    z := Partitions_of_Sets_of_Unknowns_io.Make_Partition(nvr,mhom,idz);
    if opt = 0
     then res := Multi_Projective_Transformation(lp.all,mhom,z,false);
     else res := Multi_Projective_Transformation(lp.all,mhom,z,true);
    end if;
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(res);
    return 0;
  end Job905;

  function Job906 return integer32 is -- m-homogeneous quaddobl system

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    mhom : constant natural32 := natural32(v_a(v_a'first+1));
    opt : constant natural32 := natural32(v_a(v_a'first+2));
    lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    md : constant integer32 := integer32(mhom);
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last+md);
    idz : Standard_Natural_Vectors.Vector(1..integer32(nvr));
    z : Partitions_of_Sets_of_Unknowns.Partition(1..mhom);

    use Multi_Projective_Transformations;

  begin
    Assign(nvr,b,idz);
    z := Partitions_of_Sets_of_Unknowns_io.Make_Partition(nvr,mhom,idz);
    if opt = 0
     then res := Multi_Projective_Transformation(lp.all,mhom,z,false);
     else res := Multi_Projective_Transformation(lp.all,mhom,z,true);
    end if;
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(res);
    return 0;
  end Job906;

  function Job897 return integer32 is -- add symbol passed as string

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nc : constant integer := integer(v_a(v_a'first));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    smb : constant String(1..nc)
        := C_Integer_Array_to_String(natural32(nc),v_b);

  begin
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add_String(smb);
    return 0;
  end Job897;

  function Job901 return integer32 is -- double affine transformation

    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    res : constant Standard_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last-1)
        := Affine_Transformations.Make_Affine(lp.all);

  begin
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(res);
    return 0;
  end Job901;

  function Job902 return integer32 is -- double double affine transformation

    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    res : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last-1)
        := Affine_Transformations.Make_Affine(lp.all);

  begin
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(res);
    return 0;
  end Job902;

  function Job903 return integer32 is -- quad double affine transformation

    lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    res : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last-1)
        := Affine_Transformations.Make_Affine(lp.all);

  begin
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(res);
    return 0;
  end Job903;

  function Job907 return integer32 is -- double m-hom to affine

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    mhom : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    res : constant Standard_Complex_Poly_Systems.Poly_Sys
        := Affine_Transformations.Make_Affine(lp.all,mhom);

  begin
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(res);
    return 0;
  end Job907;

  function Job908 return integer32 is -- double double m-hom to affine

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    mhom : constant natural32 := natural32(v_a(v_a'first));
    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    res : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
        := Affine_Transformations.Make_Affine(lp.all,mhom);

  begin
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(res);
    return 0;
  end Job908;

  function Job909 return integer32 is -- quad double m-hom to affine

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    mhom : constant natural32 := natural32(v_a(v_a'first));
    lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    res : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
        := Affine_Transformations.Make_Affine(lp.all,mhom);

  begin
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(res);
    return 0;
  end Job909;

  function Handle_Jobs return integer32 is

    use Standard_PolySys_Interface;
    use Standard_LaurSys_Interface;
    use DoblDobl_PolySys_Interface;
    use DoblDobl_LaurSys_Interface;
    use QuadDobl_PolySys_Interface;
    use QuadDobl_LaurSys_Interface;
    use MultPrec_PolySys_Interface;
    use MultPrec_LaurSys_Interface;

  begin
    case job is
      when 0 => return Standard_PolySys_Read(vrblvl);
      when 1 => return Standard_PolySys_Write(vrblvl);
      when 2 => return Standard_PolySys_Get_Dimension(a,vrblvl);
      when 3 => return Standard_PolySys_Set_Dimension(a,vrblvl);
      when 4 => return Standard_PolySys_Size(a,vrblvl);
      when 5 => return Standard_PolySys_Get_Term(a,b,c,vrblvl);
      when 6 => return Standard_PolySys_Add_Term(a,b,c,vrblvl);
      when 7 => return Standard_PolySys_Clear(vrblvl);
      when 8 => return Standard_PolySys_Total_Degree(a,vrblvl);
      when 9 => return Standard_PolySys_Clear_Symbols(vrblvl);
      when 10 => return Job10; -- creates a system evaluator
      when 11 => return Job11; -- creates a Jacobian matrix evaluator
     -- dropping variables from polynomials
      when 12 => return Job12; -- standard drop variable by index
      when 13 => return Job13; -- dobldobl drop variable by index
      when 14 => return Job14; -- quaddobl drop variable by index
      when 15 => return Job15; -- standard drop variable by name
      when 16 => return Job16; -- dobldobl drop variable by name
      when 17 => return Job17; -- quaddobl drop variable by name
     -- degrees of polynomials :
      when 20 => return Standard_PolySys_Degree(a,b,vrblvl);
     -- dropping variables from Laurent polynomials
      when 22 => return Job22; -- standard Laurent drop variable by index
      when 23 => return Job23; -- dobldobl Laurent drop variable by index
      when 24 => return Job24; -- quaddobl Laurent drop variable by index
      when 25 => return Job25; -- standard Laurent drop variable by name
      when 26 => return Job26; -- dobldobl Laurent drop variable by name
      when 27 => return Job27; -- quaddobl Laurent drop variable by name
     -- jobs for standard double complex Laurent polynomials :
      when 100 => return Standard_LaurSys_Read(vrblvl);
      when 101 => return Standard_LaurSys_Write(vrblvl);
      when 102 => return Standard_LaurSys_Get_Dimension(a,vrblvl);
      when 103 => return Standard_LaurSys_Set_Dimension(a,vrblvl);
      when 104 => return Standard_LaurSys_Size(a,vrblvl);
      when 105 => return Standard_LaurSys_Get_Term(a,b,c,vrblvl);
      when 106 => return Standard_LaurSys_Add_Term(a,b,c,vrblvl);
      when 107 => return Standard_LaurSys_Clear(vrblvl);
     -- jobs for double double complex Laurent polynomials :
      when 110 => return DoblDobl_LaurSys_Read(vrblvl);
      when 111 => return DoblDobl_LaurSys_Write(vrblvl);
      when 112 => return DoblDobl_LaurSys_Get_Dimension(a,vrblvl);
      when 113 => return DoblDobl_LaurSys_Set_Dimension(a,vrblvl);
      when 114 => return DoblDobl_LaurSys_Size(a,vrblvl);
      when 115 => return DoblDobl_LaurSys_Get_Term(a,b,c,vrblvl);
      when 116 => return DoblDobl_LaurSys_Add_Term(a,b,c,vrblvl);
      when 117 => return DoblDobl_LaurSys_Clear(vrblvl);
      when 118 => return DoblDobl_LaurSys_String_Save(a,b,vrblvl);
     -- jobs for quad double complex Laurent polynomials :
      when 120 => return QuadDobl_LaurSys_Read(vrblvl);
      when 121 => return QuadDobl_LaurSys_Write(vrblvl);
      when 122 => return QuadDobl_LaurSys_Get_Dimension(a,vrblvl);
      when 123 => return QuadDobl_LaurSys_Set_Dimension(a,vrblvl);
      when 124 => return QuadDobl_LaurSys_Size(a,vrblvl);
      when 125 => return QuadDobl_LaurSys_Get_Term(a,b,c,vrblvl);
      when 126 => return QuadDobl_LaurSys_Add_Term(a,b,c,vrblvl);
      when 127 => return QuadDobl_LaurSys_Clear(vrblvl);
      when 128 => return QuadDobl_LaurSys_String_Save(a,b,vrblvl);
     -- jobs for multiprecision complex Laurent polynomials :
      when 130 => return Multprec_LaurSys_Read(vrblvl);
      when 131 => return Multprec_LaurSys_Write(vrblvl);
      when 132 => return Multprec_LaurSys_Get_Dimension(a,vrblvl);
      when 133 => return Multprec_LaurSys_Set_Dimension(a,vrblvl);
      when 134 => return Multprec_LaurSys_Size(a,vrblvl);
      when 135 => return Job135; -- return a term of a polynomial
      when 136 => return Job136; -- add a term to a polynomial
      when 137 => return Multprec_LaurSys_Clear(vrblvl);
      when 138 => return Multprec_LaurSys_String_Save(a,b,vrblvl);
      when 139 => return Multprec_LaurSys_String_Load(a,b,vrblvl);
     -- jobs for double double complex polynomials
      when 200 => return DoblDobl_PolySys_Read(vrblvl);
      when 201 => return DoblDobl_PolySys_Write(vrblvl);
      when 202 => return DoblDobl_PolySys_Get_Dimension(a,vrblvl);
      when 203 => return DoblDobl_PolySys_Set_Dimension(a,vrblvl);
      when 204 => return DoblDobl_PolySys_Size(a,vrblvl);
      when 205 => return DoblDobl_PolySys_Get_Term(a,b,c,vrblvl);
      when 206 => return DoblDobl_PolySys_Add_Term(a,b,c,vrblvl);
      when 207 => return DoblDobl_PolySys_Clear(vrblvl);
      when 208 => return DoblDobl_PolySys_String_Save(a,b,vrblvl);
      when 209 => return DoblDobl_PolySys_Degree(a,b,vrblvl);
     -- jobs for quad double complex polynomials
      when 210 => return QuadDobl_PolySys_Read(vrblvl);
      when 211 => return QuadDobl_PolySys_Write(vrblvl);
      when 212 => return QuadDobl_PolySys_Get_Dimension(a,vrblvl);
      when 213 => return QuadDobl_PolySys_Set_Dimension(a,vrblvl);
      when 214 => return QuadDobl_PolySys_Size(a,vrblvl);
      when 215 => return QuadDobl_PolySys_Get_Term(a,b,c,vrblvl);
      when 216 => return QuadDobl_PolySys_Add_Term(a,b,c,vrblvl);
      when 217 => return QuadDobl_PolySys_Clear(vrblvl);
      when 218 => return QuadDobl_PolySys_String_Save(a,b,vrblvl);
      when 219 => return QuadDobl_PolySys_Degree(a,b,vrblvl);
     -- jobs for multiprecision complex polynomials
      when 220 => return Multprec_PolySys_Read(vrblvl);
      when 221 => return Multprec_PolySys_Write(vrblvl);
      when 222 => return Multprec_PolySys_Get_Dimension(a,vrblvl);
      when 223 => return Multprec_PolySys_Set_Dimension(a,vrblvl);
      when 224 => return Multprec_PolySys_Size(a,vrblvl);
      when 227 => return Multprec_PolySys_Clear(vrblvl);
      when 228 => return Multprec_PolySys_String_Save(a,b,vrblvl);
      when 229 => return Multprec_PolySys_Degree(a,b,vrblvl);
     -- jobs for interchanging polynomial as strings :
      when 67 => return Standard_PolySys_String_Load(a,b,vrblvl);
      when 68 => return DoblDobl_PolySys_String_Load(a,b,vrblvl);
      when 69 => return QuadDobl_PolySys_String_Load(a,b,vrblvl);
      when 70 => return Multprec_PolySys_String_Load(a,b,vrblvl);
      when 71 => return Job71; -- store random system in container
      when 72 => return DoblDobl_LaurSys_String_Load(a,b,vrblvl);
      when 73 => return QuadDobl_LaurSys_String_Load(a,b,vrblvl);
      when 74 => return Standard_LaurSys_String_Save(a,b,vrblvl);
      when 76 => return Standard_PolySys_String_Save(a,b,vrblvl);
      when 77 => return Standard_LaurSys_String_Load(a,b,vrblvl);
     -- random systems in double double and quad double precision
      when 78 => return Job78; -- store random system in dobldobl container 
      when 79 => return Job79; -- store random system in quaddobl container
     -- jobs to return the size limit of the string representations
      when 80 => return Standard_PolySys_String_Size(a,b,vrblvl);
      when 81 => return DoblDobl_PolySys_String_Size(a,b,vrblvl);
      when 82 => return QuadDobl_PolySys_String_Size(a,b,vrblvl);
      when 83 => return Multprec_PolySys_String_Size(a,b,vrblvl);
      when 84 => return Standard_LaurSys_String_Size(a,b,vrblvl);
      when 85 => return DoblDobl_LaurSys_String_Size(a,b,vrblvl);
      when 86 => return QuadDobl_LaurSys_String_Size(a,b,vrblvl);
      when 87 => return Multprec_LaurSys_String_Size(a,b,vrblvl);
     -- reading systems into the containers :
      when 540 => return Standard_PolySys_Read_from_File(a,b,vrblvl); 
      when 541 => return DoblDobl_PolySys_Read_from_File(a,b,vrblvl); 
      when 542 => return QuadDobl_PolySys_Read_from_File(a,b,vrblvl); 
      when 543 => return Multprec_PolySys_Read_from_File(a,b,vrblvl); 
     -- projective transformations :
      when 891 => return Job891; -- 1-homogeneous standard system
      when 892 => return Job892; -- 1-homogeneous dobldobl system
      when 893 => return Job893; -- 1-homogeneous quaddobl system
      when 904 => return Job904; -- m-homogeneous standard system
      when 905 => return Job905; -- m-homogeneous dobldobl system
      when 906 => return Job906; -- m-homogeneous quaddobl system
     -- add symbol passed as string to the table
      when 897 => return Job897; -- add symbol to the table
     -- affine transformations :
      when 901 => return Job901; -- double affine transformation
      when 902 => return Job902; -- double double affine transformation
      when 903 => return Job903; -- quad double affine transformation
      when 907 => return Job907; -- double m-hom to affine
      when 908 => return Job908; -- double double m-hom to affine
      when 909 => return Job909; -- quad double m-hom to affine
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_syscon;
