with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Random_Product_Start_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Partitions_of_Sets_of_Unknowns;    use Partitions_of_Sets_of_Unknowns;
--with Partitions_of_Sets_of_Unknowns_io; use Partitions_of_Sets_of_Unknowns_io;
with Partitions_of_Sets_Strings;
with m_Homogeneous_Bezout_Numbers;
with m_Homogeneous_Start_Systems;       use m_Homogeneous_Start_Systems;
with Set_Structure,Set_Structure_io;
with Set_Structure_Strings;
with Supporting_Set_Structure;
with Degree_Sets_Tables;
with Standard_Linear_Product_System;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with Standard_Solutions_Container;

function use_roco ( job : integer32;
                    a : C_intarrs.Pointer;
                    b : C_intarrs.Pointer;
                    c : C_dblarrs.Pointer ) return integer32 is

  function Job0 return integer32 is  -- construct supporting set structure
  
    p : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;

  begin
    Random_Product_Start_Systems.Build_Set_Structure(p.all);
    return 0;
  end Job0;

  function Job1 return integer32 is -- write supporting set structure
  begin
    Set_Structure_io.put;
    return 0;
  end Job1;

  function Job2 return integer32 is -- compute the Bezout bound

   -- r : constant natural := Set_Structure.B;
    r : constant integer32 
      := Degree_Sets_Tables.Permanent(Degree_Sets_Tables.Create);

  begin
    Assign(r,a);
    return 0;
  end Job2;

  function Job3 return integer32 is -- make random linear-product system

    n : constant natural32 := Standard_PolySys_Container.Dimension;
    q : Link_to_Poly_Sys;

  begin
    Standard_Linear_Product_System.Init(n);
    Random_Product_Start_Systems.Build_Random_Product_System(n);
    q := new Poly_Sys'(Standard_Linear_Product_System.Polynomial_System);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(q.all);
    return 0;
  end Job3;

  function Job4 return integer32 is -- solve random linear-product system

    sols : Solution_List;
    nb : natural32;

  begin
    Standard_Linear_Product_System.Solve(sols,nb);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(sols);
    return 0;
  end Job4;

  function Job5 return integer32 is -- clears the set structure
  begin
    Set_Structure.Clear;
    Standard_Linear_Product_System.Clear;
    return 0;
  end Job5;

  function Job6 return integer32 is -- return set structure string

    s : constant string := Set_Structure_Strings.to_string;
    slast : constant integer32 := integer32(s'last);
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(s);

  begin
    Assign(slast,a);
    Assign(sv,b);
    return 0;
  end Job6;

  function Job7 return integer32 is -- parse string into set structure

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant integer := integer(v_a(v_a'first));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);

  begin
    Set_Structure_Strings.Parse(s);
    return 0;
  end Job7;

  function Job8 return integer32 is -- check if set structure supports

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    val : constant integer32 := integer32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    verbose : constant boolean := (val = 1);
    sup : boolean := false;

  begin
    if lp /= null
     then sup := Supporting_Set_Structure.Is_Supporting(lp.all,verbose);
    end if;
    if sup
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  end Job8;

  function Job10 return integer32 is -- create partition for container system

    p : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    n : constant natural32 := natural32(p'last);
    z : Partition(1..n);
    bzn : natural64;
    nbz : natural32;

  begin
    m_Homogeneous_Bezout_Numbers.PB(p.all,bzn,nbz,z);
    m_Homogeneous_Bezout_Numbers.Patch(p.all,z,nbz,bzn);
    declare
      sz : constant string := Partitions_of_Sets_Strings.to_String(z(1..nbz));
      sv : constant Standard_Integer_Vectors.Vector
         := String_to_Integer_Vector(sz);
      sc : Standard_Integer_Vectors.Vector(1..2);
    begin
      sc(1) := integer32(bzn);
      sc(2) := sv'last;
      Assign(sc,a);
      Assign(sv,b);
    end;
    return 0;
  end Job10;

  function Job11 return integer32 is -- evaluate partition on system

    p : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    dim : constant natural32 := natural32(p'last);
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant integer := integer(v_a(v_a'first));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    z : constant Partition := Partitions_of_Sets_Strings.Parse(s,dim);
    mbz : natural64;

  begin
   -- put_line("the string for parsing : " & s);
   -- put("The partition parsed : "); put(z); new_line;
    mbz := m_Homogeneous_Bezout_Numbers.Bezout_Number(p.all,z);
   -- put("mbz = "); put(mbz,1); new_line;
    Assign(integer32(mbz),a);
    return 0;
  end Job11;

  function Job12 return integer32 is -- linear-product system

    p : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    dim : constant natural32 := natural32(p'last);
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant integer := integer(v_a(v_a'first));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    z : constant Partition := Partitions_of_Sets_Strings.Parse(s,dim);
    q : Link_to_Poly_Sys;

  begin
    m_Homogeneous_Start_System(p.all,z);
    q := new Poly_Sys'(Standard_Linear_Product_System.Polynomial_System);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(q.all);
    return 0;
  end Job12;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0;   -- supporting set structure
      when 1 => return Job1;   -- write supporting set structure
      when 2 => return Job2;   -- compute the Bezout bound
      when 3 => return Job3;   -- make random linear-product system
      when 4 => return Job4;   -- solve random linear-product system
      when 5 => return Job5;   -- clears the set structure
      when 6 => return Job6;   -- returns set structure string
      when 7 => return Job7;   -- parse string into set structure
      when 8 => return Job8;   -- does set structure support polynomials?
      when 10 => return Job10; -- creates partition for system in container
      when 11 => return Job11; -- evaluates partition on container system
      when 12 => return Job12; -- linear-product system based on partition 
      when others => put_line("invalid operation"); return 1;
    end case;
  exception
    when others => put("Exception raised in use_roco handling job ");
                   put(job,1); put_line(".  Will try to ignore...");
                   return 1;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_roco;
