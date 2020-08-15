with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
--with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
--with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
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
with PHCpack_Operations_io;
with Standard_PolySys_Container;
with Standard_Solutions_Container;

package body Linear_Products_Interface is

  function Linear_Products_Structure_Make
             ( vrblvl : integer32 := 0 ) return integer32 is

    p : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_Structure_Make ...");
    end if;
    Random_Product_Start_Systems.Build_Set_Structure(p.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_Structure_Make.");
      end if;
      return 110;
  end Linear_Products_Structure_Make;

  function Linear_Products_Structure_Write
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_Structure_Make ...");
    end if;
    Set_Structure_io.put;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_Structure_Write.");
      end if;
      return 111;
  end Linear_Products_Structure_Write;

  function Linear_Products_Structure_Bound
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

   -- r : constant natural := Set_Structure.B;
    r : constant integer32 
      := Degree_Sets_Tables.Permanent(Degree_Sets_Tables.Create);

  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_Structure_Bound ...");
    end if;
    Assign(r,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_Structure_Bound.");
      end if;
      return 112;
  end Linear_Products_Structure_Bound;

  function Linear_Products_System_Make
             ( vrblvl : integer32 := 0 ) return integer32 is

    n : constant natural32 := Standard_PolySys_Container.Dimension;
    q : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_System_Make ...");
    end if;
    Standard_Linear_Product_System.Init(n);
    Random_Product_Start_Systems.Build_Random_Product_System(n);
    q := new Poly_Sys'(Standard_Linear_Product_System.Polynomial_System);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_System_Make.");
      end if;
      return 113;
  end Linear_Products_System_Make;

  function Linear_Products_System_Read
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
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_System_Read ...");
    end if;
   -- put_line("opening the file " & s & " for the start system ...");
    PHCpack_Operations_io.Read_Linear_Product_Start_System(s,fail);
    if fail
     then return 163;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_System_Read.");
      end if;
      return 163;
  end Linear_Products_System_Read;

  function Linear_Products_System_Solve
             ( vrblvl : integer32 := 0 ) return integer32 is

    sols : Solution_List;
    nb : natural32;

  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_System_Solve ...");
    end if;
    Standard_Linear_Product_System.Solve(sols,nb);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_System_Solve.");
      end if;
      return 114;
  end Linear_Products_System_Solve;

  function Linear_Products_Structure_String_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is 

    s : constant string := Set_Structure_Strings.to_string;
    slast : constant integer32 := integer32(s'last);
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(s);

  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_Structure_String_Get ...");
    end if;
    Assign(slast,a);
    Assign(sv,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_Structure_String_Get.");
      end if;
      return 116;
  end Linear_Products_Structure_String_Get;

  function Linear_Products_Structure_String_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant integer := integer(v_a(v_a'first));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);

  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_Structure_String_Set ...");
    end if;
    Set_Structure_Strings.Parse(s);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_Structure_String_Set.");
      end if;
      return 117;
  end Linear_Products_Structure_String_Set;

  function Linear_Products_Structure_Check
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    val : constant integer32 := integer32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    verbose : constant boolean := (val = 1);
    sup : boolean := false;

  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_Structure_Check ...");
    end if;
    if lp /= null
     then sup := Supporting_Set_Structure.Is_Supporting(lp.all,verbose);
    end if;
    if sup
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_Structure_Check.");
      end if;
      return 118;
  end Linear_Products_Structure_Check;

  function Linear_Products_Partition_Make
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    p : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    n : constant natural32 := natural32(p'last);
    z : Partition(1..n);
    bzn : natural64;
    nbz : natural32;

  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_Partition_Make ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_Partition_Make.");
      end if;
      return 530;
  end Linear_Products_Partition_Make;

  function Linear_Products_Partition_Bound
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_Partition_Bound ...");
    end if;
   -- put_line("the string for parsing : " & s);
   -- put("The partition parsed : "); put(z); new_line;
    mbz := m_Homogeneous_Bezout_Numbers.Bezout_Number(p.all,z);
   -- put("mbz = "); put(mbz,1); new_line;
    Assign(integer32(mbz),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_Partition_Bound.");
      end if;
      return 531;
  end Linear_Products_Partition_Bound;

  function Linear_Products_Partition_System
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_Partition_System ...");
    end if;
    m_Homogeneous_Start_System(p.all,z);
    q := new Poly_Sys'(Standard_Linear_Product_System.Polynomial_System);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_Partition_System.");
      end if;
      return 532;
  end Linear_Products_Partition_System;

  function Linear_Products_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in linear_products_interface.");
      put_line("Linear_Products_Clear ...");
    end if;
    Set_Structure.Clear;
    Standard_Linear_Product_System.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in linear_products_interface.");
        put_line("Linear_Products_Clear.");
      end if;
      return 115;
  end Linear_Products_Clear;

end Linear_Products_Interface;
