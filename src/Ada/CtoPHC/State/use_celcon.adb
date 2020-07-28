--with Standard_Integer_Vectors_io;
-- use Standard_Integer_Vectors_io;

with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Supports_of_Polynomial_Systems;
with Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;    use Floating_Mixed_Subdivisions_io;
with Integer_Mixed_Subdivisions;
with Integer_Mixed_Subdivisions_io;     use Integer_Mixed_Subdivisions_io;
with Floating_Lifting_Functions;
with Floating_Integer_Convertors;
with Integer_Lifting_Utilities;         use Integer_Lifting_Utilities;
with Floating_Lifting_Utilities;        use Floating_Lifting_Utilities;
with Induced_Permutations;
with Mixed_Volume_Computation;          use Mixed_Volume_Computation;
with PHCpack_Operations;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Standard_LaurSys_Container;
with DoblDobl_LaurSys_Container;
with QuadDobl_LaurSys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Cells_Container;
with Integer_Cells_Container;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

function use_celcon ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

  function Select_Point ( L : Lists_of_Integer_Vectors.List; k : natural32 )
                        return Standard_Integer_Vectors.Link_to_Vector is

  -- DESCRIPTION :
  --   Returns the k-th point from the list L.

    use Lists_of_Integer_Vectors;

    res : Standard_Integer_Vectors.Link_to_Vector;
    tmp : List := L;

  begin
    for i in 1..(k-1) loop
      exit when Is_Null(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    if not Is_Null(tmp)
     then res := Head_Of(tmp);
    end if;
    return res;
  end Select_Point;

  function Select_Point ( L : Lists_of_Floating_Vectors.List; k : natural32 )
                        return Standard_Floating_Vectors.Link_to_Vector is

  -- DESCRIPTION :
  --   Returns the k-th point from the list L.

    use Lists_of_Floating_Vectors;

    res : Standard_Floating_Vectors.Link_to_Vector;
    tmp : List := L;

  begin
    for i in 1..(k-1) loop
      exit when Is_Null(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    if not Is_Null(tmp)
     then res := Head_Of(tmp);
    end if;
    return res;
  end Select_Point;

  function Cells_Read_Floating_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Prompts the user for a regular mixed cell configuration,
  --   induced by a floating-point lifting function,
  --   and stores the mixed cells in the container.

    use Arrays_of_Floating_Vector_Lists;
    use Floating_mixed_Subdivisions;
  
    file : file_type;
    n,m : natural32;
    r : integer32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    sub : Mixed_Subdivision;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Read_Floating_Mixed_Cells ...");
    end if;
    new_line;
    put_line("Reading the file name for a mixed-cell configuration.");
    Read_Name_and_Open_File(file);
    get(file,n,m,mix,sub);
    close(file);
    r := mix'last;
    lif := new Array_of_Lists'(Lifted_Supports(r,sub));
    Cells_Container.Initialize(mix,lif,sub);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Read_Floating_Mixed_Cells.");
      end if;
      return 80;
  end Cells_Read_Floating_Mixed_Cells;

  function Cells_Write_Floating_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Writes the mixed cells stored in the container
  --   for configurations induced by floating-point lifting.

    use Floating_mixed_Subdivisions;

    mcc : Mixed_Subdivision := Cells_Container.Retrieve;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    n,mv : natural32;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Write_Floating_Mixed_Cells ...");
    end if;
    if not Is_Null(mcc) then
      n := Cells_Container.Dimension;
      mix := Cells_Container.Type_of_Mixture;
      put(standard_output,n-1,mix.all,mcc,mv);
      put("The mixed volume is "); put(mv,1); put_line(".");
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Write_Floating_Mixed_Cells.");
      end if;
      return 81;
  end Cells_Write_Floating_Mixed_Cells;

  function Cells_Number_of_Floating_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the number of mixed cells stored in the container,
  --   induced by a floating-point lifting.

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Number_of_Floating_Mixed_Cells ...");
    end if;
    Assign(integer32(Cells_Container.Length),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Number_of_Floating_Mixed_Cells.");
      end if;
      return 82;
  end Cells_Number_of_Floating_Mixed_Cells;

  function Cells_Dimension_of_Floating_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the dimension of the mixed cells stored in the container,
  --   induced by a floating-point lifting.

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Dimension_of_Floating_Mixed_Cells ...");
    end if;
    Assign(integer32(Cells_Container.Dimension),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Dimension_of_Floating_Mixed_Cells.");
      end if;
      return 83;
  end Cells_Dimension_of_Floating_Mixed_Cells;

  function Cells_Floating_Mixture
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Integer_Vectors;
    mix : Link_to_Vector;
    r : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Floating_Mixture ...");
    end if;
    mix := Cells_Container.Type_of_Mixture;
    if mix /= null
     then r := mix'last; Assign(mix.all,b);
     else r := 0;
    end if;
    Assign(r,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Floating_Mixture.");
      end if;
      return 84;
  end Cells_Floating_Mixture;

  function Job5 return integer32 is -- return cardinalities of supports

    use Lists_of_Floating_Vectors;
    use Arrays_of_Floating_Vector_Lists;

    lif : Link_to_Array_of_Lists;
    r : integer32;

  begin
    lif := Cells_Container.Lifted_Supports;
    if lif /= null then
      r := lif'last;
      declare
        nl : Standard_Integer_Vectors.Vector(1..r);
      begin
        for i in nl'range loop
          nl(i) := integer32(Length_Of(lif(i)));
        end loop;
        Assign(nl,b);
      end;
    else
      r := 0;
    end if;
    Assign(r,a);
    return 0;
  end Job5;

  function Job6 return integer32 is -- return point in support

    use Arrays_of_Floating_Vector_Lists;

    lif : Link_to_Array_of_Lists;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    k_a : constant integer32 := integer32(v_a(v_a'first));
    k_b : constant natural32 := natural32(v_b(v_b'first));
    use Standard_Floating_Vectors;
    lpt : Link_to_Vector;

  begin
    lif := Cells_Container.Lifted_Supports;
   -- put("retrieving point "); put(k_b,1);
   -- put(" from list "); put(k_a,1); put_line(" ...");
    if lif /= null then
      lpt := Select_Point(lif(k_a),k_b);
      if lpt /= null then
        Assign(lpt.all,c);     
        return 0;
      end if;
    end if;
    return 1;
  end Job6;

  function Job7 return integer32 is -- return innner normal

    use Floating_mixed_Subdivisions;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    mic : Mixed_Cell;
    fail : boolean;

  begin
    Cells_Container.Retrieve(k,mic,fail);
    if fail
     then return 1;
     else Assign(mic.nor.all,c); return 0;
    end if;
  end Job7;

  function Job8 return integer32 is -- return #elements in cell

    use Lists_of_Floating_Vectors;
    use Floating_Mixed_Subdivisions;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    mic : Mixed_Cell;
    fail : boolean;

  begin
    Cells_Container.Retrieve(k,mic,fail);
    if not fail then
      declare
        nl : Standard_Integer_Vectors.Vector(mic.pts'range);
      begin
        for i in mic.pts'range loop
          nl(i) := integer32(Length_Of(mic.pts(i)));
        end loop;
        Assign(nl,b);
      end;
      return 0;
    end if;
    return 1;
  end Job8;

  function Job9 return integer32 is -- return point in a cell

    use Floating_Mixed_Subdivisions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    i : constant integer32 := integer32(v_b(v_b'first));
    j : constant natural32 := natural32(v_b(v_b'first+1));
    mic : Mixed_Cell;
    fail : boolean;
    use Standard_Floating_Vectors;
    lpt : Link_to_Vector;

  begin
    Cells_Container.Retrieve(k,mic,fail);
    if not fail then
      lpt := Select_Point(mic.pts(i),j);
      if lpt /= null 
       then Assign(lpt.all,c); return 0;
      end if;
    end if;
    return 1;
  end Job9;

  function Job10 return integer32 is -- return mixed volume of a cell

    use Floating_Mixed_Subdivisions;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    mix : constant Standard_Integer_Vectors.Link_to_Vector
        := Cells_Container.Type_of_Mixture;
    mic : Mixed_Cell;
    fail : boolean;
    mv : natural32;
    n : constant integer32 := integer32(Cells_Container.Dimension)-1;
    use Standard_Integer_Vectors;

  begin
    Cells_Container.Retrieve(k,mic,fail);
    if fail or mix = null then
      return 1;
    else
      Mixed_Volume(n,mix.all,mic,mv);
      Assign(integer32(mv),b);
      return 0;
    end if;
  end Job10;

  function Job11 return integer32 is -- sets type of mixture

    v : constant C_Integer_Array := C_intarrs.Value(a);
    r : constant integer32 := integer32(v(v'first));
    mix : constant Standard_Integer_Vectors.Link_to_Vector
        := new Standard_Integer_Vectors.Vector(1..r);
    r1 : constant Interfaces.C.size_t := Interfaces.C.size_t(r-1);
    mbv : constant C_Integer_Array(0..r1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(r));

  begin
    for i in 0..r-1 loop
      mix(i+1) := integer32(mbv(Interfaces.C.size_t(i)));
    end loop;
    Cells_Container.Initialize(mix);
    return 0;
  end Job11;

  function Job12 return integer32 is -- append point to a support

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    n : constant integer32 := integer32(vb(vb'first));
    x : Standard_Floating_Vectors.Vector(1..n);
    fail : boolean;

  begin
    Assign(natural32(n),c,x);
    fail := Cells_Container.Append_to_Support(k,x);
    if fail
     then return 1;
     else return 0;
    end if;
  end Job12;

  function Job13 return integer32 is -- append mixed cell to container

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    r : constant integer32 := integer32(va(va'first));
    n : constant integer32 := integer32(va(va'first+1));
    d : constant Interfaces.C.size_t := Interfaces.C.size_t(va(2));
    vb : constant C_Integer_Array(0..d-1)
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(d));
    m : constant integer32 := integer32(vb(0));
    x : Standard_Floating_Vectors.Vector(1..n);
    cnt : Standard_Integer_Vectors.Vector(1..r);
    lab : Standard_Integer_Vectors.Vector(1..m);

  begin
   -- put("  r = "); put(r,1);
   -- put("  n = "); put(n,1);
   -- put("  d = "); put(natural(d),1); new_line;
    for i in cnt'range loop
      cnt(i) := integer32(vb(Interfaces.C.size_t(i)));
    end loop;
    for i in lab'range loop
      lab(i) := integer32(vb(Interfaces.C.size_t(i+r)));
    end loop;
   -- put("the number of points in each support :"); put(cnt); new_line;
   -- put("the labels of points in each support :"); put(lab); new_line;
    Assign(natural32(n),c,x);
    Cells_Container.Append_Mixed_Cell(cnt,lab,x);
    return 0;
  end Job13;

  function Job15 return integer32 is -- retrieve a mixed cell

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    cnt,lab : Standard_Integer_Vectors.Link_to_Vector;
    normal : Standard_Floating_Vectors.Link_to_Vector;
    fail : boolean;

  begin
    Cells_Container.Retrieve_Mixed_Cell(k,fail,cnt,lab,normal);
    if fail then
      return 1;
    else
      declare
        sum : constant integer32 := Standard_Integer_Vectors.Sum(cnt.all);
        cntlab : Standard_Integer_Vectors.Vector(1..cnt'last+lab'last+1);
        ind : integer32 := 0;
      begin
        ind := ind + 1;
        cntlab(ind) := sum;
        for i in cnt'range loop
          ind := ind + 1;
          cntlab(ind) := cnt(i);
        end loop;
        for i in lab'range loop
          ind := ind + 1;
          cntlab(ind) := lab(i);
        end loop;
        Assign(cntlab,b);
      end;
      Assign(normal.all,c);
      return 0;
    end if;
  end Job15;

  function Cells_Make_Standard_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Generates a random coefficient polynomial system in double
  --   precision using the mixtures and supports in the cells container.
  --   The system is stored in the cells container.

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Make_Standard_Coefficient_System ...");
    end if;
    Cells_Container.Generate_Random_Standard_Coefficient_System;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Make_Standard_Coefficient_System.");
      end if;
      return 96;
  end Cells_Make_Standard_Coefficient_System;

  function Cells_Make_DoblDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Generates a random coefficient polynomial system in double double
  --   precision using the mixtures and supports in the cells container.
  --   The system is stored in the cells container.

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Make_DoblDobl_Coefficient_System ...");
    end if;
    Cells_Container.Generate_Random_DoblDobl_Coefficient_System;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Make_DoblDobl_Coefficient_System.");
      end if;
      return 460;
  end Cells_Make_DoblDobl_Coefficient_System;

  function Cells_Make_QuadDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Generates a random coefficient polynomial system in quad double
  --   precision using the mixtures and supports in the cells container.
  --   The system is stored in the cells container.

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Make_QuadDobl_Coefficient_System ...");
    end if;
    Cells_Container.Generate_Random_QuadDobl_Coefficient_System;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Make_QuadDobl_Coefficient_System.");
      end if;
      return 470;
  end Cells_Make_QuadDobl_Coefficient_System;

  function Cells_Read_Standard_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Prompts the user for a coefficient polynomial system
  --   in double precision and stores it in the container.

    q : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Read_Standard_Coefficient_System ...");
    end if;
    new_line;
    put_line("Reading a random coefficient polynomial system ...");
    get(q);
    Cells_Container.Initialize_Random_Standard_Coefficient_System(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Read_Standard_Coefficient_System.");
      end if;
      return 97;
  end Cells_Read_Standard_Coefficient_System;

  function Cells_Read_DoblDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Prompts the user for a coefficient polynomial system
  --   in double double precision and stores it in the container.

    q : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Read_DoblDobl_Coefficient_System ...");
    end if;
    new_line;
    put_line("Reading a random coefficient polynomial system ...");
    get(q);
    Cells_Container.Initialize_Random_DoblDobl_Coefficient_System(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Read_DoblDobl_Coefficient_System.");
      end if;
      return 461;
  end Cells_Read_DoblDobl_Coefficient_System;

  function Cells_Read_QuadDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Prompts the user for a coefficient polynomial system
  --   in quad double precision and stores it in the container.

    q : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Read_QuadDobl_Coefficient_System ...");
    end if;
    new_line;
    put_line("Reading a random coefficient polynomial system ...");
    get(q);
    Cells_Container.Initialize_Random_QuadDobl_Coefficient_System(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Read_QuadDobl_Coefficient_System.");
      end if;
      return 471;
  end Cells_Read_QuadDobl_Coefficient_System;

  function Cells_Write_Standard_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Writes the coefficient system in double precision
  --   to standard output or to the defined output file.

    q : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Cells_Container.Retrieve_Random_Standard_Coefficient_System;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Write_Standard_Coefficient_System ...");
    end if;
    if PHCpack_Operations.Is_File_Defined then
      put_line(PHCpack_Operations.output_file,q.all);
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE SOLUTIONS :");
    else
      put_line(standard_output,q.all);
      new_line(standard_output);
      put_line(standard_output,"THE SOLUTIONS :");
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Write_Standard_Coefficient_System.");
      end if;
      return 98;
  end Cells_Write_Standard_Coefficient_System;

  function Cells_Write_DoblDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Writes the coefficient system in double double precision
  --   to standard output or to the defined output file.

    q : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := Cells_Container.Retrieve_Random_DoblDobl_Coefficient_System;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Write_DoblDobl_Coefficient_System ...");
    end if;
    if PHCpack_Operations.Is_File_Defined then
      put_line(PHCpack_Operations.output_file,q.all);
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE SOLUTIONS :");
    else
      put_line(standard_output,q.all);
      new_line(standard_output);
      put_line(standard_output,"THE SOLUTIONS :");
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Write_DoblDobl_Coefficient_System.");
      end if;
      return 462;
  end Cells_Write_DoblDobl_Coefficient_System;

  function Cells_Write_QuadDobl_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Writes the coefficient system in quad double precision
  --   to standard output or to the defined output file.

    q : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := Cells_Container.Retrieve_Random_QuadDobl_Coefficient_System;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Write_QuadDobl_Coefficient_System ...");
    end if;
    if PHCpack_Operations.Is_File_Defined then
      put_line(PHCpack_Operations.output_file,q.all);
      new_line(PHCpack_Operations.output_file);
      put_line(PHCpack_Operations.output_file,"THE SOLUTIONS :");
    else
      put_line(standard_output,q.all);
      new_line(standard_output);
      put_line(standard_output,"THE SOLUTIONS :");
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Write_QuadDobl_Coefficient_System.");
      end if;
      return 472;
  end Cells_Write_QuadDobl_Coefficient_System;

  function Cells_Standard_System_into_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Copies the coefficient system in double precision
  --   into the systems container.

    q : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Cells_Container.Retrieve_Random_Standard_Coefficient_System;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Standard_System_into_Container ...");
    end if;
    Standard_PolySys_Container.Initialize(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Standard_System_into_Container.");
      end if;
      return 99;
  end Cells_Standard_System_into_Container;

  function Cells_DoblDobl_System_into_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Copies the coefficient system in double double precision
  --   into the systems container.

    q : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := Cells_Container.Retrieve_Random_DoblDobl_Coefficient_System;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_DoblDobl_System_into_Container ...");
    end if;
    DoblDobl_PolySys_Container.Initialize(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_DoblDobl_System_into_Container.");
      end if;
      return 463;
  end Cells_DoblDobl_System_into_Container;

  function Cells_QuadDobl_System_into_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Copies the coefficient system in quad double precision
  --   into the systems container.

    q : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := Cells_Container.Retrieve_Random_QuadDobl_Coefficient_System;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_QuadDobl_System_into_Container ...");
    end if;
    QuadDobl_PolySys_Container.Initialize(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_QuadDobl_System_into_Container.");
      end if;
      return 473;
  end Cells_QuadDobl_System_into_Container;

  function Cells_Standard_System_from_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Copies the system in the container for double precision
  --   into the cells container.

    q : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Standard_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Standard_System_from_Container ...");
    end if;
    Cells_Container.Initialize_Random_Standard_Coefficient_System(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Standard_System_from_Container.");
      end if;
      return 100;
  end Cells_Standard_System_from_Container;

  function Cells_DoblDobl_System_from_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Copies the system in the container for double double precision
  --   into the cells container.

    q : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := DoblDobl_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_DoblDobl_System_from_Container ...");
    end if;
    Cells_Container.Initialize_Random_DoblDobl_Coefficient_System(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_DoblDobl_System_from_Container.");
      end if;
      return 464;
  end Cells_DoblDobl_System_from_Container;

  function Cells_QuadDobl_System_from_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Copies the system in the container for quad double precision
  --   into the cells container.

    q : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := QuadDobl_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_QuadDobl_System_from_Container ...");
    end if;
    Cells_Container.Initialize_Random_QuadDobl_Coefficient_System(q.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_QuadDobl_System_from_Container.");
      end if;
      return 474;
  end Cells_QuadDobl_System_from_Container;

  function Cells_Standard_Polyhedral_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Makes a polyhedral homotopy to solve a random system
  --   in double precision.

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Standard_Polyhedral_Homotopy ...");
    end if;
    Cells_Container.Standard_Polyhedral_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Standard_Polyhedral_Homotopy.");
      end if;
      return 101;
  end Cells_Standard_Polyhedral_Homotopy;

  function Cells_DoblDobl_Polyhedral_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Makes a polyhedral homotopy to solve a random system
  --   in double double precision.

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_DoblDobl_Polyhedral_Homotopy ...");
    end if;
    Cells_Container.DoblDobl_Polyhedral_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_DoblDobl_Polyhedral_Homotopy.");
      end if;
      return 465;
  end Cells_DoblDobl_Polyhedral_Homotopy;

  function Cells_QuadDobl_Polyhedral_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Makes a polyhedral homotopy to solve a random system
  --   in quad double precision.

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_QuadDobl_Polyhedral_Homotopy ...");
    end if;
    Cells_Container.QuadDobl_Polyhedral_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_QuadDobl_Polyhedral_Homotopy.");
      end if;
      return 475;
  end Cells_QuadDobl_Polyhedral_Homotopy;

  function Cells_Standard_Start_Solve 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Solves the start system which corresponds to a mixed cell,
  --   in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to the mixed cells;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the number of solutions.

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    mv : natural32;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Standard_Start_Solve ...");
    end if;
   -- put("Entering Job22 with k = "); put(k,1); put_line(" ...");
    Cells_Container.Solve_Standard_Start_System(k,mv);
    Assign(integer32(mv),b);
   -- put("... leaving Job22 with mv = "); put(mv,1); put_line(".");
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Standard_Start_Solve.");
      end if;
      return 102;
  end Cells_Standard_Start_Solve;

  function Cells_DoblDobl_Start_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Solves the start system which corresponds to a mixed cell,
  --   in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to the mixed cells;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the number of solutions.

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    mv : natural32;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_DoblDobl_Start_Solve ...");
    end if;
   -- put("Entering Job22 with k = "); put(k,1); put_line(" ...");
    Cells_Container.Solve_DoblDobl_Start_System(k,mv);
    Assign(integer32(mv),b);
   -- put("... leaving Job22 with mv = "); put(mv,1); put_line(".");
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_DoblDobl_Start_Solve.");
      end if;
      return 466;
  end Cells_DoblDobl_Start_Solve;

  function Cells_QuadDobl_Start_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Solves the start system which corresponds to a mixed cell,
  --   in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to the mixed cells;
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   b       in b[0] is the number of solutions.

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    mv : natural32;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_QuadDobl_Start_Solve ...");
    end if;
   -- put("Entering Job22 with k = "); put(k,1); put_line(" ...");
    Cells_Container.Solve_QuadDobl_Start_System(k,mv);
    Assign(integer32(mv),b);
   -- put("... leaving Job22 with mv = "); put(mv,1); put_line(".");
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_QuadDobl_Start_Solve.");
      end if;
      return 476;
  end Cells_QuadDobl_Start_Solve;

  function Cells_Standard_Track_One_Path
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Tracks one solution path in double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to a mixed cell;
  --   b       in b[0] is the index to the start solution;
  --           in b[1] is the output code for the trackers;
  --   vrblvl  is the verbose level.

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    i : constant natural32 := natural32(vb(vb'first));
    otp : constant natural32 := natural32(vb(vb'first+1));

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Standard_Track_One_Path ...");
    end if;
    Cells_Container.Track_Standard_Solution_Path(k,i,otp);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Standard_Track_One_Path.");
      end if;
      return 103;
  end Cells_Standard_Track_One_Path;

  function Cells_DoblDobl_Track_One_Path
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Tracks one solution path in double double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to a mixed cell;
  --   b       in b[0] is the index to the start solution;
  --           in b[1] is the output code for the trackers;
  --   vrblvl  is the verbose level.

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    i : constant natural32 := natural32(vb(vb'first));
    otp : constant natural32 := natural32(vb(vb'first+1));

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_DoblDobl_Track_One_Path ...");
    end if;
    Cells_Container.Track_DoblDobl_Solution_Path(k,i,otp);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_DoblDobl_Track_One_Path.");
      end if;
      return 467;
  end Cells_DoblDobl_Track_One_Path;

  function Cells_QuadDobl_Track_One_Path
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Tracks one solution path in quad double precision.

  -- ON ENTRY :
  --   a       in a[0] is the index to a mixed cell;
  --   b       in b[0] is the index to the start solution;
  --           in b[1] is the output code for the trackers;
  --   vrblvl  is the verbose level.

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    i : constant natural32 := natural32(vb(vb'first));
    otp : constant natural32 := natural32(vb(vb'first+1));

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_QuadDobl_Track_One_Path ...");
    end if;
    Cells_Container.Track_QuadDobl_Solution_Path(k,i,otp);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_QuadDobl_Track_One_Path.");
      end if;
      return 477;
  end Cells_QuadDobl_Track_One_Path;

  function Job24 return integer32 is -- copy target solution to st container

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : Standard_Complex_Solutions.Link_to_Solution;

    use Standard_Complex_Solutions;

  begin
    ls := Cells_Container.Retrieve_Standard_Target_Solution(k,i);
    if ls /= null
     then Standard_Solutions_Container.Append(ls.all);
    end if;
    return 0;
  end Job24;

  function Job34 return integer32 is -- copy target solution to dd container

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

    use DoblDobl_Complex_Solutions;

  begin
    ls := Cells_Container.Retrieve_DoblDobl_Target_Solution(k,i);
    if ls /= null
     then DoblDobl_Solutions_Container.Append(ls.all);
    end if;
    return 0;
  end Job34;

  function Job44 return integer32 is -- copy target solution to qd container

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

    use QuadDobl_Complex_Solutions;

  begin
    ls := Cells_Container.Retrieve_QuadDobl_Target_Solution(k,i);
    if ls /= null
     then QuadDobl_Solutions_Container.Append(ls.all);
    end if;
    return 0;
  end Job44;

  function Job48 return integer32 is -- copy start solution to st container

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    ls := Cells_Container.Retrieve_Standard_Start_Solution(k,i);
    Standard_Solutions_Container.Append(ls.all);
    return 0;
  end Job48;

  function Job49 return integer32 is -- copy start solution to dd container

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    ls := Cells_Container.Retrieve_DoblDobl_Start_Solution(k,i);
    DoblDobl_Solutions_Container.Append(ls.all);
    return 0;
  end Job49;

  function Job50 return integer32 is -- copy start solution to qd container

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    ls := Cells_Container.Retrieve_QuadDobl_Start_Solution(k,i);
    QuadDobl_Solutions_Container.Append(ls.all);
    return 0;
  end Job50;

  function Job25 return integer32 is -- permute system in st container 

    use Floating_mixed_Subdivisions;

    mixsub : constant Mixed_Subdivision := Cells_Container.Retrieve;
    mix : constant Standard_Integer_Vectors.Link_to_Vector
        := Cells_Container.Type_of_Mixture;
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    lq : constant Standard_Complex_Laur_Systems.Link_to_Laur_Sys
       := Standard_LaurSys_Container.Retrieve;
    use Standard_Complex_Poly_Systems,Standard_Complex_Laur_Systems;

  begin
    if lp /= null then
      declare
        sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(lp'range)
            := Supports_of_Polynomial_Systems.Create(lp.all);
        stlb : constant double_float 
             := Floating_Lifting_Functions.Lifting_Bound(lp.all);
        fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
           := Floating_Integer_Convertors.Convert(sup);
        ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
           := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
        ip : Standard_Integer_Vectors.Vector(fs'range);
      begin
        Induced_Permutations.Remove_Artificial_Origin(ls,stlb);
        ip := Induced_Permutations.Permutation(fs,ls,mix.all);
       -- put("the permutation : "); put(ip); new_line;
       -- put_line("before permute :"); put(lp.all);
        Induced_Permutations.Permute(ip,lp.all);
       -- put_line("after permute :"); put(lp.all);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(ls);
        Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
      end;
    end if;
    if lq /= null then
      declare
        sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(lq'range)
            := Supports_of_Polynomial_Systems.Create(lq.all);
        fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
           := Floating_Integer_Convertors.Convert(sup);
        ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
           := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
        ip : Standard_Integer_Vectors.Vector(fs'range);
      begin
        ip := Induced_Permutations.Permutation(fs,ls,mix.all);
        Induced_Permutations.Permute(ip,lq.all);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(ls);
        Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
      end;
    end if;
    return 0;
  end Job25;

  function Job35 return integer32 is -- permute system in dd container 

    use Floating_mixed_Subdivisions;

    mixsub : constant Mixed_Subdivision := Cells_Container.Retrieve;
    mix : constant Standard_Integer_Vectors.Link_to_Vector
        := Cells_Container.Type_of_Mixture;
    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    lq : constant DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys
       := DoblDobl_LaurSys_Container.Retrieve;
    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Laur_Systems;

  begin
    if lp /= null then
      declare
        sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(lp'range)
            := Supports_of_Polynomial_Systems.Create(lp.all);
        stlb : constant double_float 
             := Floating_Lifting_Functions.Lifting_Bound(lp.all);
        fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
           := Floating_Integer_Convertors.Convert(sup);
        ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
           := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
        ip : Standard_Integer_Vectors.Vector(fs'range);
      begin
        Induced_Permutations.Remove_Artificial_Origin(ls,stlb);
        ip := Induced_Permutations.Permutation(fs,ls,mix.all);
       -- put("the permutation : "); put(ip); new_line;
       -- put_line("before permute :"); put(lp.all);
        Induced_Permutations.Permute(ip,lp.all);
       -- put_line("after permute :"); put(lp.all);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(ls);
        Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
      end;
    end if;
    if lq /= null then
      declare
        sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(lq'range)
            := Supports_of_Polynomial_Systems.Create(lq.all);
        fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
           := Floating_Integer_Convertors.Convert(sup);
        ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
           := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
        ip : Standard_Integer_Vectors.Vector(fs'range);
      begin
        ip := Induced_Permutations.Permutation(fs,ls,mix.all);
        Induced_Permutations.Permute(ip,lq.all);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(ls);
        Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
      end;
    end if;
    return 0;
  end Job35;

  function Job45 return integer32 is -- permute system in qd container 

    use Floating_mixed_Subdivisions;

    mixsub : constant Mixed_Subdivision := Cells_Container.Retrieve;
    mix : constant Standard_Integer_Vectors.Link_to_Vector
        := Cells_Container.Type_of_Mixture;
    lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    lq : constant QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys
       := QuadDobl_LaurSys_Container.Retrieve;
    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Laur_Systems;

  begin
    if lp /= null then
      declare
        sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(lp'range)
            := Supports_of_Polynomial_Systems.Create(lp.all);
        stlb : constant double_float 
             := Floating_Lifting_Functions.Lifting_Bound(lp.all);
        fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
           := Floating_Integer_Convertors.Convert(sup);
        ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
           := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
        ip : Standard_Integer_Vectors.Vector(fs'range);
      begin
        Induced_Permutations.Remove_Artificial_Origin(ls,stlb);
        ip := Induced_Permutations.Permutation(fs,ls,mix.all);
       -- put("the permutation : "); put(ip); new_line;
       -- put_line("before permute :"); put(lp.all);
        Induced_Permutations.Permute(ip,lp.all);
       -- put_line("after permute :"); put(lp.all);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(ls);
        Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
      end;
    end if;
    if lq /= null then
      declare
        sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(lq'range)
            := Supports_of_Polynomial_Systems.Create(lq.all);
        fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
           := Floating_Integer_Convertors.Convert(sup);
        ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
           := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
        ip : Standard_Integer_Vectors.Vector(fs'range);
      begin
        ip := Induced_Permutations.Permutation(fs,ls,mix.all);
        Induced_Permutations.Permute(ip,lq.all);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
        Arrays_of_Floating_Vector_Lists.Deep_Clear(ls);
        Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
      end;
    end if;
    return 0;
  end Job45;

  function Job46 return integer32 is -- mixed volume computation

    mv : constant natural32 := Cells_Container.Mixed_Volume;

  begin
    Assign(integer32(mv),a);
    return 0;
  end Job46;

  function Job47 return integer32 is -- initialize #distinct supports

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbr : constant natural32 := natural32(v_a(v_a'first));

  begin
    Cells_Container.Initialize_Supports(nbr);
    return 0;
  end Job47;

  function Cells_Read_Integer_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Prompts the user for a regular mixed cell configuration,
  --   induced by an integer valued lifting function,
  --   and stores the mixed cells in the container.

    use Arrays_of_Integer_Vector_Lists;
    use Integer_Mixed_Subdivisions;
  
    file : file_type;
    n,m : natural32;
    r : integer32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    sub : Mixed_Subdivision;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Read_Integer_Mixed_Cells ...");
    end if;
    new_line;
    put_line("Reading the file name for a mixed-cell configuration.");
    Read_Name_and_Open_File(file);
    get(file,n,m,mix,sub);
    close(file);
    r := mix'last;
    lif := new Array_of_Lists'(Lifted_Supports(r,sub));
    Integer_Cells_Container.Initialize(mix,lif,sub);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Read_Integer_Mixed_Cells.");
      end if;
      return 741;
  end Cells_Read_Integer_Mixed_Cells;

  function Cells_Write_Integer_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Writes the mixed cells stored in the container
  --   for configurations induced by floating-point lifting.

    use Integer_Mixed_Subdivisions;

    mcc : Mixed_Subdivision := Integer_Cells_Container.Retrieve;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    n,mv : natural32;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Write_Integer_Mixed_Cells ...");
    end if;
    if not Is_Null(mcc) then
      n := Integer_Cells_Container.Dimension;
      mix := Integer_Cells_Container.Type_of_Mixture;
      put(standard_output,n-1,mix.all,mcc,mv);
      put("The mixed volume is "); put(mv,1); put_line(".");
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Write_Integer_Mixed_Cells.");
      end if;
      return 742;
  end Cells_Write_Integer_Mixed_Cells;

  function Cells_Number_of_Integer_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the number of mixed cells stored in the container,
  --   induced by an integer-valued lifting.

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Number_of_Integer_Mixed_Cells ...");
    end if;
    Assign(integer32(Integer_Cells_Container.Length),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Number_of_Integer_Mixed_Cells.");
      end if;
      return 743;
  end Cells_Number_of_Integer_Mixed_Cells;

  function Cells_Dimension_of_Integer_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the dimension of the mixed cells stored in the container,
  --   induced by an integer-valued lifting.

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Dimension_of_Integer_Mixed_Cells ...");
    end if;
    Assign(integer32(Integer_Cells_Container.Dimension),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Dimension_of_Integer_Mixed_Cells.");
      end if;
      return 744;
  end Cells_Dimension_of_Integer_Mixed_Cells;

  function Job55 return integer32 is -- return type of mixture

    use Standard_Integer_Vectors;
    mix : Link_to_Vector;
    r : integer32;

  begin
    mix := Integer_Cells_Container.Type_of_Mixture;
    if mix /= null
     then r := mix'last; Assign(mix.all,b);
     else r := 0;
    end if;
    Assign(r,a);
    return 0;
  end Job55;

  function Job56 return integer32 is -- return cardinalities of supports

    use Lists_of_Integer_Vectors;
    use Arrays_of_Integer_Vector_Lists;

    lif : Link_to_Array_of_Lists;
    r : integer32;

  begin
    lif := Integer_Cells_Container.Lifted_Supports;
    if lif /= null then
      r := lif'last;
      declare
        nl : Standard_Integer_Vectors.Vector(1..r);
      begin
        for i in nl'range loop
          nl(i) := integer32(Length_Of(lif(i)));
        end loop;
        Assign(nl,b);
      end;
    else
      r := 0;
    end if;
    Assign(r,a);
    return 0;
  end Job56;

  function Job57 return integer32 is -- return point in support

    use Arrays_of_Integer_Vector_Lists;

    lif : Link_to_Array_of_Lists;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    k_a : constant integer32 := integer32(v_a(v_a'first));
    k_b : constant natural32 := natural32(v_b(v_b'first));
    use Standard_Integer_Vectors;
    lpt : Link_to_Vector;

  begin
    lif := Integer_Cells_Container.Lifted_Supports;
   -- put("retrieving point "); put(k_b,1);
   -- put(" from list "); put(k_a,1); put_line(" ...");
    if lif /= null then
      lpt := Select_Point(lif(k_a),k_b);
      if lpt /= null then
        declare
          fltlpt : Standard_Floating_Vectors.Vector(lpt'range);
        begin
          for k in lpt'range loop
            fltlpt(k) := double_float(lpt(k));
          end loop;
          Assign(fltlpt,c);     
        end;
        return 0;
      end if;
    end if;
    return 1;
  end Job57;

  function Job58 return integer32 is -- return innner normal

    use Integer_Mixed_Subdivisions;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    mic : Mixed_Cell;
    fail : boolean;

  begin
    Integer_Cells_Container.Retrieve(k,mic,fail);
    if fail then
      return 1;
    else
      declare
        fltnor : Standard_Floating_Vectors.Vector(mic.nor'range);
      begin
        for k in mic.nor'range loop
          fltnor(k) := double_float(mic.nor(k));
        end loop;
        Assign(fltnor,c); return 0;
      end;
    end if;
  end Job58;

  function Job59 return integer32 is -- return #elements in cell

    use Lists_of_Integer_Vectors;
    use Integer_Mixed_Subdivisions;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    mic : Mixed_Cell;
    fail : boolean;

  begin
    Integer_Cells_Container.Retrieve(k,mic,fail);
    if not fail then
      declare
        nl : Standard_Integer_Vectors.Vector(mic.pts'range);
      begin
        for i in mic.pts'range loop
          nl(i) := integer32(Length_Of(mic.pts(i)));
        end loop;
        Assign(nl,b);
      end;
      return 0;
    end if;
    return 1;
  end Job59;

  function Job60 return integer32 is -- return point in a cell

    use Integer_Mixed_Subdivisions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    i : constant integer32 := integer32(v_b(v_b'first));
    j : constant natural32 := natural32(v_b(v_b'first+1));
    mic : Mixed_Cell;
    fail : boolean;
    use Standard_Integer_Vectors;
    lpt : Link_to_Vector;

  begin
    Integer_Cells_Container.Retrieve(k,mic,fail);
    if not fail then
      lpt := Select_Point(mic.pts(i),j);
      if lpt /= null then
        declare
          fltlpt : Standard_Floating_Vectors.Vector(lpt'range);
        begin
          for k in lpt'range loop
            fltlpt(k) := double_float(lpt(k));
          end loop;
          Assign(fltlpt,c);
        end;
        return 0;
      end if;
    end if;
    return 1;
  end Job60;

  function Job61 return integer32 is -- return mixed volume of a cell

    use Integer_Mixed_Subdivisions;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    mix : constant Standard_Integer_Vectors.Link_to_Vector
        := Integer_Cells_Container.Type_of_Mixture;
    mic : Mixed_Cell;
    fail : boolean;
    mv : natural32;
    n : constant integer32 := integer32(Integer_Cells_Container.Dimension)-1;

    use Standard_Integer_Vectors;

  begin
   -- put("The dimension : "); put(n,1); new_line;
   -- put("retrieving mixed cell "); put(k,1); put_line(" ...");
    Integer_Cells_Container.Retrieve(k,mic,fail);
    if fail or mix = null then
     -- put_line("failure!");
      return 1;
    else
     -- put("The type of mixture :"); put(mix); new_line;
      Mixed_Volume(n,mix.all,mic,mv);
     -- put("The mixed volume mv is "); put(mv,1); put_line(".");
      Assign(integer32(mv),b);
      return 0;
    end if;
  end Job61;

  function Job62 return integer32 is -- initialize #distinct supports

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbr : constant natural32 := natural32(v_a(v_a'first));

  begin
    Integer_Cells_Container.Initialize_Supports(nbr);
    return 0;
  end Job62;

  function Job63 return integer32 is -- sets type of mixture

    v : constant C_Integer_Array := C_intarrs.Value(a);
    r : constant integer32 := integer32(v(v'first));
    mix : Standard_Integer_Vectors.Vector(1..r);
    lnkmix : Standard_Integer_Vectors.Link_to_Vector;
    r1 : constant Interfaces.C.size_t := Interfaces.C.size_t(r-1);
    mbv : constant C_Integer_Array(0..r1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(r));

  begin
    for i in 0..r-1 loop
      mix(i+1) := integer32(mbv(Interfaces.C.size_t(i)));
    end loop;
    lnkmix := new Standard_Integer_Vectors.Vector'(mix);
    Integer_Cells_Container.Initialize(lnkmix);
    return 0;
  end Job63;

  function Job64 return integer32 is -- append point to a support

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    n : constant integer32 := integer32(vb(vb'first));
    x : Standard_Integer_Vectors.Vector(1..n);
    fltx : Standard_Floating_Vectors.Vector(1..n);
    fail : boolean;

  begin
    Assign(natural32(n),c,fltx);
    for k in 1..n loop
      x(k) := integer32(fltx(k));
    end loop;
    fail := Integer_Cells_Container.Append_to_Support(k,x);
    if fail
     then return 1;
     else return 0;
    end if;
  end Job64;

  function Job65 return integer32 is -- append mixed cell to container

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    r : constant integer32 := integer32(va(va'first));
    n : constant integer32 := integer32(va(va'first+1));
    d : constant Interfaces.C.size_t := Interfaces.C.size_t(va(2));
    vb : constant C_Integer_Array(0..d-1)
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(d));
    m : constant integer32 := integer32(vb(0));
    x : Standard_Integer_Vectors.Vector(1..n);
    fltx : Standard_Floating_Vectors.Vector(1..n);
    cnt : Standard_Integer_Vectors.Vector(1..r);
    lab : Standard_Integer_Vectors.Vector(1..m);

  begin
   -- put("  r = "); put(r,1);
   -- put("  n = "); put(n,1);
   -- put("  d = "); put(natural(d),1); new_line;
    for i in cnt'range loop
      cnt(i) := integer32(vb(Interfaces.C.size_t(i)));
    end loop;
    for i in lab'range loop
      lab(i) := integer32(vb(Interfaces.C.size_t(i+r)));
    end loop;
   -- put("the number of points in each support :"); put(cnt); new_line;
   -- put("the labels of points in each support :"); put(lab); new_line;
    Assign(natural32(n),c,fltx);
    for k in 1..n loop
      x(k) := integer32(fltx(k));
    end loop;
    Integer_Cells_Container.Append_Mixed_Cell(cnt,lab,x);
    return 0;
  end Job65;

  function Job67 return integer32 is -- retrieve a mixed cell

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    cnt,lab : Standard_Integer_Vectors.Link_to_Vector;
    normal : Standard_Integer_Vectors.Link_to_Vector;
    fail : boolean;

  begin
    Integer_Cells_Container.Retrieve_Mixed_Cell(k,fail,cnt,lab,normal);
    if fail then
      return 1;
    else
      declare
        sum : constant integer32 := Standard_Integer_Vectors.Sum(cnt.all);
        cntlab : Standard_Integer_Vectors.Vector(1..cnt'last+lab'last+1);
        ind : integer32 := 0;
        fltnor : Standard_Floating_Vectors.Vector(normal'range);
      begin
        ind := ind + 1;
        cntlab(ind) := sum;
        for i in cnt'range loop
          ind := ind + 1;
          cntlab(ind) := cnt(i);
        end loop;
        for i in lab'range loop
          ind := ind + 1;
          cntlab(ind) := lab(i);
        end loop;
        Assign(cntlab,b);
        for k in normal'range loop
          fltnor(k) := double_float(normal(k));
        end loop;
        Assign(fltnor,c);
      end;
      return 0;
    end if;
  end Job67;

  function Job69 return integer32 is
  begin
    if Cells_Container.Is_Stable
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  end Job69;

  function Job70 return integer32 is

    nbr : constant natural32 := Cells_Container.Number_of_Original_Cells;

  begin
    Assign(integer32(nbr),a);
    return 0;
  end Job70;

  function Job71 return integer32 is

    nbr : constant natural32 := Cells_Container.Number_of_Stable_Cells;

  begin
    Assign(integer32(nbr),a);
    return 0;
  end Job71;

  function Cells_Standard_Stable_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Solve a stable start system in double precision.

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    mv : natural32;

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Standard_Stable_Solve ...");
    end if;
    Cells_Container.Solve_Stable_Standard_Start_System(k,mv);
    Assign(integer32(mv),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Standard_Stable_Solve.");
      end if;
      return 882;
  end Cells_Standard_Stable_Solve;

  function Cells_DoblDobl_Stable_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Solve a stable start system in double double precision.

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    mv : natural32;

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_DoblDobl_Stable_Solve ...");
    end if;
    Cells_Container.Solve_Stable_DoblDobl_Start_System(k,mv);
    Assign(integer32(mv),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_DoblDobl_Stable_Solve.");
      end if;
      return 883;
  end Cells_DoblDobl_Stable_Solve;

  function Cells_QuadDobl_Stable_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Solve a stable start system in quad double precision.

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    mv : natural32;

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_QuadDobl_Stable_Solve ...");
    end if;
    Cells_Container.Solve_Stable_QuadDobl_Start_System(k,mv);
    Assign(integer32(mv),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_QuadDobl_Stable_Solve.");
      end if;
      return 884;
  end Cells_QuadDobl_Stable_Solve;
 
  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Cells_Read_Floating_Mixed_Cells(vrblvl);
      when 1 => return Cells_Write_Floating_Mixed_Cells(vrblvl);
      when 2 => return Cells_Number_of_Floating_Mixed_Cells(vrblvl);
      when 3 => return Cells_Dimension_of_Floating_Mixed_Cells(vrblvl);
      when 4 => return Cells_Floating_Mixture(a,b,vrblvl);
      when 5 => return Job5;   -- return cardinalities of supports
      when 6 => return Job6;   -- return point in support 
      when 7 => return Job7;   -- return innner normal
      when 8 => return Job8;   -- return #elements in cell
      when 9 => return Job9;   -- return point in a cell
      when 10 => return Job10; -- return mixed volume of a cell
      when 11 => return Job11; -- set type of mixture
      when 12 => return Job12; -- append point to a support
      when 13 => return Job13; -- append mixed cell to container
      when 14 => Cells_Container.Clear; return 0;
      when 15 => return Job15; -- retrieve a mixed cell
      when 16 => return Cells_Make_Standard_Coefficient_System(vrblvl);
      when 17 => return Cells_Read_Standard_Coefficient_System(vrblvl);
      when 18 => return Cells_Write_Standard_Coefficient_System(vrblvl);
      when 19 => return Cells_Standard_System_into_Container(vrblvl);
      when 20 => return Cells_Standard_System_from_Container(vrblvl);
      when 21 => return Cells_Standard_Polyhedral_Homotopy(vrblvl);
      when 22 => return Cells_Standard_Start_Solve(a,b,vrblvl);
      when 23 => return Cells_Standard_Track_One_Path(a,b,vrblvl);
      when 24 => return Job24; -- copy target solution to st container
      when 25 => return Job25; -- permute a standard target system
      when 26 => return Cells_Make_DoblDobl_Coefficient_System(vrblvl);
      when 27 => return Cells_Read_DoblDobl_Coefficient_System(vrblvl);
      when 28 => return Cells_Write_DoblDobl_Coefficient_System(vrblvl);
      when 29 => return Cells_DoblDobl_System_into_Container(vrblvl);
      when 30 => return Cells_DoblDobl_System_from_Container(vrblvl);
      when 31 => return Cells_DoblDobl_Polyhedral_Homotopy(vrblvl);
      when 32 => return Cells_DoblDobl_Start_Solve(a,b,vrblvl);
      when 33 => return Cells_DoblDobl_Track_One_Path(a,b,vrblvl);
      when 34 => return Job34; -- copy target solution to dd container
      when 35 => return Job35; -- permute dobldobl target system
      when 36 => return Cells_Make_QuadDobl_Coefficient_System(vrblvl);
      when 37 => return Cells_Read_QuadDobl_Coefficient_System(vrblvl);
      when 38 => return Cells_Write_QuadDobl_Coefficient_System(vrblvl);
      when 39 => return Cells_QuadDobl_System_into_Container(vrblvl);
      when 40 => return Cells_QuadDobl_System_from_Container(vrblvl);
      when 41 => return Cells_QuadDobl_Polyhedral_Homotopy(vrblvl);
      when 42 => return Cells_QuadDobl_Start_Solve(a,b,vrblvl);
      when 43 => return Cells_QuadDobl_Track_One_Path(a,b,vrblvl);
      when 44 => return Job44; -- copy target solution to qd container
      when 45 => return Job45; -- permute quaddobl target system
      when 46 => return Job46; -- mixed volume computation
      when 47 => return Job47; -- initialize number of distinct supports
      when 48 => return Job48; -- copy start solution to st container
      when 49 => return Job49; -- copy start solution to dd container
      when 50 => return Job50; -- copy start solution to qd container
      when 51 => return Cells_Read_Integer_Mixed_Cells(vrblvl);
      when 52 => return Cells_Write_Integer_Mixed_Cells(vrblvl);
      when 53 => return Cells_Number_of_Integer_Mixed_Cells(vrblvl);
      when 54 => return Cells_Dimension_of_Integer_Mixed_Cells(vrblvl);
      when 55 => return Job55; -- return type of mixture for integer cells
      when 56 => return Job56; -- return cardinalities of supports
      when 57 => return Job57; -- return integer point in support 
      when 58 => return Job58; -- return integer innner normal
      when 59 => return Job59; -- return #elements in integer cell
      when 60 => return Job60; -- return point in a integer cell
      when 61 => return Job61; -- return mixed volume of a integer cell
      when 62 => return Job62; -- set number of supports
      when 63 => return Job63; -- set type of mixture
      when 64 => return Job64; -- append point to a support
      when 65 => return Job65; -- append mixed cell to container
      when 66 => Integer_Cells_Container.Clear; return 0;
      when 67 => return Job67; -- retrieve a mixed cell
      when 68 => Integer_Cells_Container.Make_Subdivision; return 0;
      when 69 => return Job69; -- returns Is_Stable
      when 70 => return Job70; -- return the number of original cells
      when 71 => return Job71; -- return the number of stable cells
      when 72 => return Cells_Standard_Stable_Solve(a,b,vrblvl);
      when 73 => return Cells_DoblDobl_Stable_Solve(a,b,vrblvl);
      when 74 => return Cells_QuadDobl_Stable_Solve(a,b,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  exception
    when others => put("Exception raised in use_celcon handling job ");
                   put(job,1); put_line(".  Will try to ignore...");
                   return 1;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_celcon;
