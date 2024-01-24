with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
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

package body Cells_Interface is

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
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
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
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
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

  function Cells_Get_Floating_Mixture
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
  end Cells_Get_Floating_Mixture;

  function Cells_Floating_Supports_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Lists_of_Floating_Vectors;
    use Arrays_of_Floating_Vector_Lists;

    lif : Link_to_Array_of_Lists;
    r : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Floating_Supports_Size ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Floating_Supports_Size.");
      end if;
      return 85;
  end Cells_Floating_Supports_Size;

  function Cells_Get_Floating_Support_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is 

    use Arrays_of_Floating_Vector_Lists;

    lif : Link_to_Array_of_Lists;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    k_a : constant integer32 := integer32(v_a(v_a'first));
    k_b : constant natural32 := natural32(v_b(v_b'first));
    use Standard_Floating_Vectors;
    lpt : Link_to_Vector;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Get_Floating_Support_Point ...");
    end if;
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
    return 86;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Get_Floating_Support_Point.");
      end if;
      return 86;
  end Cells_Get_Floating_Support_Point;

  function Cells_Floating_Normal
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Floating_mixed_Subdivisions;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    mic : Mixed_Cell;
    fail : boolean;

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Floating_Normal ...");
    end if;
    Cells_Container.Retrieve(k,mic,fail);
    if fail
     then return 87;
     else Assign(mic.nor.all,c); return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Floating_Normal.");
      end if;
      return 87;
  end Cells_Floating_Normal;

  function Cells_Floating_Cell_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Lists_of_Floating_Vectors;
    use Floating_Mixed_Subdivisions;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    mic : Mixed_Cell;
    fail : boolean;

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Floating_Cells_Size ...");
    end if;
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
    return 88;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Floating_Cell_Size.");
      end if;
      return 88;
  end Cells_Floating_Cell_Size;

  function Cells_Get_Floating_Cell_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Get_Floating_Cell_Point ...");
    end if;
    Cells_Container.Retrieve(k,mic,fail);
    if not fail then
      lpt := Select_Point(mic.pts(i),j);
      if lpt /= null 
       then Assign(lpt.all,c); return 0;
      end if;
    end if;
    return 89;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Get_Floating_Cell_Point.");
      end if;
      return 89;
  end Cells_Get_Floating_Cell_Point;

  function Cells_Floating_Mixed_Volume
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Floating_Mixed_Volume ...");
    end if;
    Cells_Container.Retrieve(k,mic,fail);
    if fail or mix = null then
      if vrblvl > 0 then
        if fail
         then put("failed to retrieve cell "); put(k,1); new_line;
        end if;
        if mix = null
         then put_line("failed because type of mixture mix is null.");
        end if;
      end if;
      return 90;
    else
      Mixed_Volume(n,mix.all,mic,mv);
      Assign(integer32(mv),b);
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Floating_Mixed_Volume.");
      end if;
      return 90;
  end Cells_Floating_Mixed_Volume;

  function Cells_Set_Floating_Mixture
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Set_Floating_Mixture.");
      end if;
      return 240;
  end Cells_Set_Floating_Mixture;

  function Cells_Add_Floating_Support_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    n : constant integer32 := integer32(vb(vb'first));
    x : Standard_Floating_Vectors.Vector(1..n);
    fail : boolean;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Add_Floating_Support_Point ...");
    end if;
    Assign(natural32(n),c,x);
    fail := Cells_Container.Append_to_Support(k,x);
    if fail
     then return 92;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Set_Floating_Mixture.");
      end if;
      return 92;
  end Cells_Add_Floating_Support_Point;

  function Cells_Add_Floating_Mixed_Cell
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Add_Floating_Mixed_Cell ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Add_Floating_Mixed_Cell.");
      end if;
      return 93;
  end Cells_Add_Floating_Mixed_Cell;

  function Cells_Get_Floating_Mixed_Cell
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    cnt,lab : Standard_Integer_Vectors.Link_to_Vector;
    normal : Standard_Floating_Vectors.Link_to_Vector;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Get_Floating_Mixed_Cell ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Get_Floating_Mixed_Cell.");
      end if;
      return 95;
  end Cells_Get_Floating_Mixed_Cell;

  function Cells_Make_Standard_Coefficient_System
             ( vrblvl : integer32 := 0 ) return integer32 is
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

  function Cells_Standard_TarSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : Standard_Complex_Solutions.Link_to_Solution;

    use Standard_Complex_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Standard_TarSol_into_Container ...");
    end if;
    ls := Cells_Container.Retrieve_Standard_Target_Solution(k,i);
    if ls /= null
     then Standard_Solutions_Container.Append(ls.all);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Standard_TarSol_into_Container.");
      end if;
      return 104;
  end Cells_Standard_TarSol_into_Container;

  function Cells_DoblDobl_TarSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

    use DoblDobl_Complex_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_DoblDobl_TarSol_into_Container ...");
    end if;
    ls := Cells_Container.Retrieve_DoblDobl_Target_Solution(k,i);
    if ls /= null
     then DoblDobl_Solutions_Container.Append(ls.all);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_DoblDobl_TarSol_into_Container.");
      end if;
      return 468;
  end Cells_DoblDobl_TarSol_into_Container;

  function Cells_QuadDobl_TarSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

    use QuadDobl_Complex_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_QuadDobl_TarSol_into_Container ...");
    end if;
    ls := Cells_Container.Retrieve_QuadDobl_Target_Solution(k,i);
    if ls /= null
     then QuadDobl_Solutions_Container.Append(ls.all);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_QuadDobl_TarSol_into_Container.");
      end if;
      return 478;
  end Cells_QuadDobl_TarSol_into_Container;

  function Cells_Standard_StaSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Standard_StaSol_into_Container ...");
    end if;
    ls := Cells_Container.Retrieve_Standard_Start_Solution(k,i);
    Standard_Solutions_Container.Append(ls.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_DoblDobl_StaSol_into_Container.");
      end if;
      return 597;
  end Cells_Standard_StaSol_into_Container;

  function Cells_DoblDobl_StaSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_DoblDobl_StaSol_into_Container ...");
    end if;
    ls := Cells_Container.Retrieve_DoblDobl_Start_Solution(k,i);
    DoblDobl_Solutions_Container.Append(ls.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_DoblDobl_StaSol_into_Container.");
      end if;
      return 598;
  end Cells_DoblDobl_StaSol_into_Container;

  function Cells_QuadDobl_StaSol_into_Container
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    i : constant natural32 := natural32(vb(vb'first));
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_QuadDobl_StaSol_into_Container ...");
    end if;
    ls := Cells_Container.Retrieve_QuadDobl_Start_Solution(k,i);
    QuadDobl_Solutions_Container.Append(ls.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_QuadDobl_StaSol_into_Container.");
      end if;
      return 599;
  end Cells_QuadDobl_StaSol_into_Container;

  function Cells_Standard_Permute
             ( vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_QuadDobl_Permute ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Standard_Permute.");
      end if;
      return 105;
  end Cells_Standard_Permute;

  function Cells_DoblDobl_Permute
             ( vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_DoblDobl_Permute ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_DoblDobl_Permute.");
      end if;
      return 469;
  end Cells_DoblDobl_Permute;

  function Cells_QuadDobl_Permute
             ( vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_QuadDobl_Permute ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_QuadDobl_Permute.");
      end if;
      return 479;
  end Cells_QuadDobl_Permute;

  function Cells_Floating_Mixed_Volume
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    mv : constant natural32 := Cells_Container.Mixed_Volume;

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Floating_Mixed_Volume ...");
    end if;
    Assign(integer32(mv),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Floating_Mixed_Volume.");
      end if;
      return 239;
  end Cells_Floating_Mixed_Volume;

  function Cells_Set_Floating_Number_of_Supports
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbr : constant natural32 := natural32(v_a(v_a'first));

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Set_Floating_Number_of_Supports ...");
    end if;
    Cells_Container.Initialize_Supports(nbr);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Set_Floating_Number_of_Supports.");
      end if;
      return 240;
  end Cells_Set_Floating_Number_of_Supports;

  function Cells_Read_Integer_Mixed_Cells
             ( vrblvl : integer32 := 0 ) return integer32 is

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
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
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
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
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

  function Cells_Get_Integer_Mixture
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Integer_Vectors;
    mix : Link_to_Vector;
    r : integer32;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Integer_Mixture ...");
    end if;
    mix := Integer_Cells_Container.Type_of_Mixture;
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
        put_line("Cells_Integer_Mixture.");
      end if;
      return 745;
  end Cells_Get_Integer_Mixture;

  function Cells_Integer_Supports_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Lists_of_Integer_Vectors;
    use Arrays_of_Integer_Vector_Lists;

    lif : Link_to_Array_of_Lists;
    r : integer32;

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Integer_Supports_Size ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Integer_Supports_Size.");
      end if;
      return 746;
  end Cells_Integer_Supports_Size;

  function Cells_Get_Integer_Support_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is 

    use Arrays_of_Integer_Vector_Lists;

    lif : Link_to_Array_of_Lists;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    k_a : constant integer32 := integer32(v_a(v_a'first));
    k_b : constant natural32 := natural32(v_b(v_b'first));
    use Standard_Integer_Vectors;
    lpt : Link_to_Vector;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Get_Integer_Support_Point ...");
    end if;
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
    return 747;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Integer_Supports_Size.");
      end if;
      return 747;
  end Cells_Get_Integer_Support_Point;

  function Cells_Integer_Normal
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Integer_Mixed_Subdivisions;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    mic : Mixed_Cell;
    fail : boolean;

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Integer_Normal ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Integer_Normal.");
      end if;
      return 748;
  end Cells_Integer_Normal;

  function Cells_Integer_Cell_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Lists_of_Integer_Vectors;
    use Integer_Mixed_Subdivisions;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v(v'first));
    mic : Mixed_Cell;
    fail : boolean;

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Integer_Cells_Size ...");
    end if;
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
    return 749;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Integer_Cell_Size.");
      end if;
      return 749;
  end Cells_Integer_Cell_Size;

  function Cells_Get_Integer_Cell_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Get_Integer_Cell_Point ...");
    end if;
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
    return 750;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Get_Integer_Cell_Point.");
      end if;
      return 750;
  end Cells_Get_Integer_Cell_Point;

  function Cells_Integer_Mixed_Volume
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Integer_Mixed_Volume ...");
    end if;
   -- put("The dimension : "); put(n,1); new_line;
   -- put("retrieving mixed cell "); put(k,1); put_line(" ...");
    Integer_Cells_Container.Retrieve(k,mic,fail);
    if fail or mix = null then
     -- put_line("failure!");
      return 751;
    else
     -- put("The type of mixture :"); put(mix); new_line;
      Mixed_Volume(n,mix.all,mic,mv);
     -- put("The mixed volume mv is "); put(mv,1); put_line(".");
      Assign(integer32(mv),b);
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Integer_Mixed_Volume.");
      end if;
      return 751;
  end Cells_Integer_Mixed_Volume;

  function Cells_Set_Integer_Number_of_Supports
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Sets the number of different supports,
  --   lifted with integer values.

  -- ON ENTRY :
  --   a       numbers of different supports;
  --   vrblvl  the verbose level.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbr : constant natural32 := natural32(v_a(v_a'first));

  begin
    if vrblvl > 0 then
      put("-> in cells_interface.");
      put_line("Cells_Set_Integer_Number_of_Supports ...");
    end if;
    Integer_Cells_Container.Initialize_Supports(nbr);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Set_Integer_Number_of_Supports.");
      end if;
      return 752;
  end Cells_Set_Integer_Number_of_Supports;

  function Cells_Set_Integer_Mixture
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    r : constant integer32 := integer32(v(v'first));
    mix : Standard_Integer_Vectors.Vector(1..r);
    lnkmix : Standard_Integer_Vectors.Link_to_Vector;
    r1 : constant Interfaces.C.size_t := Interfaces.C.size_t(r-1);
    mbv : constant C_Integer_Array(0..r1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(r));

  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Set_Integer_Mixture ...");
    end if;
    for i in 0..r-1 loop
      mix(i+1) := integer32(mbv(Interfaces.C.size_t(i)));
    end loop;
    lnkmix := new Standard_Integer_Vectors.Vector'(mix);
    Integer_Cells_Container.Initialize(lnkmix);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Set_Integer_Mixture.");
      end if;
      return 753;
  end Cells_Set_Integer_Mixture;

  function Cells_Add_Integer_Support_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    k : constant natural32 := natural32(va(va'first));
    n : constant integer32 := integer32(vb(vb'first));
    x : Standard_Integer_Vectors.Vector(1..n);
    fltx : Standard_Floating_Vectors.Vector(1..n);
    fail : boolean;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Add_Integer_Support_Point ...");
    end if;
    Assign(natural32(n),c,fltx);
    for k in 1..n loop
      x(k) := integer32(fltx(k));
    end loop;
    fail := Integer_Cells_Container.Append_to_Support(k,x);
    if fail
     then return 754;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Add_Integer_Support_Point.");
      end if;
      return 754;
  end Cells_Add_Integer_Support_Point;

  function Cells_Make_Integer_Subdivision
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Make_Integer_Subdivision ...");
    end if;
    Integer_Cells_Container.Make_Subdivision;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Make_Integer_Subdivision.");
      end if;
      return 758;
  end Cells_Make_Integer_Subdivision;

  function Cells_Add_Integer_Mixed_Cell
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Add_Integer_Mixed_Cell ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Add_Integer_Mixed_Cell.");
      end if;
      return 755;
  end Cells_Add_Integer_Mixed_Cell;

  function Cells_Get_Integer_Mixed_Cell
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(va(va'first));
    cnt,lab : Standard_Integer_Vectors.Link_to_Vector;
    normal : Standard_Integer_Vectors.Link_to_Vector;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Get_Integer_Mixed_Cell ...");
    end if;
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
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Get_Integer_Mixed_Cell.");
      end if;
      return 757;
  end Cells_Get_Integer_Mixed_Cell;

  function Cells_Is_Stable
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Is_Stable ...");
    end if;
    if Cells_Container.Is_Stable
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Is_Stable.");
      end if;
      return 879;
  end Cells_Is_Stable;

  function Cells_Number_of_Original_Cells
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    nbr : constant natural32 := Cells_Container.Number_of_Original_Cells;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Number_of_Original_Cells ...");
    end if;
    Assign(integer32(nbr),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Number_of_Original_Cells.");
      end if;
      return 880;
  end Cells_Number_of_Original_Cells;

  function Cells_Number_of_Stable_Cells
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    nbr : constant natural32 := Cells_Container.Number_of_Stable_Cells;

  begin
    if vrblvl > 0 then
      put_line("-> in cells_interface.Cells_Number_of_Stable_Cells ...");
    end if;
    Assign(integer32(nbr),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Number_of_Stable_Cells.");
      end if;
      return 881;
  end Cells_Number_of_Stable_Cells;

  function Cells_Standard_Stable_Solve
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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

  function Cells_Floating_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Floating_Clear ...");
    end if;
    Cells_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Floating_Clear.");
      end if;
      return 94;
  end Cells_Floating_Clear;

  function Cells_Integer_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in cells_interface.Cells_Integer_Clear ...");
    end if;
    Integer_Cells_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cells_interface.");
        put_line("Cells_Integer_Clear.");
      end if;
      return 756;
  end Cells_Integer_Clear;
 
end Cells_Interface;
