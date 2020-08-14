with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Sampling_Machine;
with DoblDobl_Sampling_Machine;
with QuadDobl_Sampling_Machine;
with Standard_Monodromy_Permutations;
with DoblDobl_Monodromy_Permutations;
with QuadDobl_Monodromy_Permutations;
with Monodromy_Partitions;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with PHCpack_Operations;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Standard_Sampling_Operations;
with DoblDobl_Sampling_Operations;
with QuadDobl_Sampling_Operations;

package body Monodromy_Interface is

  function Monodromy_Standard_Initialize_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Initialize_Sampler ...");
    end if;
    Standard_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Initialize_Sampler.");
      end if;
      return 42;
  end Monodromy_Standard_Initialize_Sampler;

  function Monodromy_DoblDobl_Initialize_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Initialize_Sampler ...");
    end if;
    DoblDobl_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Initialize_Sampler.");
      end if;
      return 632;
  end Monodromy_DoblDobl_Initialize_Sampler;

  function Monodromy_QuadDobl_Initialize_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Initialize_Sampler ...");
    end if;
    QuadDobl_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Initialize_Sampler.");
      end if;
      return 662;
  end Monodromy_QuadDobl_Initialize_Sampler;

  function Monodromy_Standard_Set_Coefficient
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Numbers;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    i : constant integer32 := integer32(va(va'first));
    j : constant integer32 := integer32(vb(vb'first));
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(2));
    re : constant double_float := double_float(vc(0));
    im : constant double_float := double_float(vc(1));
    cf : constant Complex_Number := Create(re,im);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Set_Coefficient ...");
    end if;
    Standard_Sampling_Operations.Assign_Slice(cf,i,j);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Set_Coefficient.");
      end if;
      return 43;
  end Monodromy_Standard_Set_Coefficient;

  function Monodromy_DoblDobl_Set_Coefficient
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Numbers;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    i : constant integer32 := integer32(va(va'first));
    j : constant integer32 := integer32(vb(vb'first));
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(4));
    re_hi : constant double_float := double_float(vc(0));
    re_lo : constant double_float := double_float(vc(1));
    im_hi : constant double_float := double_float(vc(2));
    im_lo : constant double_float := double_float(vc(3));
    re : constant double_double := Create(re_hi,re_lo);
    im : constant double_double := Create(im_hi,im_lo);
    cf : constant Complex_Number := Create(re,im);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Set_Coefficient ...");
    end if;
    DoblDobl_Sampling_Operations.Assign_Slice(cf,i,j);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Set_Coefficient.");
      end if;
      return 633;
  end Monodromy_DoblDobl_Set_Coefficient;

  function Monodromy_QuadDobl_Set_Coefficient
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Numbers;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    i : constant integer32 := integer32(va(va'first));
    j : constant integer32 := integer32(vb(vb'first));
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(8));
    re_hihi : constant double_float := double_float(vc(0));
    re_lohi : constant double_float := double_float(vc(1));
    re_hilo : constant double_float := double_float(vc(2));
    re_lolo : constant double_float := double_float(vc(3));
    im_hihi : constant double_float := double_float(vc(4));
    im_lohi : constant double_float := double_float(vc(5));
    im_hilo : constant double_float := double_float(vc(6));
    im_lolo : constant double_float := double_float(vc(7));
    re : constant quad_double := Create(re_hihi,re_lohi,re_hilo,re_lolo);
    im : constant quad_double := Create(im_hihi,im_lohi,im_hilo,im_lolo);
    cf : constant Complex_Number := Create(re,im);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Set_Coefficient ...");
    end if;
    QuadDobl_Sampling_Operations.Assign_Slice(cf,i,j);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Set_Coefficient.");
      end if;
      return 663;
  end Monodromy_QuadDobl_Set_Coefficient;

  function Monodromy_Standard_Store_Gamma
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Numbers;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(2));
    re : constant double_float := double_float(vc(0));
    im : constant double_float := double_float(vc(1));
    gamma : constant Complex_Number := Create(re,im);
    i : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Store_Gamma ...");
    end if;
    Standard_Sampling_Operations.Store_Gamma(gamma,i);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Store_Gamma.");
      end if;
      return 44;
  end Monodromy_Standard_Store_Gamma;

  function Monodromy_DoblDobl_Store_Gamma
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Numbers;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(4));
    re_hi : constant double_float := double_float(vc(0));
    re_lo : constant double_float := double_float(vc(1));
    im_hi : constant double_float := double_float(vc(2));
    im_lo : constant double_float := double_float(vc(3));
    re : constant double_double := create(re_hi,re_lo);
    im : constant double_double := create(im_hi,im_lo);
    gamma : constant Complex_Number := Create(re,im);
    i : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Store_Gamma ...");
    end if;
    DoblDobl_Sampling_Operations.Store_Gamma(gamma,i);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Store_Gamma.");
      end if;
      return 634;
  end Monodromy_DoblDobl_Store_Gamma;

  function Monodromy_QuadDobl_Store_Gamma
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Numbers;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(8));
    re_hihi : constant double_float := double_float(vc(0));
    re_lohi : constant double_float := double_float(vc(1));
    re_hilo : constant double_float := double_float(vc(2));
    re_lolo : constant double_float := double_float(vc(3));
    im_hihi : constant double_float := double_float(vc(4));
    im_lohi : constant double_float := double_float(vc(5));
    im_hilo : constant double_float := double_float(vc(6));
    im_lolo : constant double_float := double_float(vc(7));
    re : constant quad_double := create(re_hihi,re_lohi,re_hilo,re_lolo);
    im : constant quad_double := create(im_hihi,im_lohi,im_hilo,im_lolo);
    gamma : constant Complex_Number := Create(re,im);
    i : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Store_Gamma ...");
    end if;
    QuadDobl_Sampling_Operations.Store_Gamma(gamma,i);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Store_Gamma.");
      end if;
      return 664;
  end Monodromy_QuadDobl_Store_Gamma;

  function Monodromy_Standard_Sample
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Sample ...");
    end if;
    Standard_Sampling_Operations.Sample;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Sample.");
      end if;
      return 45;
  end Monodromy_Standard_Sample;

  function Monodromy_DoblDobl_Sample
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Sample ...");
    end if;
    DoblDobl_Sampling_Operations.Sample;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Sample.");
      end if;
      return 635;
  end Monodromy_DoblDobl_Sample;

  function Monodromy_QuadDobl_Sample
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Sample ...");
    end if;
    QuadDobl_Sampling_Operations.Sample;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Sample.");
      end if;
      return 665;
  end Monodromy_QuadDobl_Sample;

  function Monodromy_Standard_Swap_Slices
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Swap_Slices ...");
    end if;
    Standard_Sampling_Operations.Swap_Slices;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Swap_Slices.");
      end if;
      return 46;
  end Monodromy_Standard_Swap_Slices;

  function Monodromy_DoblDobl_Swap_Slices
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Swap_Slices ...");
    end if;
    DoblDobl_Sampling_Operations.Swap_Slices;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Swap_Slices.");
      end if;
      return 636;
  end Monodromy_DoblDobl_Swap_Slices;

  function Monodromy_QuadDobl_Swap_Slices
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Swap_Slices ...");
    end if;
    QuadDobl_Sampling_Operations.Swap_Slices;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Swap_Slices.");
      end if;
      return 666;
  end Monodromy_QuadDobl_Swap_Slices;

  function Monodromy_Standard_Copy_System
             ( vrblvl : integer32 := 0 ) return integer32 is

     use Standard_Complex_Poly_Systems;

     p : constant Poly_Sys := Sampling_Machine.Embedded_System;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Copy_System ...");
    end if;
    Standard_PolySys_Container.Initialize(p);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Copy_System.");
      end if;
      return 47;
  end Monodromy_Standard_Copy_System;

  function Monodromy_DoblDobl_Copy_System
             ( vrblvl : integer32 := 0 ) return integer32 is

     use DoblDobl_Complex_Poly_Systems;

     p : constant Poly_Sys := DoblDobl_Sampling_Machine.Embedded_System;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Copy_System ...");
    end if;
    DoblDobl_PolySys_Container.Initialize(p);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Copy_System.");
      end if;
      return 637;
  end Monodromy_DoblDobl_Copy_System;

  function Monodromy_QuadDobl_Copy_System
             ( vrblvl : integer32 := 0 ) return integer32 is

     use QuadDobl_Complex_Poly_Systems;

     p : constant Poly_Sys := QuadDobl_Sampling_Machine.Embedded_System;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Copy_System ...");
    end if;
    QuadDobl_PolySys_Container.Initialize(p);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Copy_System.");
      end if;
      return 667;
  end Monodromy_QuadDobl_Copy_System;

  function Monodromy_Standard_Copy_Solutions
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;

    s : constant Solution_List
      := Standard_Sampling_Operations.Retrieve_First_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Copy_Solutions ...");
    end if;
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(s);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Copy_Solutions.");
      end if;
      return 48;
  end Monodromy_Standard_Copy_Solutions;

  function Monodromy_DoblDobl_Copy_Solutions
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;

    s : constant Solution_List
      := DoblDobl_Sampling_Operations.Retrieve_First_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Copy_Solutions ...");
    end if;
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(s);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Copy_Solutions.");
      end if;
      return 638;
  end Monodromy_DoblDobl_Copy_Solutions;

  function Monodromy_QuadDobl_Copy_Solutions
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;

    s : constant Solution_List
      := QuadDobl_Sampling_Operations.Retrieve_First_Solutions;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Copy_Solutions ...");
    end if;
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(s);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Copy_Solutions.");
      end if;
      return 668;
  end Monodromy_QuadDobl_Copy_Solutions;

  function Monodromy_Standard_Grid_Solutions
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));
    s : constant Solution_List
      := Standard_Monodromy_Permutations.Retrieve(i);
    cp_s : Solution_List;  -- will be a copy of s

  begin -- since this will be traffic on the same node
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Grid_Solutions ...");
    end if;
    Copy(s,cp_s); -- we better make a copy
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(cp_s);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Grid_Solutions.");
      end if;
      return 49;
  end Monodromy_Standard_Grid_Solutions;

  function Monodromy_DoblDobl_Grid_Solutions
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));
    s : constant Solution_List
      := DoblDobl_Monodromy_Permutations.Retrieve(i);
    cp_s : Solution_List;  -- will be a copy of s

  begin -- since this will be traffic on the same node
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Grid_Solutions ...");
    end if;
    Copy(s,cp_s); -- we better make a copy
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(cp_s);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Grid_Solutions.");
      end if;
      return 639;
  end Monodromy_DoblDobl_Grid_Solutions;

  function Monodromy_QuadDobl_Grid_Solutions
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));
    s : constant Solution_List
      := QuadDobl_Monodromy_Permutations.Retrieve(i);
    cp_s : Solution_List;  -- will be a copy of s

  begin -- since this will be traffic on the same node
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Grid_Solutions ...");
    end if;
    Copy(s,cp_s); -- we better make a copy
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(cp_s);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Grid_Solutions.");
      end if;
      return 669;
  end Monodromy_QuadDobl_Grid_Solutions;

  function Monodromy_Standard_Init_Permutations
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(va(va'first));
    d : constant integer32 := integer32(vb(0));
    k : constant integer32 := integer32(vb(1));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Init_Permutations ...");
    end if;
    Standard_Monodromy_Permutations.Initialize(n,d,k);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Init_Permutations.");
      end if;
      return 42;
  end Monodromy_Standard_Init_Permutations;

  function Monodromy_DoblDobl_Init_Permutations
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(va(va'first));
    d : constant integer32 := integer32(vb(0));
    k : constant integer32 := integer32(vb(1));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Init_Permutations ...");
    end if;
    DoblDobl_Monodromy_Permutations.Initialize(n,d,k);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Init_Permutations.");
      end if;
      return 640;
  end Monodromy_DoblDobl_Init_Permutations;

  function Monodromy_QuadDobl_Init_Permutations
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(va(va'first));
    d : constant integer32 := integer32(vb(0));
    k : constant integer32 := integer32(vb(1));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Init_Permutations ...");
    end if;
    QuadDobl_Monodromy_Permutations.Initialize(n,d,k);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Init_Permutations.");
      end if;
      return 670;
  end Monodromy_QuadDobl_Init_Permutations;

  function Monodromy_Standard_Perm_Solutions
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;

    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Perm_Solutions ...");
    end if;
    Standard_Monodromy_Permutations.Store(sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Perm_Solutions.");
      end if;
      return 51;
  end Monodromy_Standard_Perm_Solutions;

  function Monodromy_DoblDobl_Perm_Solutions
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;

    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Perm_Solutions ...");
    end if;
    DoblDobl_Monodromy_Permutations.Store(sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Perm_Solutions.");
      end if;
      return 641;
  end Monodromy_DoblDobl_Perm_Solutions;

  function Monodromy_QuadDobl_Perm_Solutions
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;

    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Perm_Solutions ...");
    end if;
    QuadDobl_Monodromy_Permutations.Store(sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Perm_Solutions.");
      end if;
      return 671;
  end Monodromy_QuadDobl_Perm_Solutions;

  function Monodromy_Standard_Permutation
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    perm : constant Standard_Natural_Vectors.Vector
         := Standard_Monodromy_Permutations.Permutation;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Permutation ...");
    end if;
    Assign(perm,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Permutation.");
      end if;
      return 52;
  end Monodromy_Standard_Permutation;

  function Monodromy_DoblDobl_Permutation
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    perm : constant Standard_Natural_Vectors.Vector
         := DoblDobl_Monodromy_Permutations.Permutation;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Permutation ...");
    end if;
    Assign(perm,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Permutation.");
      end if;
      return 642;
  end Monodromy_DoblDobl_Permutation;

  function Monodromy_QuadDobl_Permutation
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    perm : constant Standard_Natural_Vectors.Vector
         := QuadDobl_Monodromy_Permutations.Permutation;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Permutation ...");
    end if;
    Assign(perm,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Permutation.");
      end if;
      return 672;
  end Monodromy_QuadDobl_Permutation;

  function Monodromy_Standard_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    p : Standard_Natural_Vectors.Vector(1..n);
    nf : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Update ...");
    end if;
    Assign(natural32(n),b,p);
    Standard_Monodromy_Permutations.Update_Decomposition(p,nf(1),nf(2));
    Assign(nf,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Update.");
      end if;
      return 53;
  end Monodromy_Standard_Update;

  function Monodromy_DoblDobl_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    p : Standard_Natural_Vectors.Vector(1..n);
    nf : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Update ...");
    end if;
    Assign(natural32(n),b,p);
    DoblDobl_Monodromy_Permutations.Update_Decomposition(p,nf(1),nf(2));
    Assign(nf,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Update.");
      end if;
      return 643;
  end Monodromy_DoblDobl_Update;

  function Monodromy_QuadDobl_Update
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    p : Standard_Natural_Vectors.Vector(1..n);
    nf : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Update ...");
    end if;
    Assign(natural32(n),b,p);
    QuadDobl_Monodromy_Permutations.Update_Decomposition(p,nf(1),nf(2));
    Assign(nf,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Update.");
      end if;
      return 673;
  end Monodromy_QuadDobl_Update;

  function Monodromy_Standard_Write
             ( vrblvl : integer32 := 0 ) return integer32 is

    deco : constant Standard_Natural_VecVecs.Link_to_VecVec
         := Standard_Monodromy_Permutations.Decomposition;

    use Standard_Natural_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Write ...");
    end if;
    if deco /= null then
      if PHCpack_Operations.Is_File_Defined then
        Monodromy_Partitions.Write_Factors
          (PHCpack_Operations.output_file,deco.all);
      else
        Monodromy_Partitions.Write_Factors(standard_output,deco.all);
      end if;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Write.");
      end if;
      return 54;
  end Monodromy_Standard_Write;

  function Monodromy_DoblDobl_Write
             ( vrblvl : integer32 := 0 ) return integer32 is

    deco : constant Standard_Natural_VecVecs.Link_to_VecVec
         := DoblDobl_Monodromy_Permutations.Decomposition;

    use Standard_Natural_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Write ...");
    end if;
    if deco /= null then
      if PHCpack_Operations.Is_File_Defined then
        Monodromy_Partitions.Write_Factors
          (PHCpack_Operations.output_file,deco.all);
      else
        Monodromy_Partitions.Write_Factors(standard_output,deco.all);
      end if;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Write.");
      end if;
      return 644;
  end Monodromy_DoblDobl_Write;

  function Monodromy_QuadDobl_Write
             ( vrblvl : integer32 := 0 ) return integer32 is

    deco : constant Standard_Natural_VecVecs.Link_to_VecVec
         := QuadDobl_Monodromy_Permutations.Decomposition;

    use Standard_Natural_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Write ...");
    end if;
    if deco /= null then
      if PHCpack_Operations.Is_File_Defined then
        Monodromy_Partitions.Write_Factors
          (PHCpack_Operations.output_file,deco.all);
      else
        Monodromy_Partitions.Write_Factors(standard_output,deco.all);
      end if;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Write.");
      end if;
      return 674;
  end Monodromy_QuadDobl_Write;

  function Monodromy_Standard_Trace_Test
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    done : constant boolean
         := Standard_Monodromy_Permutations.Certify_with_Linear_Trace;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Trace_Test ...");
    end if;
    if done
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Trace_Test.");
      end if;
      return 55;
  end Monodromy_Standard_Trace_Test;

  function Monodromy_DoblDobl_Trace_Test
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    done : constant boolean
         := DoblDobl_Monodromy_Permutations.Certify_with_Linear_Trace;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Trace_Test ...");
    end if;
    if done
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Trace_Test.");
      end if;
      return 645;
  end Monodromy_DoblDobl_Trace_Test;

  function Monodromy_QuadDobl_Trace_Test
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    done : constant boolean
         := QuadDobl_Monodromy_Permutations.Certify_with_Linear_Trace;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Trace_Test ...");
    end if;
    if done
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Trace_Test.");
      end if;
      return 675;
  end Monodromy_QuadDobl_Trace_Test;

end Monodromy_Interface;
