with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Sampling_Machine;
with Sampling_Laurent_Machine;
with DoblDobl_Sampling_Machine;
with DoblDobl_Sampling_Laurent_Machine;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sampling_Laurent_Machine;
with Standard_Monodromy_Permutations;
with DoblDobl_Monodromy_Permutations;
with QuadDobl_Monodromy_Permutations;
with Monodromy_Partitions;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with PHCpack_Operations;
with Standard_PolySys_Container;
with Standard_LaurSys_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_LaurSys_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_LaurSys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Standard_Sampling_Operations;
with DoblDobl_Sampling_Operations;
with QuadDobl_Sampling_Operations;

with Standard_Integer_Numbers_io;
 use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors_io;
 use Standard_Natural_Vectors_io;

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
         := Standard_Monodromy_Permutations.Permutation(vrblvl);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Permutation ...");
      put("vrblvl : "); put(vrblvl,1); new_line;
      put("the permutation :"); put(perm); new_line;
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

  function Monodromy_Standard_Diagnostics
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Numbers;

    err,dis : double_float;
    ada_c : Complex_Number;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Diagnostics ...");
    end if;
    Standard_Monodromy_Permutations.Trace_Grid_Diagnostics(err,dis);
    ada_c := Create(err,dis);  -- a complex number is an array
    Assign(ada_c,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Diagnostics.");
      end if;
      return 56;
  end Monodromy_Standard_Diagnostics;

  function Monodromy_DoblDobl_Diagnostics
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Numbers; -- return hi_parts of double doubles

    dd_err,dd_dis : double_double;
    st_err,st_dis : double_float;
    ada_c : Complex_Number;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Diagnostics ...");
    end if;
    DoblDobl_Monodromy_Permutations.Trace_Grid_Diagnostics(dd_err,dd_dis);
    st_err := to_double(dd_err);
    st_dis := to_double(dd_dis);
    ada_c := Create(st_err,st_dis);  -- a complex number is an array
    Assign(ada_c,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Diagnostics.");
      end if;
      return 646;
  end Monodromy_DoblDobl_Diagnostics;

  function Monodromy_QuadDobl_Diagnostics
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Numbers; -- return hi_parts of quad doubles

    qd_err,qd_dis : quad_double;
    st_err,st_dis : double_float;
    ada_c : Complex_Number;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Diagnostics ...");
    end if;
    QuadDobl_Monodromy_Permutations.Trace_Grid_Diagnostics(qd_err,qd_dis);
    st_err := to_double(qd_err);
    st_dis := to_double(qd_dis);
    ada_c := Create(st_err,st_dis);  -- a complex number is an array
    Assign(ada_c,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Diagnostics.");
      end if;
      return 676;
  end Monodromy_QuadDobl_Diagnostics;

  function Monodromy_Standard_Trace_Sum
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    f : Standard_Natural_Vectors.Vector(1..n);
    d : double_float;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Trace_Sum ...");
    end if;
    Assign(natural32(n),b,f);
    d := Standard_Monodromy_Permutations.Trace_Sum_Difference(f);
    Assign(d,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Trace_Sum.");
      end if;
      return 57;
  end Monodromy_Standard_Trace_Sum;

  function Monodromy_DoblDobl_Trace_Sum
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    f : Standard_Natural_Vectors.Vector(1..n);
    dd_d : double_double;
    st_d : double_float;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Trace_Sum ...");
    end if;
    Assign(natural32(n),b,f);
    dd_d := DoblDobl_Monodromy_Permutations.Trace_Sum_Difference(f);
    st_d := to_double(dd_d);
    Assign(st_d,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Trace_Sum.");
      end if;
      return 647;
  end Monodromy_DoblDobl_Trace_Sum;

  function Monodromy_QuadDobl_Trace_Sum
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    f : Standard_Natural_Vectors.Vector(1..n);
    qd_d : quad_double;
    st_d : double_float;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Trace_Sum ...");
    end if;
    Assign(natural32(n),b,f);
    qd_d := QuadDobl_Monodromy_Permutations.Trace_Sum_Difference(f);
    st_d := to_double(qd_d);
    Assign(st_d,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Trace_Sum.");
      end if;
      return 677;
  end Monodromy_QuadDobl_Trace_Sum;

  function Monodromy_Standard_Initialize_Slices
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    nb : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Initialize_Slices ...");
    end if;
    Standard_Sampling_Operations.Initialize_Slices(nb);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Initialize_Slices.");
      end if;
      return 59;
  end Monodromy_Standard_Initialize_Slices;

  function Monodromy_DoblDobl_Initialize_Slices
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    nb : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Initialize_Slices ...");
    end if;
    DoblDobl_Sampling_Operations.Initialize_Slices(nb);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Initialize_Slices.");
      end if;
      return 649;
  end Monodromy_DoblDobl_Initialize_Slices;

  function Monodromy_QuadDobl_Initialize_Slices
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    nb : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Initialize_Slices ...");
    end if;
    QuadDobl_Sampling_Operations.Initialize_Slices(nb);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Initialize_Slices.");
      end if;
      return 679;
  end Monodromy_QuadDobl_Initialize_Slices;

  function Monodromy_Standard_Index
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    label : constant integer32 := integer32(va(0));
    slice : constant integer32 := integer32(va(1));
    result : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in monodromy_interface.Monodromy_Standard_Index ...");
    end if;
    result := Standard_Monodromy_Permutations.In_Slice(label,slice);
    Assign(result,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Index.");
      end if;
      return 58;
  end Monodromy_Standard_Index;

  function Monodromy_DoblDobl_Index
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    label : constant integer32 := integer32(va(0));
    slice : constant integer32 := integer32(va(1));
    result : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in monodromy_interface.Monodromy_DoblDobl_Index ...");
    end if;
    result := DoblDobl_Monodromy_Permutations.In_Slice(label,slice);
    Assign(result,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Index.");
      end if;
      return 648;
  end Monodromy_DoblDobl_Index;

  function Monodromy_QuadDobl_Index
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    label : constant integer32 := integer32(va(0));
    slice : constant integer32 := integer32(va(1));
    result : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in monodromy_interface.Monodromy_QuadDobl_Index ...");
    end if;
    result := QuadDobl_Monodromy_Permutations.In_Slice(label,slice);
    Assign(result,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Index.");
      end if;
      return 678;
  end Monodromy_QuadDobl_index;

  function Monodromy_Standard_Set_Target
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Set_Target ...");
    end if;
    Standard_Sampling_Operations.Set_Target_Slices(i);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Set_Target.");
      end if;
      return 62;
  end Monodromy_Standard_Set_Target;

  function Monodromy_DoblDobl_Set_Target
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Set_Target ...");
    end if;
    DoblDobl_Sampling_Operations.Set_Target_Slices(i);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Set_Target.");
      end if;
      return 652;
  end Monodromy_DoblDobl_Set_Target;

  function Monodromy_QuadDobl_Set_Target
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Set_Target ...");
    end if;
    QuadDobl_Sampling_Operations.Set_Target_Slices(i);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Set_Target.");
      end if;
      return 682;
  end Monodromy_QuadDobl_Set_Target;

  function Monodromy_Standard_Loop
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    start_slice : constant integer32 := integer32(va(0));
    target_slice : constant integer32 := integer32(va(1));
    start_label : constant integer32 := integer32(vb(vb'first));
    target_label : integer32;
    tol : constant double_float := 1.0E-8;
    sls,tls : Link_to_Solution;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Loop ...");
    end if;
    if start_slice = 0 then
      sls := Standard_Monodromy_Permutations.Retrieve
               (start_label,start_slice);
    else
      sls := Standard_Monodromy_Permutations.Retrieve
               (start_label,start_slice+2);
                         -- +2 for trace grid
    end if;
    tls := Standard_Sampling_Operations.Sample_Loop
             (start_slice,target_slice,sls);
    if target_slice = 0 then
      target_label := Standard_Monodromy_Permutations.Match
                        (tls,target_slice,tol);
    else
      target_label := Standard_Monodromy_Permutations.Match
                        (tls,target_slice+2,tol);
    end if;
    Assign(target_label,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Loop.");
      end if;
      return 63;
  end Monodromy_Standard_Loop;

  function Monodromy_DoblDobl_Loop
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    start_slice : constant integer32 := integer32(va(0));
    target_slice : constant integer32 := integer32(va(1));
    start_label : constant integer32 := integer32(vb(vb'first));
    target_label : integer32;
    tol : constant double_float := 1.0E-8;
    sls,tls : Link_to_Solution;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Loop ...");
    end if;
    if start_slice = 0 then
      sls := DoblDobl_Monodromy_Permutations.Retrieve
               (start_label,start_slice);
    else
      sls := DoblDobl_Monodromy_Permutations.Retrieve
               (start_label,start_slice+2);
                         -- +2 for trace grid
    end if;
    tls := DoblDobl_Sampling_Operations.Sample_Loop
             (start_slice,target_slice,sls);
    if target_slice = 0 then
      target_label := DoblDobl_Monodromy_Permutations.Match
                        (tls,target_slice,tol);
    else
      target_label := DoblDobl_Monodromy_Permutations.Match
                        (tls,target_slice+2,tol);
    end if;
    Assign(target_label,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Loop.");
      end if;
      return 653;
  end Monodromy_DoblDobl_Loop;

  function Monodromy_QuadDobl_Loop
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    start_slice : constant integer32 := integer32(va(0));
    target_slice : constant integer32 := integer32(va(1));
    start_label : constant integer32 := integer32(vb(vb'first));
    target_label : integer32;
    tol : constant double_float := 1.0E-8;
    sls,tls : Link_to_Solution;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Loop ...");
    end if;
    if start_slice = 0 then
      sls := QuadDobl_Monodromy_Permutations.Retrieve
               (start_label,start_slice);
    else
      sls := QuadDobl_Monodromy_Permutations.Retrieve
               (start_label,start_slice+2);
                         -- +2 for trace grid
    end if;
    tls := QuadDobl_Sampling_Operations.Sample_Loop
             (start_slice,target_slice,sls);
    if target_slice = 0 then
      target_label := QuadDobl_Monodromy_Permutations.Match
                        (tls,target_slice,tol);
    else
      target_label := QuadDobl_Monodromy_Permutations.Match
                        (tls,target_slice+2,tol);
    end if;
    Assign(target_label,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Loop.");
      end if;
      return 683;
  end Monodromy_QuadDobl_Loop;

  function Monodromy_Standard_Factor_Count
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    f : constant natural32
      := Standard_Monodromy_Permutations.Number_of_Irreducible_Factors;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Factor_Count ...");
    end if;
    Assign(integer32(f),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Factor_Count.");
      end if;
      return 68;
  end Monodromy_Standard_Factor_Count;

  function Monodromy_DoblDobl_Factor_Count
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    f : constant natural32
      := DoblDobl_Monodromy_Permutations.Number_of_Irreducible_Factors;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Factor_Count ...");
    end if;
    Assign(integer32(f),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Factor_Count.");
      end if;
      return 656;
  end Monodromy_DoblDobl_Factor_Count;

  function Monodromy_QuadDobl_Factor_Count
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    f : constant natural32
      := QuadDobl_Monodromy_Permutations.Number_of_Irreducible_Factors;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Factor_Count ...");
    end if;
    Assign(integer32(f),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Factor_Count.");
      end if;
      return 686;
  end Monodromy_QuadDobl_Factor_Count;

  function Monodromy_Standard_Get_Factor
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant integer32 := integer32(v_a(v_a'first));
    f : constant Standard_Natural_Vectors.Link_to_Vector
      := Standard_Monodromy_Permutations.Component(k);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Get_Factor ...");
    end if;
    Assign(f'last,a);
    Assign(f.all,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Get_Factor.");
      end if;
      return 69;
  end Monodromy_Standard_Get_Factor;

  function Monodromy_DoblDobl_Get_Factor
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant integer32 := integer32(v_a(v_a'first));
    f : constant Standard_Natural_Vectors.Link_to_Vector
      := DoblDobl_Monodromy_Permutations.Component(k);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Get_Factor ...");
    end if;
    Assign(f'last,a);
    Assign(f.all,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Get_Factor.");
      end if;
      return 657;
  end Monodromy_DoblDobl_Get_Factor;

  function Monodromy_QuadDobl_Get_Factor
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant integer32 := integer32(v_a(v_a'first));
    f : constant Standard_Natural_Vectors.Link_to_Vector
      := QuadDobl_Monodromy_Permutations.Component(k);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Get_Factor ...");
    end if;
    Assign(f'last,a);
    Assign(f.all,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Get_Factor.");
      end if;
      return 687;
  end Monodromy_QuadDobl_Get_Factor;

  function Monodromy_Standard_Set_Silent
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Set_Silent ...");
    end if;
    Standard_Monodromy_Permutations.stay_silent := true;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Set_Silent.");
      end if;
      return 39;
  end Monodromy_Standard_Set_Silent;

  function Monodromy_DoblDobl_Set_Silent
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Set_Silent ...");
    end if;
    DoblDobl_Monodromy_Permutations.stay_silent := true;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Set_Silent.");
      end if;
      return 658;
  end Monodromy_DoblDobl_Set_Silent;

  function Monodromy_QuadDobl_Set_Silent
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Set_Silent ...");
    end if;
    QuadDobl_Monodromy_Permutations.stay_silent := true;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Set_Silent.");
      end if;
      return 688;
  end Monodromy_QuadDobl_Set_Silent;

  function Monodromy_Standard_Set_Verbose
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Set_Verbose ...");
    end if;
    Standard_Monodromy_Permutations.stay_silent := false;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Set_Verbose.");
      end if;
      return 630;
  end Monodromy_Standard_Set_Verbose;

  function Monodromy_DoblDobl_Set_Verbose
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Set_Verbose ...");
    end if;
    DoblDobl_Monodromy_Permutations.stay_silent := false;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Set_Verbose.");
      end if;
      return 660;
  end Monodromy_DoblDobl_Set_Verbose;

  function Monodromy_QuadDobl_Set_Verbose
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Set_Verbose ...");
    end if;
    QuadDobl_Monodromy_Permutations.stay_silent := false;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Set_Verbose.");
      end if;
      return 690;
  end Monodromy_QuadDobl_Set_Verbose;

  function Monodromy_Standard_Random
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Numbers;

    res : constant Complex_Number := Standard_Random_Numbers.Random1;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Random ...");
    end if;
    Assign(res,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Random.");
      end if;
      return 280;
  end Monodromy_Standard_Random;

  function Monodromy_DoblDobl_Random
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Numbers;

    res : constant Complex_Number := DoblDobl_Random_Numbers.Random1;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Random ...");
    end if;
    Assign(res,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Random.");
      end if;
      return 659;
  end Monodromy_DoblDobl_Random;

  function Monodromy_QuadDobl_Random
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Numbers;

    res : constant Complex_Number := QuadDobl_Random_Numbers.Random1;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Random ...");
    end if;
    Assign(res,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Random.");
      end if;
      return 689;
  end Monodromy_QuadDobl_Random;

  function Convert_to_Hyperplanes
             ( v : Standard_Floating_Vectors.Vector; k,n : integer32 )
             return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Uses the numbers in v to create k hyperplanes in n-space.

    res : Standard_Complex_VecVecs.VecVec(1..k);
    ind : integer32 := v'first;

  begin
    for i in 1..k loop 
      declare
        hyp : Standard_Complex_Vectors.Vector(0..n);
      begin
        for j in 0..n loop
          hyp(j) := Standard_Complex_Numbers.Create(v(ind),v(ind+1));
          ind := ind+2;
        end loop;
        res(i) := new Standard_Complex_Vectors.Vector'(hyp);
      end;
    end loop;
    return res;
  end Convert_to_Hyperplanes;

  function Convert_to_Coefficients
             ( n : integer32; v : Standard_Complex_VecVecs.VecVec )
             return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 1..n with the complex coefficients of v
  --   stored as sequences of real and imaginary parts.

    res : Standard_Floating_Vectors.Vector(1..n);
    lv : Standard_Complex_Vectors.Link_to_Vector;
    ind : integer32 := 0;

  begin
    for i in v'range loop
      lv := v(i);
      for j in lv'range loop
        ind := ind + 1; res(ind) := Standard_Complex_Numbers.REAL_PART(lv(j));
        ind := ind + 1; res(ind) := Standard_Complex_Numbers.IMAG_PART(lv(j));
      end loop;
    end loop;
    return res;
  end Convert_to_Coefficients;

  function Monodromy_Standard_Add_Slice
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    nb_cff : constant integer32 := integer32(va(0));
    k : constant integer32 := integer32(va(1));
    n : constant integer32 := integer32(va(2));
    cff : Standard_Floating_Vectors.Vector(1..nb_cff);
    v : Standard_Complex_VecVecs.VecVec(1..k);

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Add_Slice ...");
    end if;
    Assign(natural32(nb_cff),c,cff);
    v := Convert_to_Hyperplanes(cff,k,n);
    Standard_Sampling_Operations.Add_Slices(v);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Add_Slice.");
      end if;
      return 60;
  end Monodromy_Standard_Add_Slice;

  function Monodromy_Standard_Get_Slice
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    vb : constant C_Integer_Array := C_intarrs.Value(b);
    i : constant integer32 := integer32(vb(vb'first));
    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    nb_cff : constant integer32 := integer32(va(0));
   -- dim : constant natural32 := natural32(va(1));
   -- n : constant natural32 := natural32(va(2));
    v : Standard_Complex_VecVecs.Link_to_VecVec;
    cff : Standard_Floating_Vectors.Vector(1..nb_cff);
    use Standard_Complex_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Get_Slice ...");
    end if;
    v := Standard_Sampling_Operations.Retrieve_Slices(i);
    if v /= null then
      cff := Convert_to_Coefficients(nb_cff,v.all);
      Assign(cff,c);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Get_Slice.");
      end if;
      return 61;
  end Monodromy_Standard_Get_Slice;

  function Monodromy_Standard_Init_Laurent_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Init_Laurent_Sampler ...");
    end if;
    Standard_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Init_Sampler.");
      end if;
      return 804;
  end Monodromy_Standard_Init_Laurent_sampler;

  function Monodromy_DoblDobl_Init_Laurent_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Init_Laurent_Sampler ...");
    end if;
    DoblDobl_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Init_Sampler.");
      end if;
      return 805;
  end Monodromy_DoblDobl_Init_Laurent_Sampler;

  function Monodromy_QuadDobl_Init_Laurent_Sampler
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Init_Laurent_Sampler ...");
    end if;
    QuadDobl_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Init_Sampler.");
      end if;
      return 806;
  end Monodromy_QuadDobl_Init_Laurent_Sampler;

  function Monodromy_Standard_Copy_Laurent_System
             ( vrblvl : integer32 := 0 ) return integer32 is

     use Standard_Complex_Laur_Systems;

     p : constant Laur_Sys := Sampling_Laurent_Machine.Embedded_System;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Copy_Laurent_System ...");
    end if;
    Standard_LaurSys_Container.Initialize(p);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Copy_System.");
      end if;
      return 807;
  end Monodromy_Standard_Copy_Laurent_System;

  function Monodromy_DoblDobl_Copy_Laurent_System
             ( vrblvl : integer32 := 0 ) return integer32 is

     use DoblDobl_Complex_Laur_Systems;

     p : constant Laur_Sys
       := DoblDobl_Sampling_Laurent_Machine.Embedded_System;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Copy_Laurent_System ...");
    end if;
    DoblDobl_LaurSys_Container.Initialize(p);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Copy_System.");
      end if;
      return 808;
  end Monodromy_DoblDobl_Copy_Laurent_System;

  function Monodromy_QuadDobl_Copy_Laurent_System
             ( vrblvl : integer32 := 0 ) return integer32 is

     use QuadDobl_Complex_Laur_Systems;

     p : constant Laur_Sys
       := QuadDobl_Sampling_Laurent_Machine.Embedded_System;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Copy_Laurent_System ...");
    end if;
    QuadDobl_LaurSys_Container.Initialize(p);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Copy_System.");
      end if;
      return 808;
  end Monodromy_QuadDobl_Copy_Laurent_System;

end Monodromy_Interface;
