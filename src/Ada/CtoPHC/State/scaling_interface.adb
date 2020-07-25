with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with Standard_Scaling;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;
with DoblDobl_Scaling;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_cv;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;
with QuadDobl_Scaling;
with Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Solutions;
with Multprec_PolySys_Container;
with Multprec_Solutions_Container;
with Multprec_Scaling;
with Assignments_in_Ada_and_C;           use Assignments_in_Ada_and_C;

package body Scaling_Interface is

  function Scale_Standard_System
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    stp : constant integer32 := integer32(v_a(v_a'first));
    rco : double_float;
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    scf : Standard_Complex_Vectors.Vector(1..2*lp'last+1);

  begin
    if vrblvl > 0 then
      put_line("-> in scaling_interface.Scale_Standard_System ...");
    end if;
    if stp = 0 then
      Standard_Scaling.Scale(lp.all);
    else
      if stp = 1
       then Standard_Scaling.Scale(lp.all,10,false,rco,scf(1..2*lp'last));
       else Standard_Scaling.Scale(lp.all,10,true,rco,scf(1..2*lp'last));
      end if;
      scf(scf'last) := Standard_Complex_Numbers.Create(rco);
      Assign(scf,c);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in scaling_interface.");
        put_line("Scale_Multprec_System.");
      end if;
      return 590;
  end Scale_Standard_System;

  function Scale_DoblDobl_System
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    stp : constant integer32 := integer32(v_a(v_a'first));
    rco : double_double;
    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    scf : DoblDobl_Complex_Vectors.Vector(1..2*lp'last+1);

  begin
    if vrblvl > 0 then
      put_line("-> in scaling_interface.Scale_DoblDobl_System ...");
    end if;
    if stp = 0 then
      DoblDobl_Scaling.Scale(lp.all);
    else
      if stp = 1
       then DoblDobl_Scaling.Scale(lp.all,10,false,rco,scf(1..2*lp'last));
       else DoblDobl_Scaling.Scale(lp.all,10,true,rco,scf(1..2*lp'last));
      end if;
      scf(scf'last) := DoblDobl_Complex_Numbers.Create(rco);
      Assign(scf,c);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in scaling_interface.");
        put_line("Scale_DoblDobl_System.");
      end if;
      return 591;
  end Scale_DoblDobl_System;

  function Scale_QuadDobl_System
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    stp : constant integer32 := integer32(v_a(v_a'first));
    rco : quad_double;
    lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    scf : QuadDobl_Complex_Vectors.Vector(1..2*lp'last+1);

  begin
    if vrblvl > 0 then
      put_line("-> in scaling_interface.Scale_QuadDobl_System ...");
    end if;
    if stp = 0 then
      QuadDobl_Scaling.Scale(lp.all);
    else
      if stp = 1
       then QuadDobl_Scaling.Scale(lp.all,10,false,rco,scf(1..2*lp'last));
       else QuadDobl_Scaling.Scale(lp.all,10,true,rco,scf(1..2*lp'last));
      end if;
      scf(scf'last) := QuadDobl_Complex_Numbers.Create(rco);
      Assign(scf,c);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in scaling_interface.");
        put_line("Scale_QuadDobl_System.");
      end if;
      return 592;
  end Scale_QuadDobl_System;

  function Scale_Multprec_System
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Vectors_cv;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    stp : constant integer32 := integer32(v_a(v_a'first));
    rco : Floating_Number;
    lp : constant Multprec_Complex_Poly_Systems.Link_to_Poly_Sys
       := Multprec_PolySys_Container.Retrieve;
    scf : Multprec_Complex_Vectors.Vector(1..2*lp'last+1);
    qd_scf : QuadDobl_Complex_Vectors.Vector(scf'range);

  begin
    if vrblvl > 0 then
      put_line("-> in scaling_interface.Scale_Multprec_System ...");
    end if;
    if stp = 0 then
      Multprec_Scaling.Scale(lp.all);
    else
      if stp = 1
       then Multprec_Scaling.Scale(lp.all,10,false,rco,scf(1..2*lp'last));
       else Multprec_Scaling.Scale(lp.all,10,true,rco,scf(1..2*lp'last));
      end if;
      scf(scf'last) := Multprec_Complex_Numbers.Create(rco);
      qd_scf := Multprec_to_QuadDobl_Complex(scf);
      Assign(qd_scf,c);
      Multprec_Complex_Vectors.Clear(scf);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in scaling_interface.");
        put_line("Scale_Multprec_System.");
      end if;
      return 593;
  end Scale_Multprec_System;

  function Scale_Standard_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    bas : constant natural32 := natural32(v_b(v_b'first));
    sols : Solution_List := Standard_Solutions_Container.Retrieve;
    scf : Standard_Complex_Vectors.Vector(1..dim/2); 
 
  begin
    if vrblvl > 0 then
      put_line("-> in scaling_interface.Scale_Standard_Solutions ...");
    end if;
    Assign(natural32(dim),c,scf);
    Standard_Scaling.Scale(bas,scf,sols);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in scaling_interface.");
        put_line("Scale_Standard_Solutions.");
      end if;
      return 594;
  end Scale_Standard_Solutions;

  function Scale_DoblDobl_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    bas : constant natural32 := natural32(v_b(v_b'first));
    sols : Solution_List := DoblDobl_Solutions_Container.Retrieve;
    scf : DoblDobl_Complex_Vectors.Vector(1..dim/4); 
 
  begin
    if vrblvl > 0 then
      put_line("-> in scaling_interface.Scale_DoblDobl_Solutions ...");
    end if;
    Assign(natural32(dim),c,scf);
    DoblDobl_Scaling.Scale(bas,scf,sols);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in scaling_interface.");
        put_line("Scale_DoblDobl_Solutions.");
      end if;
      return 595;
  end Scale_DoblDobl_Solutions;

  function Scale_QuadDobl_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    bas : constant natural32 := natural32(v_b(v_b'first));
    sols : Solution_List := QuadDobl_Solutions_Container.Retrieve;
    scf : QuadDobl_Complex_Vectors.Vector(1..dim/8); 
 
  begin
    if vrblvl > 0 then
      put_line("-> in scaling_interface.Scale_QuadDobl_Solutions ...");
    end if;
    Assign(natural32(dim),c,scf);
    QuadDobl_Scaling.Scale(bas,scf,sols);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in scaling_interface.");
        put_line("Scale_QuadDobl_Solutions.");
      end if;
      return 596;
  end Scale_QuadDobl_Solutions;

  function Scale_Multprec_Solutions
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Vectors_cv;
    use Multprec_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    bas : constant natural32 := natural32(v_b(v_b'first));
    sols : Solution_List := Multprec_Solutions_Container.Retrieve;
    scf : QuadDobl_Complex_Vectors.Vector(1..dim/8); 
    mp_scf : Multprec_Complex_Vectors.Vector(scf'range); 
 
  begin
    if vrblvl > 0 then
      put_line("-> in scaling_interface.Scale_Multprec_Solutions ...");
    end if;
    Assign(natural32(dim),c,scf);
    mp_scf := QuadDobl_Complex_to_Multprec(scf);
    Multprec_Scaling.Scale(bas,mp_scf,sols);
    Multprec_Complex_Vectors.Clear(mp_scf);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in scaling_interface.");
        put_line("Scale_Multprec_Solutions.");
      end if;
      return -597; -- not exported?
  end Scale_Multprec_Solutions;

end Scaling_Interface;
