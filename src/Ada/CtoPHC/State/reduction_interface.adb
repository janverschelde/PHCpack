with Interfaces.C;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;
with Dobldobl_Complex_Poly_Systems;
with Quaddobl_Complex_Poly_Systems;
with Standard_PolySys_Container;
with Dobldobl_PolySys_Container;
with Quaddobl_PolySys_Container;
with Reduction_of_Polynomial_Systems;    use Reduction_of_Polynomial_Systems;
with Assignments_in_Ada_and_C;           use Assignments_in_Ada_and_C;

package body Reduction_Interface is

  function Reduction_Standard_Linear
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is 

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    diagonalize : constant integer32 := integer32(v_a(v_a'first));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    diagonal,inconsistent,infinite : boolean := false;

  begin
    if vrblvl > 0
     then put_line("-> in reduction_interface.Reduction_Standard_Linear ...");
    end if;
    if diagonalize = 1 
     then Sparse_Reduce(lp.all,inconsistent,infinite);
     else Reduce(lp.all,diagonal,inconsistent,infinite);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in reduction_interface.");
        put_line("Reduction_Standard_Linear.");
      end if;
      return 707;
  end Reduction_Standard_Linear;

  function Reduction_DoblDobl_Linear
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is 

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    diagonalize : constant integer32 := integer32(v_a(v_a'first));
    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    diagonal,inconsistent,infinite : boolean := false;

  begin
    if vrblvl > 0
     then put_line("-> in reduction_interface.Reduction_DoblDobl_Linear ...");
    end if;
    if diagonalize = 2 
     then Sparse_Reduce(lp.all,inconsistent,infinite);
     else Reduce(lp.all,diagonal,inconsistent,infinite);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in reduction_interface.");
        put_line("Reduction_DoblDobl_Linear.");
      end if;
      return 708;
  end Reduction_DoblDobl_Linear;

  function Reduction_QuadDobl_Linear
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is 

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    diagonalize : constant integer32 := integer32(v_a(v_a'first));
    lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    diagonal,inconsistent,infinite : boolean := false;

  begin
    if vrblvl > 0
     then put_line("-> in reduction_interface.Reduction_QuadDobl_Linear ...");
    end if;
    if diagonalize = 2 
     then Sparse_Reduce(lp.all,inconsistent,infinite);
     else Reduce(lp.all,diagonal,inconsistent,infinite);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in reduction_interface.");
        put_line("Reduction_QuadDobl_Linear.");
      end if;
      return 709;
  end Reduction_QuadDobl_Linear;

  function Reduction_Standard_Nonlinear
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    equrmax : constant natural32 := natural32(v_a(v_a'first));
    spolmax : constant natural32 := natural32(v_a(v_a'first+1));
    rpolmax : constant natural32 := natural32(v_a(v_a'first+2));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    eqcnt,spcnt,rpcnt : natural32 := 0;
    res : Standard_Complex_Poly_Systems.Poly_Sys(lp'range);
    counts : Standard_Integer_Vectors.Vector(1..3);

  begin
    if vrblvl > 0 then
      put_line("-> in reduction_interface.Reduction_Standard_Nonlinear ...");
    end if;
    Standard_Complex_Poly_Systems.Copy(lp.all,res);
    Reduce(lp.all,res,eqcnt,equrmax,spcnt,spolmax,rpcnt,rpolmax);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(res);
    counts(1) := integer32(eqcnt);
    counts(2) := integer32(spcnt);
    counts(3) := integer32(rpcnt);
    Assign(counts,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in reduction_interface.");
        put_line("Reduction_QuadDobl_Linear.");
      end if;
      return 710;
  end Reduction_Standard_Nonlinear;

end Reduction_Interface;
