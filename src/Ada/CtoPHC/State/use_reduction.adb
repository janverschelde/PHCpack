with Interfaces.C;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;
with Dobldobl_Complex_Poly_Systems;
with Quaddobl_Complex_Poly_Systems;
with Standard_PolySys_Container;
with Dobldobl_PolySys_Container;
with Quaddobl_PolySys_Container;
with Reduction_of_Polynomial_Systems;    use Reduction_of_Polynomial_Systems;
with Assignments_in_Ada_and_C;           use Assignments_in_Ada_and_C;

function use_reduction ( job : integer32;
                         a : C_intarrs.Pointer;
                         b : C_intarrs.Pointer;
                         c : C_dblarrs.Pointer ) return integer32 is

  function Job1 return integer32 is -- reduce standard system

  -- DESCRIPTION :
  --   Extracts from the parameter a the diagonalize flag,
  --   retrieves the system from the container,
  --   and applies linear reduction in standard double precision.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    diagonalize : constant integer32 := integer32(v_a(v_a'first));
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    diagonal,inconsistent,infinite : boolean := false;

  begin
    if diagonalize = 1 
     then Sparse_Reduce(lp.all,inconsistent,infinite);
     else Reduce(lp.all,diagonal,inconsistent,infinite);
    end if;
    return 0;
  exception
    when others => put_line("Exception raised in standard linear reduction.");
                   raise;
  end Job1;

  function Job2 return integer32 is -- reduce dobldobl system

  -- DESCRIPTION :
  --   Extracts from the parameter a the diagonalize flag,
  --   retrieves the system from the container,
  --   and applies linear reduction in double double precision.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    diagonalize : constant integer32 := integer32(v_a(v_a'first));
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    diagonal,inconsistent,infinite : boolean := false;

  begin
    if diagonalize = 2 
     then Sparse_Reduce(lp.all,inconsistent,infinite);
     else Reduce(lp.all,diagonal,inconsistent,infinite);
    end if;
    return 0;
  exception
    when others => put_line("Exception raised in dobldobl linear reduction.");
                   raise;
  end Job2;

  function Job3 return integer32 is -- reduce quaddobl system

  -- DESCRIPTION :
  --   Extracts from the parameter a the diagonalize flag,
  --   retrieves the system from the container,
  --   and applies linear reduction in quad double precision.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    diagonalize : constant integer32 := integer32(v_a(v_a'first));
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    diagonal,inconsistent,infinite : boolean := false;

  begin
    if diagonalize = 2 
     then Sparse_Reduce(lp.all,inconsistent,infinite);
     else Reduce(lp.all,diagonal,inconsistent,infinite);
    end if;
    return 0;
  exception
    when others => put_line("Exception raised in quaddobl linear reduction.");
                   raise;
  end Job3;

  function Job4 return integer32 is -- standard nonlinear reduction

  -- DESCRIPTION :
  --   Extracts from the parameter a the diagonalize flag,
  --   retrieves the system from the container,
  --   and applies linear reduction in standard double precision.

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    equrmax : constant natural32 := natural32(v_a(v_a'first));
    spolmax : constant natural32 := natural32(v_a(v_a'first+1));
    rpolmax : constant natural32 := natural32(v_a(v_a'first+2));
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    eqcnt,spcnt,rpcnt : natural32 := 0;
    res : Standard_Complex_Poly_Systems.Poly_Sys(lp'range);
    counts : Standard_Integer_Vectors.Vector(1..3);

  begin
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
      put_line("Exception raised in standard nonlinear reduction.");
      raise;
  end Job4;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 1 => return Job1;  -- reduce system in the standard container
      when 2 => return Job2;  -- reduce system in the dobldobl container
      when 3 => return Job3;  -- reduce system in the quaddobl container
      when 4 => return Job4;  -- standard nonlinear reduction
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  exception
    when others => put("Exception raised in use_reduction handling job ");
                   put(job+706,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_reduction;
