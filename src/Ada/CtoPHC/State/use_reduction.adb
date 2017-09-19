with Interfaces.C;
with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
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
    when others => return 1;
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
    when others => return 2;
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
    when others => return 3;
  end Job3;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 1 => return Job1;  -- scale system in the standard container
      when 2 => return Job2;  -- scale system in the dobldobl container
      when 3 => return Job3;  -- scale system in the quaddobl container
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  exception
    when others => put("Exception raised in use_scaling handling job ");
                   put(job+706,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_reduction;
