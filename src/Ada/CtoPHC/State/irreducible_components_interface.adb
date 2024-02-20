with Interfaces.C;
with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Natural_VecVecs;
with Standard_Integer_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;    use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;    use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;    use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;    use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;    use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;    use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Standard_LaurSys_Container;
with DoblDobl_LaurSys_Container;
with QuadDobl_LaurSys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Embeddings_and_Cascades;             use Embeddings_and_Cascades;
with Standard_Witness_Solutions;
with DoblDobl_Witness_Solutions;
with QuadDobl_Witness_Solutions;
with Store_Witness_Solutions;             use Store_Witness_Solutions;
with Write_Witness_Solutions;             use Write_Witness_Solutions;
with Path_Counts_Table;
with Assignments_in_Ada_and_C;            use Assignments_in_Ada_and_C;

package body Irreducible_Components_Interface is

  procedure extract_solver_options
              ( a : in C_intarrs.Pointer; b : in C_intarrs.Pointer;
                nbtasks,topdim : out natural32;
                filter,factor,verbose : out boolean ) is

  -- DESCRIPTION :
  --   Extracts the options from the arguments a and b.

  -- ON RETURN :
  --   nbtasks  number of tasks,
  --   topdim   top dimension to start the cascade,
  --   filter   if the witness supersets need filtering,
  --   factor   if the witness sets will be factored,
  --   verbose  for intermediate output.

    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    vrb : constant integer32 := integer32(v_b(v_b'first));
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));
    use Interfaces.C;

  begin
    verbose := (vrb = 1);
    nbtasks := natural32(v_a(v_a'first));
    topdim := natural32(v_a(v_a'first+1));
    filter := (natural32(v_a(v_a'first+2)) = 1);
    factor := (natural32(v_a(v_a'first+3)) = 1);
    if verbose then
      put("The number of tasks : "); put(nbtasks,1); new_line;
      put("The top dimension : "); put(topdim,1); new_line;
      if filter
       then put_line("The witness supersets will be filtered.");
       else put_line("The witness supersets will not be filtered.");
      end if;
      if factor
       then put_line("The witness supersets will be factored.");
       else put_line("The witness supersets will not be factored.");
      end if;
    end if;
  end extract_solver_options;

  procedure Store_Factors
               ( b : in C_intarrs.Pointer;
                 idxfac : in Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                 vrblvl : in integer32 := 0 ) is

  -- DESCRTIPTION :
  --   Stores the factors in idxfac for later retrieval.
  --   Assigns to the pointer b the number of characters
  --   in the string representation of the decomposition.

    strdeco : constant string
            := Path_Counts_Table.Decomposition_String(idxfac.all);

  begin
    if vrblvl > 0 then
      put_line("The irreducible factors in the decomposition :");
      put_line(strdeco);
    end if;
    Path_Counts_Table.Store_Decomposition(idxfac);
    Assign(integer32(strdeco'last),b);
  end Store_Factors;

  function Standard_Polynomial_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("Standard_Polynomial_Solver ...");
    end if;
    extract_solver_options(a,b,nbtasks,topdim,filter,factor,verbose);
    if lp /= null then
      nq := natural32(lp'last);
      nv := Number_of_Unknowns(lp(lp'first));
      lowdim := Lower_Dimension(nq,nv);
    end if;
    if verbose then
      if lp = null then
        put_line("No polynomial system in the container!?");
      else
        put_line("The polynomial system on input :"); put(lp.all);
        put("Lower bound on the dimension : "); put(lowdim,1); new_line;
      end if;
    end if;
    if lp /= null then
      Standard_Witness_Solutions.Initialize(topdim);
      Standard_Solve_with_Callback
        (nbtasks,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
      if verbose then
        Standard_Write(topdim,lowdim);
        Write_Counts(filter,factor,pc,fc,idxfac);
      end if;
      if factor
       then Store_Factors(b,idxfac,vrblvl);
       else Assign(0,b);
      end if;
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("Standard_Polynomial_Solver.");
      end if;
    return 845;
  end Standard_Polynomial_Solver;

  function Standard_Laurent_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("Standard_Laurent_Solver ...");
    end if;
    extract_solver_options(a,b,nbtasks,topdim,filter,factor,verbose);
    if lp /= null then
      nq := natural32(lp'last);
      nv := Number_of_Unknowns(lp(lp'first));
      lowdim := Lower_Dimension(nq,nv);
    end if;
    if verbose then
      if lp = null then
        put_line("No polynomial system in the container!?");
      else
        put_line("The polynomial system on input :"); put(lp.all);
        put("Lower bound on the dimension : "); put(lowdim,1); new_line;
      end if;
    end if;
    if lp /= null then
      Standard_Witness_Solutions.Initialize(topdim);
      Standard_Solve_with_Callback
        (nbtasks,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
      if verbose then
        Standard_Write(topdim,lowdim);
        Write_Counts(filter,factor,pc,fc,idxfac);
      end if;
      if factor
       then Store_Factors(b,idxfac,vrblvl);
       else Assign(0,b);
      end if;
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("Standard_Laurent_Solver.");
      end if;
    return 846;
  end Standard_Laurent_Solver;

  function DoblDobl_Polynomial_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("DoblDobl_Polynomial_Solver ...");
    end if;
    extract_solver_options(a,b,nbtasks,topdim,filter,factor,verbose);
    if lp /= null then
      nq := natural32(lp'last);
      nv := Number_of_Unknowns(lp(lp'first));
      lowdim := Lower_Dimension(nq,nv);
    end if;
    if verbose then
      if lp = null then
        put_line("No polynomial system in the container!?");
      else
        put_line("The polynomial system on input :"); put(lp.all);
        put("Lower bound on the dimension : "); put(lowdim,1); new_line;
      end if;
    end if;
    if lp /= null then
      DoblDobl_Witness_Solutions.Initialize(topdim);
      DoblDobl_Solve_with_Callback
        (nbtasks,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
      if verbose then
        DoblDobl_Write(topdim,lowdim);
        Write_Counts(filter,factor,pc,fc,idxfac);
      end if;
      if factor
       then Store_Factors(b,idxfac,vrblvl);
       else Assign(0,b);
      end if;
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("DoblDobl_Polynomial_Solver.");
      end if;
    return 847;
  end DoblDobl_Polynomial_Solver;

  function DoblDobl_Laurent_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Laur_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("DoblDobl_Laurent_Solver ...");
    end if;
    extract_solver_options(a,b,nbtasks,topdim,filter,factor,verbose);
    if lp /= null then
      nq := natural32(lp'last);
      nv := Number_of_Unknowns(lp(lp'first));
      lowdim := Lower_Dimension(nq,nv);
    end if;
    if verbose then
      if lp = null then
        put_line("No polynomial system in the container!?");
      else
        put_line("The polynomial system on input :"); put(lp.all);
        put("Lower bound on the dimension : "); put(lowdim,1); new_line;
      end if;
    end if;
    if lp /= null then
      DoblDobl_Witness_Solutions.Initialize(topdim);
      DoblDobl_Solve_with_Callback
        (nbtasks,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
      if verbose then
        DoblDobl_Write(topdim,lowdim);
        Write_Counts(filter,factor,pc,fc,idxfac);
      end if;
      if factor
       then Store_Factors(b,idxfac,vrblvl);
       else Assign(0,b);
      end if;
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("DoblDobl_Laurent_Solver.");
      end if;
    return 848;
  end DoblDobl_Laurent_Solver;

  function QuadDobl_Polynomial_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("QuadDobl_Polynomial_Solver ...");
    end if;
    extract_solver_options(a,b,nbtasks,topdim,filter,factor,verbose);
    if lp /= null then
      nq := natural32(lp'last);
      nv := Number_of_Unknowns(lp(lp'first));
      lowdim := Lower_Dimension(nq,nv);
    end if;
    if verbose then
      if lp = null then
        put_line("No polynomial system in the container!?");
      else
        put_line("The polynomial system on input :"); put(lp.all);
        put("Lower bound on the dimension : "); put(lowdim,1); new_line;
      end if;
    end if;
    if lp /= null then
      QuadDobl_Witness_Solutions.Initialize(topdim);
      QuadDobl_Solve_with_Callback
        (nbtasks,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
      if verbose then
        QuadDobl_Write(topdim,lowdim);
        Write_Counts(filter,factor,pc,fc,idxfac);
      end if;
      if factor
       then Store_Factors(b,idxfac,vrblvl);
       else Assign(0,b);
      end if;
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("QuadDobl_Polynomial_Solver.");
      end if;
    return 849;
  end QuadDobl_Polynomial_Solver;

  function QuadDobl_Laurent_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("QuadDobl_Laurent_Solver ...");
    end if;
    extract_solver_options(a,b,nbtasks,topdim,filter,factor,verbose);
    if lp /= null then
      nq := natural32(lp'last);
      nv := Number_of_Unknowns(lp(lp'first));
      lowdim := Lower_Dimension(nq,nv);
    end if;
    if verbose then
      if lp = null then
        put_line("No polynomial system in the container!?");
      else
        put_line("The polynomial system on input :"); put(lp.all);
        put("Lower bound on the dimension : "); put(lowdim,1); new_line;
      end if;
    end if;
    if lp /= null then
      QuadDobl_Witness_Solutions.Initialize(topdim);
      QuadDobl_Solve_with_Callback
        (nbtasks,topdim,lowdim,lp.all,filter,factor,pc,fc,idxfac,Store'access);
      if verbose then
        QuadDobl_Write(topdim,lowdim);
        Write_Counts(filter,factor,pc,fc,idxfac);
      end if;
      if factor
       then Store_Factors(b,idxfac,vrblvl);
       else Assign(0,b);
      end if;
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("QuadDobl_Laurent_Solver.");
      end if;
    return 850;
  end QuadDobl_Laurent_Solver;

  function Standard_Polynomial_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys
       := Standard_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := Standard_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("Standard_Polynomial_WitSet_Copy ...");
    end if;
    if lp /= null then
      Standard_PolySys_Container.Clear;
      Standard_PolySys_Container.Initialize(lp.all);
    end if;
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(ws);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("Standard_Polynomial_WitSet_Copy.");
      end if;
    return 851;
  end Standard_Polynomial_WitSet_Copy;

  function Standard_Laurent_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Laur_Sys
       := Standard_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := Standard_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("Standard_Laurent_WitSet_Copy ...");
    end if;
    if lp /= null then
      Standard_LaurSys_Container.Clear;
      Standard_LaurSys_Container.Initialize(lp.all);
    end if;
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(ws);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("Standard_Laurent_WitSet_Copy.");
      end if;
    return 852;
  end Standard_Laurent_WitSet_Copy;

  function DoblDobl_Polynomial_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys
       := DoblDobl_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := DoblDobl_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("DoblDobl_Polynomial_WitSet_Copy ...");
    end if;
    if lp /= null then
      DoblDobl_PolySys_Container.Clear;
      DoblDobl_PolySys_Container.Initialize(lp.all);
    end if;
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(ws);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("DoblDobl_Polynomial_WitSet_Copy.");
      end if;
    return 853;
  end DoblDobl_Polynomial_WitSet_Copy;

  function DoblDobl_Laurent_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Laur_Sys
       := DoblDobl_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := DoblDobl_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("DoblDobl_Laurent_WitSet_Copy ...");
    end if;
    if lp /= null then
      DoblDobl_LaurSys_Container.Clear;
      DoblDobl_LaurSys_Container.Initialize(lp.all);
    end if;
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(ws);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("DoblDobl_Laurent_WitSet_Copy.");
      end if;
    return 854;
  end DoblDobl_Laurent_WitSet_Copy;

  function QuadDobl_Polynomial_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys
       := QuadDobl_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := QuadDobl_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("QuadDobl_Polynomial_WitSet_Copy ...");
    end if;
    if lp /= null then
      QuadDobl_PolySys_Container.Clear;
      QuadDobl_PolySys_Container.Initialize(lp.all);
    end if;
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(ws);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("QuadDobl_Polynomial_WitSet_Copy.");
      end if;
    return 855;
  end QuadDobl_Polynomial_WitSet_Copy;

  function QuadDobl_Laurent_WitSet_Copy
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Laur_Sys
       := QuadDobl_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := QuadDobl_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("QuadDobl_Laurent_WitSet_Copy ...");
    end if;
    if lp /= null then
      QuadDobl_LaurSys_Container.Clear;
      QuadDobl_LaurSys_Container.Initialize(lp.all);
    end if;
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(ws);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("QuadDobl_Laurent_WitSet_Copy.");
      end if;
    return 856;
  end QuadDobl_Laurent_WitSet_Copy;

  function Standard_WitSet_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("Standard_WitSet_Clear ...");
    end if;
    Standard_Witness_Solutions.Clear;
    Path_Counts_Table.Clear_Decomposition;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("Standard_WitSet_Clear.");
      end if;
    return 857;
  end Standard_WitSet_Clear;

  function DoblDobl_WitSet_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("DoblDobl_WitSet_Clear ...");
    end if;
    DoblDobl_Witness_Solutions.Clear;
    Path_Counts_Table.Clear_Decomposition;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("DoblDobl_WitSet_Clear.");
      end if;
    return 858;
  end DoblDobl_WitSet_Clear;

  function QuadDobl_WitSet_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("QuadDobl_Witness_Solutions_Clear ...");
    end if;
    QuadDobl_Witness_Solutions.Clear;
    Path_Counts_Table.Clear_Decomposition;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("QuadDobl_WitSet_Clear.");
      end if;
    return 859;
  end QuadDobl_WitSet_Clear;

  function Irreducible_Factor_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    size : constant integer32 := integer32(v_b(v_b'first));
    deco : constant Standard_Natural_VecVecs.Link_to_Array_of_VecVecs
         := Path_Counts_Table.Get_Decomposition;

    use Standard_Natural_VecVecs;

  begin
    if vrblvl > 0 then
      put("-> in irreducible_components_interface.");
      put_line("Irreducible_Factor_String ...");
    end if;
    if deco = null then
      if size /= 0 then
        put("Retrieval of factors failed.");
        put_line("  Cannot make the factor string!");
      end if;
    else
      declare
        strdeco : constant string
                := Path_Counts_Table.Decomposition_String(deco.all);
        sv : constant Standard_Integer_Vectors.Vector
           := String_to_Integer_Vector(strdeco);
      begin
        Assign(integer32(sv'last),a);
        Assign(sv,b);
      end;
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in irreducible_components_interface.");
        put_line("Irreducible_Factor_String.");
      end if;
    return 993;
  end Irreducible_Factor_String;

end Irreducible_Components_Interface;
