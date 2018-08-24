with Interfaces.C;
with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Natural_VecVecs;
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

function use_witsols ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32 is

  procedure extract_solver_options
              ( nbtasks,topdim : out natural32;
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

  function Job0 return integer32 is

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a polynomial system, in standard
  --   double precision, eventually followed by a filter and a factor.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    extract_solver_options(nbtasks,topdim,filter,factor,verbose);
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
    end if;
    return 0;
  end Job0;

  function Job1 return integer32 is

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a Laurent system, in standard
  --   double precision, eventually followed by a filter and a factor.

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    extract_solver_options(nbtasks,topdim,filter,factor,verbose);
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
    end if;
    return 0;
  end Job1;

  function Job2 return integer32 is

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a polynomial system, in double
  --   double precision, eventually followed by a filter and a factor.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    extract_solver_options(nbtasks,topdim,filter,factor,verbose);
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
    end if;
    return 0;
  end Job2;

  function Job3 return integer32 is

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a Laurent system, in double
  --   double precision, eventually followed by a filter and a factor.

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Laur_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    extract_solver_options(nbtasks,topdim,filter,factor,verbose);
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
    end if;
    return 0;
  end Job3;

  function Job4 return integer32 is

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a polynomial system, in quad
  --   double precision, eventually followed by a filter and a factor.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    extract_solver_options(nbtasks,topdim,filter,factor,verbose);
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
    end if;
    return 0;
  end Job4;

  function Job5 return integer32 is

  -- DESCRIPTION :
  --   Runs the cascade homotopies on a Laurent system, in quad
  --   double precision, eventually followed by a filter and a factor.

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Systems;

    nbtasks,topdim,lowdim,nq,nv : natural32;
    filter,factor,verbose : boolean;
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    pc,fc : Standard_Natural_VecVecs.Link_to_VecVec;
    idxfac : Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;

  begin
    extract_solver_options(nbtasks,topdim,filter,factor,verbose);
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
    end if;
    return 0;
  end Job5;

  function Job6 return integer32 is

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a polynomial system in
  --   standard double precision and copies it into the containers.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys
       := Standard_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := Standard_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if lp /= null then
      Standard_PolySys_Container.Clear;
      Standard_PolySys_Container.Initialize(lp.all);
    end if;
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(ws);
    return 0;
  end Job6;

  function Job7 return integer32 is

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a Laurent polynomial system in
  --   standard double precision and copies it into the containers.

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Laur_Sys
       := Standard_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := Standard_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if lp /= null then
      Standard_LaurSys_Container.Clear;
      Standard_LaurSys_Container.Initialize(lp.all);
    end if;
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(ws);
    return 0;
  end Job7;

  function Job8 return integer32 is

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a polynomial system in
  --   double double precision and copies it into the containers.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys
       := DoblDobl_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := DoblDobl_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if lp /= null then
      DoblDobl_PolySys_Container.Clear;
      DoblDobl_PolySys_Container.Initialize(lp.all);
    end if;
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(ws);
    return 0;
  end Job8;

  function Job9 return integer32 is

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a Laurent polynomial system in
  --   double double precision and copies it into the containers.

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Laur_Sys
       := DoblDobl_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := DoblDobl_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if lp /= null then
      DoblDobl_LaurSys_Container.Clear;
      DoblDobl_LaurSys_Container.Initialize(lp.all);
    end if;
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(ws);
    return 0;
  end Job9;

  function Job10 return integer32 is

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a polynomial system in
  --   quad double precision and copies it into the containers.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys
       := QuadDobl_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := QuadDobl_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if lp /= null then
      QuadDobl_PolySys_Container.Clear;
      QuadDobl_PolySys_Container.Initialize(lp.all);
    end if;
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(ws);
    return 0;
  end Job10;

  function Job11 return integer32 is

  -- DESCRIPTION :
  --   Retrieves a witness set stored for a Laurent polynomial system in
  --   quad double precision and copies it into the containers.

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Laur_Sys
       := QuadDobl_Witness_Solutions.Load_Embedded_System(dim);
    ws : constant Solution_List
       := QuadDobl_Witness_Solutions.Load_Witness_Points(dim);

  begin
    if lp /= null then
      QuadDobl_LaurSys_Container.Clear;
      QuadDobl_LaurSys_Container.Initialize(lp.all);
    end if;
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(ws);
    return 0;
  end Job11;

  function Job12 return integer32 is

  -- DESCRIPTION :
  --   Deallocates the witness solutions in standard double precision.

  begin
    Standard_Witness_Solutions.Clear;
    return 0;
  end Job12;

  function Job13 return integer32 is

  -- DESCRIPTION :
  --   Deallocates the witness solutions in double double precision.

  begin
    DoblDobl_Witness_Solutions.Clear;
    return 0;
  end Job13;

  function Job14 return integer32 is

  -- DESCRIPTION :
  --   Deallocates the witness solutions in quad double precision.

  begin
    QuadDobl_Witness_Solutions.Clear;
    return 0;
  end Job14;

  function do_jobs return integer32 is
  begin
    case job is
      when  0 => return Job0; -- solve polynomial system with standard doubles
      when  1 => return Job1; -- solve Laurent system with standard doubles
      when  2 => return Job2; -- solve polynomial system with double doubles
      when  3 => return Job3; -- solve Laurent system with double doubles
      when  4 => return Job4; -- solve polynomial system with quad doubles
      when  5 => return Job5; -- solve Laurent system with quad doubles
      when  6 => return Job6; -- copy standard polysys witness set
      when  7 => return Job7; -- copy standard laursys witness set
      when  8 => return Job8; -- copy dobldobl polysys witness set
      when  9 => return Job9; -- copy dobldobl laursys witness set
      when 10 => return Job10; -- copy quaddobl polysys witness set
      when 11 => return Job11; -- copy quaddobl laursys witness set
      when 12 => return Job12; -- clear standard witness solutions
      when 13 => return Job13; -- clear dobldobl witness solutions
      when 14 => return Job14; -- clear quaddobl witness solutions
      when others => return -1;
    end case;
  end do_jobs;

begin
  return do_jobs;
end use_witsols;
