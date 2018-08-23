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
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Standard_LaurSys_Container;
with DoblDobl_LaurSys_Container;
with QuadDobl_LaurSys_Container;
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

  function do_jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- solve polynomial system with standard doubles
      when 1 => return Job1; -- solve Laurent system with standard doubles
      when 2 => return Job2; -- solve polynomial system with double doubles
      when 3 => return Job3; -- solve Laurent system with double doubles
      when 4 => return Job4; -- solve polynomial system with quad doubles
      when 5 => return Job5; -- solve Laurent system with quad doubles
      when others => return -1;
    end case;
  end do_jobs;

begin
  return do_jobs;
end use_witsols;
