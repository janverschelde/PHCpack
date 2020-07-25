with Interfaces.C;
with text_io;                              use text_io;
with Standard_Natural_Numbers;             use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;          use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors_io;          use Standard_Complex_Vectors_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  -- for debugging
 use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Standard_Homotopy;
with Standard_Complex_Solutions;
with DoblDobl_Homotopy;
with DoblDobl_Complex_Solutions;
with QuadDobl_Homotopy;
with QuadDobl_Complex_Solutions;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_Vectors_io;
with Standard_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors_io;
with DoblDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors_io;
with QuadDobl_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io; -- for debugging
with Series_and_Solutions;
with Power_Series_Methods;                 use Power_Series_Methods;
with Standard_Pade_Approximants;
with Standard_Pade_Approximants_io;
with DoblDobl_Pade_Approximants;
with DoblDobl_Pade_Approximants_io;
with QuadDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants_io;
with Homotopy_Pade_Approximants;
with Standard_PolySys_Container;
with Standard_Systems_Pool;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_Systems_Pool;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_Systems_Pool;
with QuadDobl_Solutions_Container;

package body Power_Series_Interface is

  procedure extract_options
              ( a : in C_intarrs.Pointer; b : in C_intarrs.Pointer;
                idx,maxdeg,nbrit : out integer32; verbose : out boolean ) is

  -- DESCRIPTION :
  --   Extracts the options from the arguments a and b.

    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant integer32 := integer32(v_b(v_b'first));
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    use Interfaces.C;

  begin
    verbose := (vrb = 1);
    idx := integer32(v_a(v_a'first));
    maxdeg := integer32(v_a(v_a'first+1));
    nbrit := integer32(v_a(v_a'first+2));
    if verbose then
      put("The index of the series parameter : "); put(idx,1); new_line;
      put("The maximal degree of the series : "); put(maxdeg,1); new_line;
      put("The number of Newton steps : "); put(nbrit,1); new_line;
    end if;
  end extract_options;

  procedure Load_Series_Solutions
              ( s : out Standard_Complex_Series_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Loads the series s from the systems in the systems pool,
  --   converting every polynomial in the pool to a series.

    len : constant integer32 := integer32(Standard_Systems_Pool.Size);
    res : Standard_Complex_Series_VecVecs.VecVec(1..len);

  begin
    for k in 1..len loop
      declare
        lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
           := Standard_Systems_Pool.Retrieve(k);
        sv : constant Standard_Complex_Series_Vectors.Vector
           := Complex_Series_and_Polynomials.System_to_Series_Vector(lp.all);
      begin
        res(k) := new Standard_Complex_Series_Vectors.Vector'(sv);
      end;
    end loop;
    s := new Standard_Complex_Series_VecVecs.VecVec'(res);
  end Load_Series_Solutions;

  procedure Load_Series_Solutions
              ( s : out DoblDobl_Complex_Series_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Loads the series s from the systems in the systems pool,
  --   converting every polynomial in the pool to a series.

    len : constant integer32 := integer32(DoblDobl_Systems_Pool.Size);
    res : DoblDobl_Complex_Series_VecVecs.VecVec(1..len);

  begin
    for k in 1..len loop
      declare
        lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
           := DoblDobl_Systems_Pool.Retrieve(k);
        sv : constant DoblDobl_Complex_Series_Vectors.Vector
           := Complex_Series_and_Polynomials.System_to_Series_Vector(lp.all);
      begin
        res(k) := new DoblDobl_Complex_Series_Vectors.Vector'(sv);
      end;
    end loop;
    s := new DoblDobl_Complex_Series_VecVecs.VecVec'(res);
  end Load_Series_Solutions;

  procedure Load_Series_Solutions
              ( s : out QuadDobl_Complex_Series_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Loads the series s from the systems in the systems pool,
  --   converting every polynomial in the pool to a series.

    len : constant integer32 := integer32(QuadDobl_Systems_Pool.Size);
    res : QuadDobl_Complex_Series_VecVecs.VecVec(1..len);

  begin
    for k in 1..len loop
      declare
        lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
           := QuadDobl_Systems_Pool.Retrieve(k);
        sv : constant QuadDobl_Complex_Series_Vectors.Vector
           := Complex_Series_and_Polynomials.System_to_Series_Vector(lp.all);
      begin
        res(k) := new QuadDobl_Complex_Series_Vectors.Vector'(sv);
      end;
    end loop;
    s := new QuadDobl_Complex_Series_VecVecs.VecVec'(res);
  end Load_Series_Solutions;

  procedure Store_Series_Solutions
              ( s : in Standard_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Stores the solutions in s in the systems pool,
  --   after converting every series vector to a polynomial system.

  begin
    Standard_Systems_Pool.Initialize(s'last);
    for k in s'range loop
      declare
        p : constant Standard_Complex_Poly_Systems.Poly_Sys
          := Complex_Series_and_Polynomials.Series_Vector_to_System(s(k).all);
      begin
        Standard_Systems_Pool.Initialize(k,p);
      end;
    end loop;
  end Store_Series_Solutions;

  procedure Store_Series_Solutions
              ( s : in DoblDobl_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Stores the solutions in s in the systems pool,
  --   after converting every series vector to a polynomial system.

  begin
    DoblDobl_Systems_Pool.Initialize(s'last);
    for k in s'range loop
      declare
        p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
          := Complex_Series_and_Polynomials.Series_Vector_to_System(s(k).all);
      begin
        DoblDobl_Systems_Pool.Initialize(k,p);
      end;
    end loop;
  end Store_Series_Solutions;

  procedure Store_Series_Solutions
              ( s : in QuadDobl_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Stores the solutions in s in the systems pool,
  --   after converting every series vector to a polynomial system.

  begin
    QuadDobl_Systems_Pool.Initialize(s'last);
    for k in s'range loop
      declare
        p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
          := Complex_Series_and_Polynomials.Series_Vector_to_System(s(k).all);
      begin
        QuadDobl_Systems_Pool.Initialize(k,p);
      end;
    end loop;
  end Store_Series_Solutions;

  procedure Run_Newton
              ( nq,idx,dim,maxdeg,nbrit : in integer32;
                echelon,verbose : in boolean;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in standard double precision.

  -- ON ENTRY :
  --   nq       number of equations in p;
  --   idx      index to the series parameter;
  --   dim      the number of coordinates in the series;
  --   maxdeg   the maximal degree of the series;
  --   nbrit    the number of Newton steps;
  --   echelon  flag to indicate if echelon form needs to be used;
  --   verbose  flag to print message to screen;
  --   p        a polynomial of nq equations in nv unknowns;
  --   s        start terms in a solution series.

  begin
    if echelon then
      if verbose
       then put_line("Echelon Newton will be applied.");
      end if;
      Run_Echelon_Newton(maxdeg,nbrit,p,s,verbose);
    else
      if nq = dim then
        if verbose
         then put_line("LU Newton will be applied.");
        end if;
        Run_LU_Newton(maxdeg,nbrit,p,s,verbose);
      else
        if verbose
         then put_line("QR Newton will be applied.");
        end if;
        Run_QR_Newton(maxdeg,nbrit,p,s,verbose);
      end if;
    end if;
  end Run_Newton;

  procedure Run_Newton
              ( nq,idx,dim,maxdeg,nbrit : in integer32;
                echelon,verbose : in boolean;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in double double precision.

  -- ON ENTRY :
  --   nq       number of equations in p;
  --   idx      index to the series parameter;
  --   dim      the number of coordinates in the series;
  --   maxdeg   the maximal degree of the series;
  --   nbrit    the number of Newton steps;
  --   echelon  flag to indicate whether to use echelon Newton;
  --   verbose  flag for printing intermediate output;
  --   p        a polynomial of nq equations in nv unknowns;
  --   s        start terms in a solution series.

  begin
    if echelon then
      if verbose
       then put_line("Echelon Newton will be applied.");
      end if;
      Run_Echelon_Newton(maxdeg,nbrit,p,s,verbose);
    else
      if nq = dim then
        if verbose
         then put_line("LU Newton will be applied.");
        end if;
        Run_LU_Newton(maxdeg,nbrit,p,s,verbose);
      else
        if verbose
         then put_line("QR Newton will be applied.");
        end if;
        Run_QR_Newton(maxdeg,nbrit,p,s,verbose);
      end if;
    end if;
  end Run_Newton;

  procedure Run_Newton
              ( nq,idx,dim,maxdeg,nbrit : in integer32;
                echelon,verbose : in boolean;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in quad double precision.

  -- ON ENTRY :
  --   nq       number of equations in p;
  --   idx      index to the series parameter;
  --   dim      the number of coordinates in the series;
  --   maxdeg   the maximal degree of the series;
  --   nbrit    the number of Newton steps;
  --   echelon  flag to indicate whether to use echelon Newton;
  --   verbose  flag for intermediate output;
  --   p        a polynomial of nq equations in nv unknowns;
  --   s        start terms in a solution series.

  begin
    if echelon then
      if verbose
       then put_line("Echelon Newton will be applied.");
      end if;
      Run_Echelon_Newton(maxdeg,nbrit,p,s,verbose);
    else
      if nq = dim then
        if verbose
         then put_line("LU Newton will be applied.");
        end if;
        Run_LU_Newton(maxdeg,nbrit,p,s,verbose);
      else
        if verbose
         then put_line("QR Newton will be applied.");
        end if;
        Run_QR_Newton(maxdeg,nbrit,p,s,verbose);
      end if;
    end if;
  end Run_Newton;

  procedure Run_Newton
              ( nq,idx,dim,maxdeg,nbrit : in integer32; verbose : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in standard double precision.

  -- ON ENTRY :
  --   nq       number of equations in p;
  --   idx      index to the series parameter;
  --   dim      the number of coordinates in the series;
  --   maxdeg   the maximal degree of the series;
  --   nbrit    the number of Newton steps;
  --   verbose  for additional output;
  --   p        a polynomial of nq equations in nv unknowns;
  --   s        a list of solutions.

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : Standard_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : Standard_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(nq,idx,dim,maxdeg,nbrit,true,verbose,srp,srv);
    Store_Series_Solutions(srv);
    Standard_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton;

  procedure Run_Newton
              ( nq,idx,dim,maxdeg,nbrit : in integer32; verbose : in boolean;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in double double precision.

  -- ON ENTRY :
  --   nq       number of equations in p;
  --   idx      index to the series parameter;
  --   dim      the number of coordinates in the series;
  --   maxdeg   the maximal degree of the series;
  --   nbrit    the number of Newton steps;
  --   verbose  for additional output;
  --   p        a polynomial of nq equations in nv unknowns;
  --   s        a list of solutions.

    use DoblDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : DoblDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(nq,idx,dim,maxdeg,nbrit,true,verbose,srp,srv);
    Store_Series_Solutions(srv);
    DoblDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton;

  procedure Run_Newton
              ( nq,idx,dim,maxdeg,nbrit : in integer32; verbose : in boolean;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in quad double precision.

  -- ON ENTRY :
  --   nq       number of equations in p;
  --   idx      index to the series parameter;
  --   dim      the number of coordinates in the series;
  --   maxdeg   the maximal degree of the series;
  --   nbrit    the number of Newton steps;
  --   verbose  for additional output;
  --   p        a polynomial of nq equations in nv unknowns;
  --   s        a list of solutions.

    use QuadDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : QuadDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(nq,idx,dim,maxdeg,nbrit,true,verbose,srp,srv);
    Store_Series_Solutions(srv);
    QuadDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton;

  function Series_Standard_Newton_at_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    nq : constant integer32 := lp'last;
    nv : constant integer32 := Head_Of(sols).n;
    idx,maxdeg,nbr,dim : integer32;
    verbose : boolean;

  begin
    if vrblvl > 0 then
      put("-> in power_series_interface.");
      put_line("Series_Standard_Newton_at_Point ...");
    end if;
    extract_options(a,b,idx,maxdeg,nbr,verbose);
    dim := (if idx = 0 then nv else nv-1);
    if verbose then
      put("Number of equations in the system : "); put(nq,1); new_line;
      put("The dimension of the series : "); put(dim,1); new_line;
    end if;
    Run_Newton(nq,idx,dim,maxdeg,nbr,verbose,lp.all,sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in power_series_interface.");
        put_line("Series_Standard_Newton_at_Point.");
      end if;
      return 691;
  end Series_Standard_Newton_at_Point;

  function Series_DoblDobl_Newton_at_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    nq : constant integer32 := lp'last;
    nv : constant integer32 := Head_Of(sols).n;
    idx,maxdeg,nbr,dim : integer32;
    verbose : boolean;

  begin
    if vrblvl > 0 then
      put("-> in power_series_interface.");
      put_line("Series_DoblDobl_Newton_at_Point ...");
    end if;
    extract_options(a,b,idx,maxdeg,nbr,verbose);
    dim := (if idx = 0 then nv else nv-1);
    if verbose then
      put("Number of equations in the system : "); put(nq,1); new_line;
      put("The dimension of the series : "); put(dim,1); new_line;
    end if;
    Run_Newton(nq,idx,dim,maxdeg,nbr,verbose,lp.all,sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in power_series_interface.");
        put_line("Series_DoblDobl_Newton_at_Point.");
      end if;
      return 692;
  end Series_DoblDobl_Newton_at_Point;

  function Series_QuadDobl_Newton_at_Point
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    nq : constant integer32 := lp'last;
    nv : constant integer32 := Head_Of(sols).n;
    idx,maxdeg,nbr,dim : integer32;
    verbose : boolean;

  begin
    if vrblvl > 0 then
      put("-> in power_series_interface.");
      put_line("Series_QuadDobl_Newton_at_Point ...");
    end if;
    extract_options(a,b,idx,maxdeg,nbr,verbose);
    dim := (if idx = 0 then nv else nv-1);
    if verbose then
      put("Number of equations in the system : "); put(nq,1); new_line;
      put("The dimension of the series : "); put(dim,1); new_line;
    end if;
    Run_Newton(nq,idx,dim,maxdeg,nbr,verbose,lp.all,sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in power_series_interface.");
        put_line("Series_QuadDobl_Newton_at_Point.");
      end if;
      return 693;
  end Series_QuadDobl_Newton_at_Point;

  function Series_Standard_Newton_at_Series
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    nq : constant integer32 := lp'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(lp(lp'first)));
    idx,maxdeg,nbr,dim : integer32;
    verbose : boolean;
    srv : Standard_Complex_Series_VecVecs.Link_to_VecVec;
    srp : Standard_CSeries_Poly_Systems.Poly_Sys(lp'range);

  begin
    if vrblvl > 0 then
      put("-> in power_series_interface.");
      put_line("Series_Standard_Newton_at_Series ...");
    end if;
    Load_Series_Solutions(srv);
    extract_options(a,b,idx,maxdeg,nbr,verbose);
    srp := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    dim := (if idx = 0 then nv else nv-1);
    if verbose then
      put("Number of equations in the system : "); put(nq,1); new_line;
      put("The dimension of the series : "); put(dim,1); new_line;
      put("The index of the parameter : "); put(idx,1); new_line;
      put_line("The polynomials in the system :"); put(lp.all);
      put_line("The system converted to series :");
      Complex_Series_and_Polynomials_io.put(srp);
    end if;
    Run_Newton(nq,idx,dim,maxdeg,nbr,true,verbose,srp,srv.all);
    Store_Series_Solutions(srv.all);
    Standard_CSeries_Poly_Systems.Clear(srp);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in power_series_interface.");
        put_line("Series_Standard_Newton_at_Series.");
      end if;
      return 694;
  end Series_Standard_Newton_at_Series;

  function Series_DoblDobl_Newton_at_Series
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    nq : constant integer32 := lp'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(lp(lp'first)));
    idx,maxdeg,nbr,dim : integer32;
    verbose : boolean;
    srv : DoblDobl_Complex_Series_VecVecs.Link_to_VecVec;
    srp : DoblDobl_CSeries_Poly_Systems.Poly_Sys(lp'range);

  begin
    if vrblvl > 0 then
      put("-> in power_series_interface.");
      put_line("Series_DoblDobl_Newton_at_Series ...");
    end if;
    Load_Series_Solutions(srv);
    extract_options(a,b,idx,maxdeg,nbr,verbose);
    srp := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    dim := (if idx = 0 then nv else nv-1);
    if verbose then
      put("Number of equations in the system : "); put(nq,1); new_line;
      put("The dimension of the series : "); put(dim,1); new_line;
    end if;
    Run_Newton(nq,idx,dim,maxdeg,nbr,true,verbose,srp,srv.all);
    Store_Series_Solutions(srv.all);
    DoblDobl_CSeries_Poly_Systems.Clear(srp);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in power_series_interface.");
        put_line("Series_Standard_Newton_at_Series.");
      end if;
      return 695;
  end Series_DoblDobl_Newton_at_Series;

  function Series_QuadDobl_Newton_at_Series
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    nq : constant integer32 := lp'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(lp(lp'first)));
    idx,maxdeg,nbr,dim : integer32;
    verbose : boolean;
    srv : QuadDobl_Complex_Series_VecVecs.Link_to_VecVec;
    srp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(lp'range);

  begin
    if vrblvl > 0 then
      put("-> in power_series_interface.");
      put_line("Series_QuadDobl_Newton_at_Series ...");
    end if;
    Load_Series_Solutions(srv);
    extract_options(a,b,idx,maxdeg,nbr,verbose);
    srp := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    dim := (if idx = 0 then nv else nv-1);
    if verbose then
      put("Number of equations in the system : "); put(nq,1); new_line;
      put("The dimension of the series : "); put(dim,1); new_line;
    end if;
    Run_Newton(nq,idx,dim,maxdeg,nbr,true,verbose,srp,srv.all);
    Store_Series_Solutions(srv.all);
    QuadDobl_CSeries_Poly_Systems.Clear(srp);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in power_series_interface.");
        put_line("Series_QuadDobl_Newton_at_Series.");
      end if;
      return 696;
  end Series_QuadDobl_Newton_at_Series;

  procedure extract_options_for_pade
              ( a : in C_intarrs.Pointer; b : in C_intarrs.Pointer;
                idx,numdeg,dendeg,nbrit : out integer32; 
                verbose : out boolean ) is

  -- DESCRIPTION :
  --   Extracts the options from the arguments a and b.

    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant integer32 := integer32(v_b(v_b'first));
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));
    use Interfaces.C;

  begin
    verbose := (vrb = 1);
    idx := integer32(v_a(v_a'first));
    numdeg := integer32(v_a(v_a'first+1));
    dendeg := integer32(v_a(v_a'first+2));
    nbrit := integer32(v_a(v_a'first+3));
    if verbose then
      put("The index of the series parameter : "); put(idx,1); new_line;
      put("The degree of the numerator : "); put(numdeg,1); new_line;
      put("The degree of the denominator : "); put(dendeg,1); new_line;
      put("The number of Newton steps : "); put(nbrit,1); new_line;
    end if;
  end extract_options_for_pade;

  function Coordinates ( s : Standard_Complex_Solutions.Link_to_Solution;
                         idx : integer32 )
                       return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns s.v if idx equals zero, otherwise return s.v with the
  --   coordinate with index idx removed.

  -- REQUIRED : s /= null;

  begin
    if idx = 0 then
      return s.v;
    elsif idx = s.v'last then
      return s.v(s.v'first..s.v'last-1);
    else
      declare
        res : Standard_Complex_Vectors.Vector(s.v'first..s.v'last-1);
      begin
        for k in 1..(idx-1) loop
          res(k) := s.v(k);
        end loop;
        for k in (idx+1)..s.v'last loop
          res(k-1) := s.v(k);
        end loop;
        return res;
      end;
    end if;
  end Coordinates;

  function Coordinates ( s : DoblDobl_Complex_Solutions.Link_to_Solution;
                         idx : integer32 )
                       return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns s.v if idx equals zero, otherwise return s.v with the
  --   coordinate with index idx removed.

  -- REQUIRED : s /= null;

  begin
    if idx = 0 then
      return s.v;
    elsif idx = s.v'last then
      return s.v(s.v'first..s.v'last-1);
    else
      declare
        res : DoblDobl_Complex_Vectors.Vector(s.v'first..s.v'last-1);
      begin
        for k in 1..(idx-1) loop
          res(k) := s.v(k);
        end loop;
        for k in (idx+1)..s.v'last loop
          res(k-1) := s.v(k);
        end loop;
        return res;
      end;
    end if;
  end Coordinates;

  function Coordinates ( s : QuadDobl_Complex_Solutions.Link_to_Solution;
                         idx : integer32 )
                       return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns s.v if idx equals zero, otherwise return s.v with the
  --   coordinate with index idx removed.

  -- REQUIRED : s /= null;

  begin
    if idx = 0 then
      return s.v;
    elsif idx = s.v'last then
      return s.v(s.v'first..s.v'last-1);
    else
      declare
        res : QuadDobl_Complex_Vectors.Vector(s.v'first..s.v'last-1);
      begin
        for k in 1..(idx-1) loop
          res(k) := s.v(k);
        end loop;
        for k in (idx+1)..s.v'last loop
          res(k-1) := s.v(k);
        end loop;
        return res;
      end;
    end if;
  end Coordinates;

  function Series_Standard_Pade
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    nq : constant integer32 := lp'last;
    nv : constant integer32 := Head_Of(sols).n;
    idx,numdeg,dendeg,nbr,dim : integer32;
    verbose : boolean;
    tmp : Solution_List;
    cnt : integer32 := 0;

  begin
    if vrblvl > 0 then
      put_line("-> in power_series_interface.Series_Standard_Pade ...");
    end if;
    extract_options_for_pade(a,b,idx,numdeg,dendeg,nbr,verbose);
    dim := (if idx = 0 then nv else nv-1);
    if verbose then
      put("Number of equations in the system : "); put(nq,1); new_line;
      put("Number of variables in the solutions : "); put(nv,1); new_line;
      put("The dimension of the series : "); put(dim,1); new_line;
    end if;
    Standard_Homotopy.Create(lp.all,idx);
    Standard_Systems_Pool.Initialize(integer32(Length_Of(sols)));
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
        sol : constant Link_to_Solution := Head_Of(tmp);
        solvec : constant Standard_Complex_Vectors.Vector
               := Coordinates(sol,idx);
        nbt : constant natural32 := natural32(numdeg+dendeg+1);
        nit : constant natural32 := 4*nbt;
        srv : Standard_Complex_Series_Vectors.Vector(1..dim);
        eva : Standard_Complex_Series_Vectors.Vector(1..nq);
        pv : Standard_Pade_Approximants.Pade_Vector(srv'range);
      begin
        if verbose then
          put_line("Calling the Pade approximant constructor ...");
          put_line("The solution vector :"); put_line(sol.v);
        end if;
        Homotopy_Pade_Approximants.Standard_Pade_Approximant
          (solvec,idx,nq,numdeg,dendeg,nit,srv,eva,pv);
        if verbose then
          put_line("The solution series :");
          Standard_Complex_Series_Vectors_io.put(srv);
          put_line("The evaluated solution series :");
          Standard_Complex_Series_Vectors_io.put(eva);
          put_line("The Pade approximant :");
          for i in pv'range loop
            put_line(Standard_Pade_Approximants_io.Write(pv(i)));
          end loop;
        end if;
        cnt := cnt + 1;
        if verbose
         then put("Storing Pade approximant "); put(cnt,1); put_line(" ...");
        end if;
        declare
          spv : constant Standard_Complex_Poly_Systems.Poly_Sys
              := Standard_Pade_Approximants_io.to_System(pv);
        begin
          Standard_Systems_Pool.Initialize(cnt,spv);
        end;
        if verbose
         then put("Done storing Pade approximant "); put(cnt,1); new_line;
        end if;
        tmp := Tail_Of(tmp);
      end;
    end loop;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in power_series_interface.");
        put_line("Series_QuadDobl_Pade.");
      end if;
      return 704;
  end Series_Standard_Pade;

  function Series_DoblDobl_Pade
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    nq : constant integer32 := lp'last;
    nv : constant integer32 := Head_Of(sols).n;
    idx,numdeg,dendeg,nbr,dim : integer32;
    verbose : boolean;
    tmp : Solution_List;
    cnt : integer32 := 0;

  begin
    if vrblvl > 0 then
      put_line("-> in power_series_interface.Series_DoblDobl_Pade ...");
    end if;
    extract_options_for_pade(a,b,idx,numdeg,dendeg,nbr,verbose);
    dim := (if idx = 0 then nv else nv-1);
    if verbose then
      put("Number of equations in the system : "); put(nq,1); new_line;
      put("Number of variables in the solutions : "); put(nv,1); new_line;
      put("The dimension of the series : "); put(dim,1); new_line;
    end if;
    DoblDobl_Homotopy.Create(lp.all,idx);
    DoblDobl_Systems_Pool.Initialize(integer32(Length_Of(sols)));
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
        sol : constant Link_to_Solution := Head_Of(tmp);
        solvec : constant DoblDobl_Complex_Vectors.Vector
               := Coordinates(sol,idx);
        nbt : constant natural32 := natural32(numdeg+dendeg+1);
        nit : constant natural32 := 4*nbt;
        srv : DoblDobl_Complex_Series_Vectors.Vector(1..dim);
        eva : DoblDobl_Complex_Series_Vectors.Vector(1..nq);
        pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range);
      begin
        Homotopy_Pade_Approximants.DoblDobl_Pade_Approximant
          (solvec,idx,nq,numdeg,dendeg,nit,srv,eva,pv);
        if verbose then
          put_line("The solution series :");
          DoblDobl_Complex_Series_Vectors_io.put(srv);
          put_line("The evaluated solution series :");
          DoblDobl_Complex_Series_Vectors_io.put(eva);
          put_line("The Pade approximant :");
          for i in pv'range loop
            put_line(DoblDobl_Pade_Approximants_io.Write(pv(i)));
          end loop;
        end if;
        cnt := cnt + 1;
        declare
          spv : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
              := DoblDobl_Pade_Approximants_io.to_System(pv);
        begin
          DoblDobl_Systems_Pool.Initialize(cnt,spv);
        end;
        tmp := Tail_Of(tmp);
      end;
    end loop;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in power_series_interface.");
        put_line("Series_DoblDobl_Pade.");
      end if;
      return 705;
  end Series_DoblDobl_Pade;

  function Series_QuadDobl_Pade
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    nq : constant integer32 := lp'last;
    nv : constant integer32 := Head_Of(sols).n;
    idx,numdeg,dendeg,nbr,dim : integer32;
    verbose : boolean;
    tmp : Solution_List;
    cnt : integer32 := 0;

  begin
    if vrblvl > 0 then
      put_line("-> in power_series_interface.Series_QuadDobl_Pade ...");
    end if;
    extract_options_for_pade(a,b,idx,numdeg,dendeg,nbr,verbose);
    dim := (if idx = 0 then nv else nv-1);
    if verbose then
      put("Number of equations in the system : "); put(nq,1); new_line;
      put("Number of variables in the solutions : "); put(nv,1); new_line;
      put("The dimension of the series : "); put(dim,1); new_line;
    end if;
    QuadDobl_Homotopy.Create(lp.all,idx);
    QuadDobl_Systems_Pool.Initialize(integer32(Length_Of(sols)));
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
        sol : constant Link_to_Solution := Head_Of(tmp);
        solvec : constant QuadDobl_Complex_Vectors.Vector
               := Coordinates(sol,idx);
        nbt : constant natural32 := natural32(numdeg+dendeg+1);
        nit : constant natural32 := 4*nbt;
        srv : QuadDobl_Complex_Series_Vectors.Vector(1..dim);
        eva : QuadDobl_Complex_Series_Vectors.Vector(1..nq);
        pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range);
      begin
        Homotopy_Pade_Approximants.QuadDobl_Pade_Approximant
          (solvec,idx,nq,numdeg,dendeg,nit,srv,eva,pv);
        if verbose then
          put_line("The solution series :");
          QuadDobl_Complex_Series_Vectors_io.put(srv);
          put_line("The evaluated solution series :");
          QuadDobl_Complex_Series_Vectors_io.put(eva);
          put_line("The Pade approximant :");
          for i in pv'range loop
            put_line(QuadDobl_Pade_Approximants_io.Write(pv(i)));
          end loop;
        end if;
        cnt := cnt + 1;
        declare
          spv : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
              := QuadDobl_Pade_Approximants_io.to_System(pv);
        begin
          QuadDobl_Systems_Pool.Initialize(cnt,spv);
        end;
        tmp := Tail_Of(tmp);
      end;
    end loop;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in power_series_interface.");
        put_line("Series_QuadDobl_Pade.");
      end if;
      return 706;
  end Series_QuadDobl_Pade;

end Power_Series_Interface;
