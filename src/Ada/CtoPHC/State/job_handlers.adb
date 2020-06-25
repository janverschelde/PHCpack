with Interfaces.C;
with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Laur_Poly_Convertors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Floating_Mixed_Subdivisions;
with Black_Mixed_Volume_Computations;
with Black_Box_Solvers;
with Greeting_Banners;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with Standard_LaurSys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_LaurSys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_LaurSys_Container;
with QuadDobl_Solutions_Container;
with Cells_Container;

package body Job_Handlers is

  function Version_String ( a : C_intarrs.Pointer;
                            b : C_intarrs.Pointer;
                            vrblvl : integer32 := 0 ) return integer32 is

    s : constant string := Greeting_Banners.Version;
    n : constant integer32 := integer32(s'last);
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(s);

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Version_String");
    end if;
    Assign(n,a);
    Assign(sv,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in job_handlers.Version_String.");
      end if;
      return 999;
  end Version_String;

  function Get_Seed ( a : C_intarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

    seed : constant integer32 := Standard_Random_Numbers.Get_Seed;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Get_Seed");
    end if;
    Assign(seed,a);
    return 0;
  exception
    when others =>
      if vrblvl > 0
       then put_line("Exception raised in job_handlers.Get_Seed.");
      end if;
      return 997;
  end Get_Seed;

  function Set_Seed ( a : C_intarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    fixed_seed : constant natural32 := natural32(v_a(v_a'first));

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Set_Seed");
    end if;
    Standard_Random_Numbers.Set_Seed(fixed_seed);
    return 0;
  exception
    when others =>
      if vrblvl > 0
       then put_line("Exception raised in job_handlers.Set_Seed.");
      end if;
      return 998;
  end Set_Seed;

  function Standard_Polynomial_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32 is

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;
    use Interfaces.C;

   -- n : constant natural := Standard_PolySys_Container.Dimension;
    v_b : constant C_Integer_Array(0..1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc,nr : natural32;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    sols : Solution_List;
    lsroco : Link_to_String;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Standard_Polynomial_Solver");
    end if;
   -- put("Dimension of the system in the container : "); put(n,1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 77;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 77;
    end if;
   -- Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
    if silent then
      if ntasks = 0
       then Black_Box_Solvers.Solve(lp.all,silent,true,rc,sols);
       else Black_Box_Solvers.Solve(ntasks,lp.all,silent,true,rc,sols);
      end if;
    else
      if ntasks = 0
       then Black_Box_Solvers.Solve(lp.all,true,rc,lsroco,sols);
       else Black_Box_Solvers.Solve(ntasks,lp.all,true,rc,lsroco,sols);
      end if;
      if lsroco = null then
        nr := 0;
      else
        nr := natural32(lsroco'last);
        declare
          sv : constant Standard_Integer_Vectors.Vector
             := String_to_Integer_Vector(lsroco.all);
        begin
          Assign(sv,b);
        end;
        Clear(lsroco);
      end if;
    end if;
   -- Assign(integer32(rc),a);
   -- put("nr = "); put(integer32(nr),1); new_line;
    rcnr(1) := integer32(rc);
    rcnr(2) := integer32(nr);
    Assign(rcnr,a);
    Standard_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.Standard_Polynomial_Solver.");
      end if;
      return 77;
  end Standard_Polynomial_Solver;

  function Standard_Laurent_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32 is

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;
    use Standard_Complex_Laur_Systems;
    use Interfaces.C;

    v_b : constant C_Integer_Array(0..1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc,nr : natural32;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    lsroco : Link_to_String;
    sols : Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Standard_Laurent_Solver");
    end if;
   -- put("nv = "); put(integer32(nv),1);
   -- put("  nq = "); put(integer32(nq),1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 75;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 75;
    end if;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
     -- Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
      if silent then
        if ntasks = 0 -- patch for multitasking and deflation
         then Black_Box_Solvers.Solve(lp.all,silent,rc,sols);
         else Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
        end if;
      else
        if ntasks = 0 -- patch for multitasking and deflation
         then Black_Box_Solvers.Solve(lp.all,rc,lsroco,sols);
         else Black_Box_Solvers.Solve(ntasks,lp.all,rc,lsroco,sols);
        end if;
        if lsroco = null then
          nr := 0;
        else
          nr := natural32(lsroco'last);
          declare
            sv : constant Standard_Integer_Vectors.Vector
               := String_to_Integer_Vector(lsroco.all);
          begin
            Assign(sv,b);
          end;
          Clear(lsroco);
        end if;
      end if;
    else
      declare
        use Standard_Laur_Poly_Convertors;
        p : Poly_Sys(lp'range) := Positive_Laurent_Polynomial_System(lp.all);
      begin
       -- Black_Box_Solvers.Solve(ntasks,p,silent,rc,sols);
        if silent then
          if ntasks = 0 -- patch for deflation with multitasking
           then Black_Box_Solvers.Solve(p,silent,true,rc,sols);
           else Black_Box_Solvers.Solve(ntasks,p,silent,true,rc,sols);
          end if;
        else
          if ntasks = 0 -- patch for deflation with multitasking
           then Black_Box_Solvers.Solve(p,true,rc,lsroco,sols);
           else Black_Box_Solvers.Solve(ntasks,p,true,rc,lsroco,sols);
          end if;
          if lsroco = null then
            nr := 0;
          else
            nr := natural32(lsroco'last);
            declare
              sv : constant Standard_Integer_Vectors.Vector
                 := String_to_Integer_Vector(lsroco.all);
            begin
              Assign(sv,b);
            end;
            Clear(lsroco);
          end if;
        end if;
        Standard_Complex_Poly_Systems.Clear(p);
      end;
    end if;
   -- Assign(integer32(rc),a);
    rcnr(1) := integer32(rc);
    rcnr(2) := integer32(nr);
    Assign(rcnr,a);
    Standard_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.Standard_Laurent_Solver.");
      end if;
      return 75;
  end Standard_Laurent_Solver;

  function DoblDobl_Polynomial_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32 is

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;
    use Interfaces.C;

    v_b : constant C_Integer_Array(0..1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc,nr : natural32;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    lsroco : Link_to_String;
    sols : Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.DoblDobl_Polynomial_Solver");
    end if;
   -- put("Dimension of the system in the container : "); put(n,1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 700;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 700;
    end if;
    if silent then
      if ntasks = 0 -- patch for multitasking and deflation
       then Black_Box_Solvers.Solve(lp.all,silent,rc,sols);
       else Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
      end if;
    else
      if ntasks = 0 -- patch for multitasking and deflation
       then Black_Box_Solvers.Solve(lp.all,rc,lsroco,sols);
       else Black_Box_Solvers.Solve(ntasks,lp.all,rc,lsroco,sols);
      end if;
      if lsroco = null then
        nr := 0;
      else
        nr := natural32(lsroco'last); 
        declare
          sv : constant Standard_Integer_Vectors.Vector
             := String_to_Integer_Vector(lsroco.all);
        begin
          Assign(sv,b);
        end;
        Clear(lsroco);
      end if;
    end if;
    rcnr(1) := integer32(rc);
    rcnr(2) := integer32(nr);
    Assign(rcnr,a);
    DoblDobl_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.DoblDobl_Polynomial_Solver.");
      end if;
      return 700;
  end DoblDobl_Polynomial_Solver;

  function DoblDobl_Laurent_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32 is

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Laur_Systems;
    use Interfaces.C;

    v_b : constant C_Integer_Array(0..1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc,nr : natural32;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    lsroco : Link_to_String;
    sols : Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.DoblDobl_Laurent_Solver");
    end if;
   -- put("nv = "); put(integer32(nv),1);
   -- put("  nq = "); put(integer32(nq),1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 701;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 701;
    end if;
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
      if silent then
        if ntasks = 0 -- patch for multitasking and deflation
         then Black_Box_Solvers.Solve(lp.all,silent,rc,sols);
         else Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
        end if;
      else
        if ntasks = 0 -- patch for multitasking and deflation
         then Black_Box_Solvers.Solve(lp.all,rc,lsroco,sols);
         else Black_Box_Solvers.Solve(ntasks,lp.all,rc,lsroco,sols);
        end if;
        if lsroco = null then
          nr := 0;
        else
          nr := natural32(lsroco'last);
          declare
            sv : constant Standard_Integer_Vectors.Vector
               := String_to_Integer_Vector(lsroco.all);
          begin
            Assign(sv,b);
          end;
          Clear(lsroco);
        end if;
      end if;
    else
      declare
        use DoblDobl_Laur_Poly_Convertors;
        p : constant Poly_Sys := Positive_Laurent_Polynomial_System(lp.all);
      begin
        if silent then
          if ntasks = 0 -- patch for multitasking and deflation
           then Black_Box_Solvers.Solve(p,silent,rc,sols);
           else Black_Box_Solvers.Solve(ntasks,p,silent,rc,sols);
          end if;
        else
          if ntasks = 0 -- patch for multitasking and deflation
           then Black_Box_Solvers.Solve(p,rc,lsroco,sols);
           else Black_Box_Solvers.Solve(ntasks,p,rc,lsroco,sols);
          end if;
          if lsroco = null then
            nr := 0;
          else
            nr := natural32(lsroco'last);
            declare
              sv : constant Standard_Integer_Vectors.Vector
                 := String_to_Integer_Vector(lsroco.all);
            begin
              Assign(sv,b);
            end;
            Clear(lsroco);
          end if;
        end if;
      end;
    end if;
    rcnr(1) := integer32(rc);
    rcnr(2) := integer32(nr);
    Assign(rcnr,a);
    DoblDobl_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.DoblDobl_Laurent_Solver.");
      end if;
      return 701;
  end DoblDobl_Laurent_Solver;

  function QuadDobl_Polynomial_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32 is

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;
    use Interfaces.C;

    v_b : constant C_Integer_Array(0..1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc,nr : natural32;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    lsroco : Link_to_String;
    sols : Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.QuadDobl_Polynomial_Solver");
    end if;
   -- put("Dimension of the system in the container : "); put(n,1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 702;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 702;
    end if;
    if silent then
      if ntasks = 0 -- patch for multitasking and deflation
       then Black_Box_Solvers.Solve(lp.all,silent,rc,sols);
       else Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
      end if;
    else
      if ntasks = 0 -- patch for multitasking and deflation
       then Black_Box_Solvers.Solve(lp.all,rc,lsroco,sols);
       else Black_Box_Solvers.Solve(ntasks,lp.all,rc,lsroco,sols);
      end if;
      if lsroco = null then
        nr := 0;
      else
        nr := natural32(lsroco'last);
        declare
          sv : constant Standard_Integer_Vectors.Vector
             := String_to_Integer_Vector(lsroco.all);
        begin
          Assign(sv,b);
        end;
        Clear(lsroco);
      end if;
    end if;
    rcnr(1) := integer32(rc);
    rcnr(2) := integer32(nr);
    Assign(rcnr,a);
    QuadDobl_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.QuadDobl_Polynomial_Solver.");
      end if;
      return 702;
  end QuadDobl_Polynomial_Solver;

  function QuadDobl_Laurent_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32 is

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Laur_Systems;
    use Interfaces.C;

    v_b : constant C_Integer_Array(0..1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc,nr : natural32;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    lsroco : Link_to_String;
    sols : Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.QuadDobl_Laurent_Solver");
    end if;
   -- put("nv = "); put(integer32(nv),1);
   -- put("  nq = "); put(integer32(nq),1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 703;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 703;
    end if;
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
      if silent then
        if ntasks = 0 -- patch for multitasking and deflation
         then Black_Box_Solvers.Solve(lp.all,silent,rc,sols);
         else Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
        end if;
      else
        if ntasks = 0 -- patch for multitasking and deflation
         then Black_Box_Solvers.Solve(lp.all,rc,lsroco,sols);
         else Black_Box_Solvers.Solve(ntasks,lp.all,rc,lsroco,sols);
        end if;
        if lsroco = null then
          nr := 0;
        else
          nr := natural32(lsroco'last);
          declare
            sv : constant Standard_Integer_Vectors.Vector
               := String_to_Integer_Vector(lsroco.all);
          begin
            Assign(sv,b);
          end;
          Clear(lsroco);
        end if;
      end if;
    else
      declare
        use QuadDobl_Laur_Poly_Convertors;
        p : constant Poly_Sys := Positive_Laurent_Polynomial_System(lp.all);
      begin
        if silent then
          if ntasks = 0 -- patch for multitasking and deflation
           then Black_Box_Solvers.Solve(p,silent,rc,sols);
           else Black_Box_Solvers.Solve(ntasks,p,silent,rc,sols);
          end if;
        else
          if ntasks = 0 -- patch for multitasking and deflation
           then Black_Box_Solvers.Solve(p,rc,lsroco,sols);
           else Black_Box_Solvers.Solve(ntasks,p,rc,lsroco,sols);
          end if;
          if lsroco = null then
            nr := 0;
          else
            nr := natural32(lsroco'last);
            declare
              sv : constant Standard_Integer_Vectors.Vector
                 := String_to_Integer_Vector(lsroco.all);
            begin
              Assign(sv,b);
            end;
            Clear(lsroco);
          end if;
        end if;
      end;
    end if;
    rcnr(1) := integer32(rc);
    rcnr(2) := integer32(nr);
    Assign(rcnr,a);
    QuadDobl_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.QuadDobl_Laurent_Solver.");
      end if;
      return 703;
  end QuadDobl_Laurent_Solver;

  function Mixed_Volume
             ( a : C_intarrs.Pointer; vrblvl : integer := 0 )
             return integer32 is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;
    use Black_Mixed_Volume_Computations;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
    mixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    mv : natural32;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Mixed_Volume");
    end if;
    if lp /= null then
      Black_Box_Mixed_Volume_Computation
        (lp.all,mix,perm,iprm,lifsup,mixsub,mv);
    else
      declare
        lq : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
      begin
        Black_Box_Mixed_Volume_Computation
          (lq.all,mix,perm,iprm,lifsup,mixsub,mv);
      end;
    end if;
    Assign(integer32(mv),a);
    Cells_Container.Initialize(mix,lifsup,mixsub);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.Mixed_Volume.");
      end if;
      return 78;
  end Mixed_Volume;

  function Stable_Mixed_Volume
             ( a,b : C_intarrs.Pointer; vrblvl : integer := 0 )
             return integer32 is

    use Standard_Complex_Poly_Systems;
    use Black_Mixed_Volume_Computations;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
    mixsub,orgmcc,stbmcc : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    mv,smv,totmv,orgcnt,stbcnt : natural32;
    stlb : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Stable_Mixed_Volume");
    end if;
    Black_Box_Mixed_Volume_Computation
      (lp.all,mix,perm,iprm,stlb,lifsup,mixsub,orgmcc,stbmcc,mv,smv,totmv,
       orgcnt,stbcnt);
    Assign(integer32(mv),a);
    Assign(integer32(smv),b);
    Cells_Container.Initialize(stlb,mix,lifsup,mixsub,orgmcc,stbmcc);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.Stable_Mixed_Volume.");
      end if;
      return 79;
  end Stable_Mixed_Volume;

end Job_Handlers;
