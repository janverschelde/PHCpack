with Interfaces.C;
with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Timing_Package;                    use Timing_Package;
with Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;
with Standard_Laur_Poly_Convertors;
with Standard_Poly_Laur_Convertors;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Laur_Poly_Convertors;
with DoblDobl_Poly_Laur_Convertors;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Laur_Poly_Convertors;
with QuadDobl_Poly_Laur_Convertors;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Condition_Tables;
with Standard_Coefficient_Circuits;
with Standard_Circuit_Makers;
with Standard_Refiner_Circuits;
with Floating_Mixed_Subdivisions;
with Black_Mixed_Volume_Computations;
with Black_Box_Solvers;
with Black_Box_Polyhedral_Solvers;
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
with Number_of_Cores;
with PHCpack_Operations;

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

  function Get_Core_Count ( a : C_intarrs.Pointer;
                            vrblvl : integer32 := 0 ) return integer32 is

    nbr : constant integer32 := Number_of_Cores;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Get_Core_Count");
    end if;
    Assign(nbr,a);
    return 0;
  exception
    when others =>
      if vrblvl > 0
       then put_line("Exception raised in job_handlers.Get_Core_Count.");
      end if;
      return 994;
  end Get_Core_Count;

  function Standard_Polynomial_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32 is

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;
    use Interfaces.C;

   -- n : constant natural := Standard_PolySys_Container.Dimension;
    v_b : constant C_Integer_Array(0..2)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(3));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    mvfocus : constant natural32 := natural32(v_b(v_b'first+2));
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc,nr : natural32 := 0;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    sols : Solution_List;
    lsroco : Link_to_String;
    gamma : Standard_Complex_Numbers.Complex_Number
          := Standard_Complex_Numbers.Create(0.0);

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
    declare
      q : Poly_Sys(lp'range);
      qsols : Standard_Complex_Solutions.Solution_List;
    begin
      if silent then
        if mvfocus = 0 then
          if ntasks = 0 then
            Black_Box_Solvers.Solve
              (lp.all,silent,true,rc,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,lp.all,silent,true,rc,gamma,q,qsols,sols,vrblvl-1);
          end if;
        else
          if ntasks = 0 then
            Black_Box_Polyhedral_Solvers.Solve
              (lp.all,silent,true,rc,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Polyhedral_Solvers.Solve
              (ntasks,lp.all,silent,true,rc,gamma,q,qsols,sols,vrblvl-1);
          end if;
        end if;
      else
        if mvfocus = 0 then
          if ntasks = 0 then
            Black_Box_Solvers.Solve
              (lp.all,true,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,lp.all,true,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
          end if;
        else
          if ntasks = 0 then
            Black_Box_Polyhedral_Solvers.Solve
              (lp.all,true,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Polyhedral_Solvers.Solve
              (ntasks,lp.all,true,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
          end if;
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
      PHCpack_Operations.Store_Start_System(q);
      PHCpack_Operations.Store_Start_Solutions(qsols);
      PHCpack_Operations.Store_Gamma_Constant(gamma);
      Standard_Complex_Poly_Systems.Clear(q);
      Standard_Complex_Solutions.Deep_Clear(qsols);
    end;
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

    v_b : constant C_Integer_Array(0..2)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(3));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    mvfocus : constant natural32 := natural32(v_b(v_b'first+2));
    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc,nr : natural32 := 0; -- must be added in case silent
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    lsroco : Link_to_String;
    sols : Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number
          := Standard_Complex_Numbers.Create(0.0);

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
      declare
        q : Laur_Sys(lp'range);
        qsols : Standard_Complex_Solutions.Solution_List;
      begin
       -- Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
        if silent then
          if ntasks = 0 then -- patch for multitasking and deflation
            Black_Box_Solvers.Solve
              (lp.all,silent,rc,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,lp.all,silent,rc,gamma,q,qsols,sols,vrblvl-1);
          end if;
        else
          if ntasks = 0 then -- patch for multitasking and deflation
            Black_Box_Solvers.Solve
              (lp.all,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
             (ntasks,lp.all,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
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
        PHCpack_Operations.Store_Start_System(q);
        PHCpack_Operations.Store_Start_Solutions(qsols);
        PHCpack_Operations.Store_Gamma_Constant(gamma);
        Standard_Complex_Laur_Systems.Clear(q);
        Standard_Complex_Solutions.Deep_Clear(qsols);
      end;
    else
      declare
        use Standard_Laur_Poly_Convertors;
        use Standard_Poly_Laur_Convertors;
        p : Poly_Sys(lp'range) := Positive_Laurent_Polynomial_System(lp.all);
        q : Poly_Sys(p'range);
        qq : Laur_Sys(p'range);
        qsols : Solution_List;
      begin
       -- Black_Box_Solvers.Solve(ntasks,p,silent,rc,sols);
        if silent then
          if mvfocus = 0 then
            if ntasks = 0 then -- patch for deflation with multitasking
              Black_Box_Solvers.Solve
                (p,silent,true,rc,gamma,q,qsols,sols,vrblvl-1);
            else
              Black_Box_Solvers.Solve
                (ntasks,p,silent,true,rc,gamma,q,qsols,sols,vrblvl-1);
            end if;
          else 
            if ntasks = 0 then -- patch for deflation with multitasking
              Black_Box_Polyhedral_Solvers.Solve
                (p,silent,true,rc,gamma,q,qsols,sols,vrblvl-1);
            else
              Black_Box_Polyhedral_Solvers.Solve
                (ntasks,p,silent,true,rc,gamma,q,qsols,sols,vrblvl-1);
            end if;
          end if;
        else
          if mvfocus = 0 then
            if ntasks = 0 then -- patch for deflation with multitasking
              Black_Box_Solvers.Solve
                (p,true,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
            else
              Black_Box_Solvers.Solve
                (ntasks,p,true,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
            end if;
          else
            if ntasks = 0 then -- patch for deflation with multitasking
              Black_Box_Polyhedral_Solvers.Solve
                (p,true,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
            else
              Black_Box_Polyhedral_Solvers.Solve
                (ntasks,p,true,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
            end if;
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
        qq := Polynomial_to_Laurent_System(q);
        PHCpack_Operations.Store_Start_System(qq);
        PHCpack_Operations.Store_Start_Solutions(qsols);
        PHCpack_Operations.Store_Gamma_Constant(gamma);
        Standard_Complex_Poly_Systems.Clear(q);
        Standard_Complex_Laur_Systems.Clear(qq);
        Standard_Complex_Solutions.Deep_Clear(qsols);
      end;
     -- Assign(integer32(rc),a);
    end if;
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
    rc,nr : natural32 := 0;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    lsroco : Link_to_String;
    sols : Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Complex_Numbers.Create(integer(0));

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
    declare
      q : Poly_Sys(lp'range);
      qsols : Solution_List;
    begin
      if silent then
        if ntasks = 0 then -- patch for multitasking and deflation
          Black_Box_Solvers.Solve
            (lp.all,silent,rc,gamma,q,qsols,sols,vrblvl-1);
        else
          Black_Box_Solvers.Solve
            (ntasks,lp.all,silent,rc,gamma,q,qsols,sols,vrblvl-1);
        end if;
      else
        if ntasks = 0 then -- patch for multitasking and deflation
          Black_Box_Solvers.Solve
            (lp.all,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
        else
          Black_Box_Solvers.Solve
            (ntasks,lp.all,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
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
      PHCpack_Operations.Store_Start_System(q);
      PHCpack_Operations.Store_Start_Solutions(qsols);
      PHCpack_Operations.Store_Gamma_Constant(gamma);
      DoblDobl_Complex_Poly_Systems.Clear(q);
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
    end;
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
    rc,nr : natural32 := 0;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    lsroco : Link_to_String;
    sols : Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Complex_Numbers.Create(integer(0));

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
      declare
        q : Laur_Sys(lp'range);
        qsols : Solution_List;
      begin
        if silent then
          if ntasks = 0 then -- patch for multitasking and deflation
            Black_Box_Solvers.Solve
              (lp.all,silent,rc,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,lp.all,silent,rc,gamma,q,qsols,sols,vrblvl-1);
          end if;
        else
          if ntasks = 0 then -- patch for multitasking and deflation
            Black_Box_Solvers.Solve
              (lp.all,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,lp.all,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
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
        PHCpack_Operations.Store_Start_System(q);
        PHCpack_Operations.Store_Start_Solutions(qsols);
        PHCpack_Operations.Store_Gamma_Constant(gamma);
        DoblDobl_Complex_Laur_Systems.Clear(q);
        DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      end;
    else
      declare
        use DoblDobl_Laur_Poly_Convertors;
        use DoblDobl_Poly_Laur_Convertors;
        p : constant Poly_Sys := Positive_Laurent_Polynomial_System(lp.all);
        q : Poly_Sys(p'range);
        qq : Laur_Sys(p'range);
        qsols : Solution_List;
      begin
        if silent then
          if ntasks = 0 then -- patch for multitasking and deflation
            Black_Box_Solvers.Solve(p,silent,rc,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,p,silent,rc,gamma,q,qsols,sols,vrblvl-1);
          end if;
        else
          if ntasks = 0 then -- patch for multitasking and deflation
            Black_Box_Solvers.Solve(p,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,p,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
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
        qq := Polynomial_to_Laurent_System(q);
        PHCpack_Operations.Store_Start_System(qq);
        PHCpack_Operations.Store_Start_Solutions(qsols);
        PHCpack_Operations.Store_Gamma_Constant(gamma);
        DoblDobl_Complex_Poly_Systems.Clear(q);
        DoblDobl_Complex_Laur_Systems.Clear(qq);
        DoblDobl_Complex_Solutions.Deep_Clear(qsols);
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
    rc,nr : natural32 := 0;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    lsroco : Link_to_String;
    sols : Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Complex_Numbers.Create(integer(0));

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
    declare
      q : Poly_Sys(lp'range);
      qsols : Solution_List;
    begin
      if silent then
        if ntasks = 0 then -- patch for multitasking and deflation
          Black_Box_Solvers.Solve
            (lp.all,silent,rc,gamma,q,qsols,sols,vrblvl-1);
        else
          Black_Box_Solvers.Solve
            (ntasks,lp.all,silent,rc,gamma,q,qsols,sols,vrblvl-1);
        end if;
      else
        if ntasks = 0 then -- patch for multitasking and deflation
          Black_Box_Solvers.Solve
            (lp.all,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
        else
          Black_Box_Solvers.Solve
            (ntasks,lp.all,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
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
      PHCpack_Operations.Store_Start_System(q);
      PHCpack_Operations.Store_Start_Solutions(qsols);
      PHCpack_Operations.Store_Gamma_Constant(gamma);
      QuadDobl_Complex_Poly_Systems.Clear(q);
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
    end;
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
    rc,nr : natural32 := 0;
    rcnr : Standard_Integer_Vectors.Vector(1..2);
    lsroco : Link_to_String;
    sols : Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Complex_Numbers.Create(integer(0));

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
      declare
        q : Laur_Sys(lp'range);
        qsols : Solution_List;
      begin
        if silent then
          if ntasks = 0 then -- patch for multitasking and deflation
            Black_Box_Solvers.Solve
              (lp.all,silent,rc,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,lp.all,silent,rc,gamma,q,qsols,sols,vrblvl-1);
          end if;
        else
          if ntasks = 0 then -- patch for multitasking and deflation
            Black_Box_Solvers.Solve
              (lp.all,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,lp.all,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
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
        PHCpack_Operations.Store_Start_System(q);
        PHCpack_Operations.Store_Start_Solutions(qsols);
        PHCpack_Operations.Store_Gamma_Constant(gamma);
        QuadDobl_Complex_Laur_Systems.Clear(q);
        QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      end;
    else
      declare
        use QuadDobl_Laur_Poly_Convertors;
        use QuadDobl_Poly_Laur_Convertors;
        p : constant Poly_Sys := Positive_Laurent_Polynomial_System(lp.all);
        q : Poly_Sys(p'range);
        qq : Laur_Sys(q'range);
        qsols : Solution_List;
      begin
        if silent then
          if ntasks = 0 then -- patch for multitasking and deflation
            Black_Box_Solvers.Solve(p,silent,rc,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,p,silent,rc,gamma,q,qsols,sols,vrblvl-1);
          end if;
        else
          if ntasks = 0 then -- patch for multitasking and deflation
            Black_Box_Solvers.Solve(p,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
          else
            Black_Box_Solvers.Solve
              (ntasks,p,rc,lsroco,gamma,q,qsols,sols,vrblvl-1);
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
        qq := Polynomial_to_Laurent_System(q);
        PHCpack_Operations.Store_Start_System(qq);
        PHCpack_Operations.Store_Start_Solutions(qsols);
        PHCpack_Operations.Store_Gamma_Constant(gamma);
        QuadDobl_Complex_Poly_Systems.Clear(q);
        QuadDobl_Complex_Laur_Systems.Clear(qq);
        QuadDobl_Complex_Solutions.Deep_Clear(qsols);
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

  function Get_Gamma_Constant
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    gamma : Standard_Floating_Vectors.Vector(1..2);
    regamma,imgamma : double_float;
    ddregamma,ddimgamma : double_double;
    qdregamma,qdimgamma : quad_double;
    dgamma : Standard_Complex_Numbers.Complex_Number;
    ddgamma : DoblDobl_Complex_Numbers.Complex_Number;
    qdgamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Get_Gamma_Constant ...");
    end if;
    case prc is
      when 1 => 
        PHCpack_Operations.Retrieve_Gamma_Constant(dgamma);
        regamma := Standard_Complex_Numbers.REAL_PART(dgamma);
        imgamma := Standard_Complex_Numbers.IMAG_PART(dgamma);
      when 2 => 
        PHCpack_Operations.Retrieve_Gamma_Constant(ddgamma);
        ddregamma := DoblDobl_Complex_Numbers.REAL_PART(ddgamma);
        ddimgamma := DoblDobl_Complex_Numbers.IMAG_PART(ddgamma);
        regamma := hi_part(ddregamma);
        imgamma := hi_part(ddimgamma);
      when 4 => 
        PHCpack_Operations.Retrieve_Gamma_Constant(qdgamma);
        qdregamma := QuadDobl_Complex_Numbers.REAL_PART(qdgamma);
        qdimgamma := QuadDobl_Complex_Numbers.IMAG_PART(qdgamma);
        regamma := hihi_part(qdregamma);
        imgamma := hihi_part(qdimgamma);
      when others =>
        if vrblvl > 0 then
          put("Value for precision "); put(prc,1);
          put_line(" is not valid.");
        end if;
    end case;
    gamma(1) := regamma;
    gamma(2) := imgamma;
    Assign(gamma,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.Get_Gamma_Constant.");
      end if;
      return 995;
  end Get_Gamma_Constant;

  function Set_Gamma_Constant
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    v_c : constant C_Double_Array
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(2));
    regamma : constant double_float := double_float(v_c(v_c'first));
    imgamma : constant double_float := double_float(v_c(v_c'first+1));
    ddregamma,ddimgamma : double_double;
    qdregamma,qdimgamma : quad_double;
    dgamma : Standard_Complex_Numbers.Complex_Number;
    ddgamma : DoblDobl_Complex_Numbers.Complex_Number;
    qdgamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Set_Gamma_Constant ...");
    end if;
    case prc is
      when 1 => 
        dgamma := Standard_Complex_Numbers.Create(regamma,imgamma);
        PHCpack_Operations.Store_Gamma_Constant(dgamma);
      when 2 =>
        ddregamma := create(regamma); ddimgamma := create(imgamma);
        ddgamma := DoblDobl_Complex_Numbers.Create(ddregamma,ddimgamma);
        PHCpack_Operations.Store_Gamma_Constant(ddgamma);
      when 4 =>
        qdregamma := create(regamma); qdimgamma := create(imgamma);
        qdgamma := QuadDobl_Complex_Numbers.Create(qdregamma,qdimgamma);
        PHCpack_Operations.Store_Gamma_Constant(qdgamma);
      when others =>
        if vrblvl > 0 then
          put("Value for precision "); put(prc,1);
          put_line(" is not valid.");
        end if;
    end case;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.Set_Gamma_Constant.");
      end if;
      return 996;
  end Set_Gamma_Constant;

  function Mixed_Volume
             ( a : C_intarrs.Pointer; vrblvl : integer32 := 0 )
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
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
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

  function Standard_Condition_Report
             ( a,b : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    maxit : constant natural32 := natural32(v_a(v_a'first));
    vrbval : constant natural32 := natural32(v_a(v_a'first+1));
    verbose : constant boolean := (vrbval = 1);
    nbrchar : constant natural32 := natural32(v_a(v_a'first+2));
    v_c : constant C_Double_Array
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(3));
    tolres : constant double_float := double_float(v_c(v_c'first));
    tolerr : constant double_float := double_float(v_c(v_c'first+1));
    tolsing : constant double_float := double_float(v_c(v_c'first+2));
    cnt : Standard_Natural_Vectors.Vector(1..6) := (1..6 => 0);
    tab : Standard_Natural_Vectors.Vector(1..48) := (1..48 => 0);
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15);
    p : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Standard_PolySys_Container.Retrieve;
    sols : constant Standard_Complex_Solutions.Solution_List 
         := Standard_Solutions_Container.Retrieve;
    s : Standard_Coefficient_Circuits.Link_to_System
      := Standard_Circuit_Makers.Make_Coefficient_System(p,false);
    timer : Timing_Widget;
    idx : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in job_handlers.Standard_Condition_Report");
    end if;
    if verbose then
      new_line;
      Standard_Refiner_Circuits.Show_Parameters(maxit,tolres,tolerr,tolsing);
    end if;
    if nbrchar = 0 then
      Standard_Refiner_Circuits.Inlined_Run
        (s,sols,maxit,tolres,tolerr,tolsing,cnt(1),cnt(2),cnt(3),
         cnt(4),cnt(5),cnt(6),t_err,t_rco,t_res,verbose,vrblvl);
      if verbose then
        Standard_Condition_Tables.Write_Tables
          (standard_output,t_err,t_res,t_rco);
      end if;
    else
      declare
        n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nbrchar-1);
        v_b : constant C_Integer_Array(0..n1)
            := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbrchar));
        name : constant String(1..integer(nbrchar))
             := C_Integer_Array_to_String(nbrchar,v_b);
        file : file_type;
      begin
        if verbose then
          new_line;
          put_line("Creating the file " & name & " ...");
        end if;
        Communications_with_User.Create_Output_File(file,name);
        put(file,natural32(p'last),p.all);
        new_line(file);
        tstart(timer);
        Standard_Refiner_Circuits.Inlined_Run
          (file,s,sols,maxit,tolres,tolerr,tolsing,cnt(1),cnt(2),cnt(3),
           cnt(4),cnt(5),cnt(6),t_err,t_rco,t_res,verbose,vrblvl);
        tstop(timer);
        put(file,"number of regular solutions   : ");
        put(file,cnt(4),1); new_line(file);
        put(file,"number of singular solutions  : ");
        put(file,cnt(5),1); new_line(file);
        put(file,"number of real solutions      : ");
        put(file,cnt(2),1); new_line(file);
        put(file,"number of complex solutions   : ");
        put(file,cnt(3),1); new_line(file);
        put(file,"number of clustered solutions : ");
        put(file,cnt(6),1); new_line(file);
        put(file,"number of failures            : ");
        put(file,cnt(1),1); new_line(file);
        Standard_Complex_Solutions_io.put_bar(file);
        Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
        new_line(file);
        print_times(file,timer,"Newton with condition table report");
        Standard_Coefficient_Circuits.Clear(s);
        close(file);
      end;
    end if;
    idx := 0;
    for i in t_err'range loop
      idx := idx + 1; tab(idx) := t_err(i);
    end loop;
    for i in t_err'range loop
      idx := idx + 1; tab(idx) := t_rco(i);
    end loop;
    for i in t_err'range loop
      idx := idx + 1; tab(idx) := t_res(i);
    end loop;
    Assign(cnt,a);
    Assign(tab,b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised");
        put_line(" in job_handlers.Standard_Condition_Report.");
      end if;
      return 920;
  end Standard_Condition_Report;

end Job_Handlers;
