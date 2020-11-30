with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_VecMats_io;        use Standard_Complex_VecMats_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
--with Standard_Matrix_Inversion;          use Standard_Matrix_Inversion;
with Symbol_Table;                       use Symbol_Table;
with Matrix_Indeterminates;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Homotopy;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Brackets,Brackets_io;               use Brackets,Brackets_io;
with Bracket_Monomials;                  use Bracket_Monomials;
with Standard_Bracket_Systems;           use Standard_Bracket_Systems;
with Specialization_of_Planes;           use Specialization_of_Planes;
with Symbolic_Minor_Equations;           use Symbolic_Minor_Equations;
with Numeric_Minor_Equations;            use Numeric_Minor_Equations;
with Determinantal_Systems;              use Determinantal_Systems;
with Plane_Representations;              use Plane_Representations;
with Localization_Posets;                use Localization_Posets;
with Localization_Posets_io;             use Localization_Posets_io;
with Deformation_Posets;                 use Deformation_Posets;
with Make_Input_Planes;

package body Main_Pieri_Homotopies is

-- AUXILIARIES IN THE SOLVERS :

  procedure Add_t_Symbol is

  -- DESCRIPTION :
  --   Adds the symbol for the continuation parameter t to the symbol table.

    tsb : Symbol;

  begin
    Symbol_Table.Enlarge(1);
    tsb(1) := 't';
    for i in 2..tsb'last loop
      tsb(i) := ' ';
    end loop;
    Symbol_Table.Add(tsb);
  end Add_t_Symbol;

  procedure Set_Parameters ( file : in file_type; report : out boolean ) is

  -- DESCRIPTION :
  --   Interactive determination of the continuation and output parameters.

    oc : natural32;

  begin
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := not (oc = 0);
    new_line;
    put_line("No more input expected.  See output file for results...");
    new_line;
    new_line(file);
  end Set_Parameters;

  function First_Standard_Plane
             ( m,p : natural32 ) return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the m-plane spanned by the first m standard basis vectors.

    res : Standard_Complex_Matrices.Matrix(1..integer32(m+p),1..integer32(m));

  begin
    for j in res'range(2)loop          -- assign j-th column
      for i in res'range(1) loop
        res(i,j) := Create(0.0);       -- initialize with zeros
      end loop;
      res(j,j) := Create(1.0);         -- j-th column = j-th basis vector
    end loop;
    return res;
  end First_Standard_Plane;

  function Last_Standard_Plane
             ( m,p : natural32 ) return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the m-plane spanned by the last m standard basis vectors.

    res : Standard_Complex_Matrices.Matrix(1..integer32(m+p),1..integer32(m));

  begin
    for j in res'range(2) loop          -- assign j-th column
      for i in res'range(1) loop
        res(i,j) := Create(0.0);        -- initialize with zeros
      end loop;                         -- j-th vector = (p+j)-th basis vector
      res(integer32(p)+j,j) := Create(1.0); 
    end loop;
    return res;
  end Last_Standard_Plane;

--  procedure Basis_Change
--              ( n : in natural; vm : in VecMat; transvm : out VecMat;
--                trans : out Standard_Complex_Matrices.Matrix ) is
--
--  -- DESCRIPTION :
--  --   Changes basis so that the last two planes in vm are spanned by
--  --   the first and last standard basis vectors respectively.
--
--    wrk : Standard_Complex_Matrices.Matrix(1..n,1..n);
--   -- penuplane : constant Standard_Complex_Matrices.Matrix
--   --   := vm(vm'last-1).all;
--    frstplane : constant Standard_Complex_Matrices.Matrix := vm(vm'first).all;
--    lastplane : constant Standard_Complex_Matrices.Matrix := vm(vm'last).all;
--    colind : natural := 0;
--    use Standard_Complex_Matrices;
--    ranflt : double_float;
--
--  begin
--    for j in frstplane'range(2) loop              -- fill up the work matrix
--      colind := colind + 1;
--      for i in frstplane'range(1) loop
--         wrk(i,colind) := frstplane(i,j);
--      end loop;
--    end loop;
--    for j in colind+1..(n-lastplane'length(2)) loop  -- random spacing
--      for i in 1..n loop
--        ranflt := Random;
--        wrk(i,j) := Create(ranflt);
--      end loop;
--    end loop;
--    colind := n+1;
--    for j in reverse lastplane'range(2) loop
--      colind := colind - 1;
--      for i in lastplane'range(1) loop
--         wrk(i,colind) := lastplane(i,j);
--      end loop;
--    end loop;
--    trans := Inverse(wrk);                        -- transformation = inverse
--    for i in vm'first+1..vm'last-1 loop
--      transvm(i-1) := new Standard_Complex_Matrices.Matrix'
--                           (trans*vm(i).all);
--    end loop;
--    transvm(transvm'last-1)
--      := new Standard_Complex_Matrices.Matrix'(trans*vm(vm'first).all);
--    transvm(transvm'last)
--      := new Standard_Complex_Matrices.Matrix'(trans*vm(vm'last).all);
--  end Basis_Change;

  function Solution_Plane ( locmap : Standard_Natural_Matrices.Matrix;
                            mat : Standard_Complex_Matrices.Matrix )
                          return Solution is

  -- DESCRIPTION :
  --   Returns the representation of the solution plane as a vector.

    solloc : constant Standard_Complex_Vectors.Vector
           := Vector_Rep(locmap,Localize(locmap,mat));
    sol : Solution(solloc'length);

  begin
    sol.m := 1;
    sol.t := Create(0.0);
    sol.res := 0.0;
    sol.err := 0.0;
    sol.rco := 0.0;
    sol.v := solloc;
    return sol;
  end Solution_Plane;

  function Solution_Planes ( locmap : Standard_Natural_Matrices.Matrix;
                             vm : VecMat ) return Solution_List is

  -- DESCRIPTION :
  --   Returns the representation of the vector of planes as a solution list.

    res,res_last : Solution_List;

  begin
    for i in vm'range loop
      Append(res,res_last,Solution_Plane(locmap,vm(i).all));
    end loop;
    return res;
  end Solution_Planes;

  procedure Write_Solution_Planes
                ( file : in file_type; sols : in Solution_List;
                  locmap : in Standard_Natural_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the solution vectors as matrices, according to the 
  --   patterns of ones and zeros in the localization map.

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      put(file,Matrix_Rep(locmap,ls.v));
      new_line(file);
      tmp := Tail_Of(tmp);
    end loop;   
  end Write_Solution_Planes;

  function Create_Hypersurface_System
             ( n : natural32; locmap : Standard_Natural_Matrices.Matrix;
               xpm : Standard_Complex_Poly_Matrices.Matrix; planes : VecMat )
             return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the polynomial system that collects the hypersurface
  --   intersection conditions for meeting the given m-planes.
  --   The system is localized according to the given localization map.

    wrksys : Poly_Sys(1..integer32(n))
           := Polynomial_Equations(planes(1..integer32(n)),xpm);
    res : constant Poly_Sys(1..integer32(n)) := Localize(locmap,wrksys);

  begin
    Clear(wrksys);
    return res;
  end Create_Hypersurface_System;

  function Create_General_System
             ( locmap : Standard_Natural_Matrices.Matrix;
               xpm : Standard_Complex_Poly_Matrices.Matrix; planes : VecMat )
             return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the polynomial system that collects the general
  --   intersection conditions for meeting the given (m+1-k(i))-planes.
  --   The system is localized according to the given localization map.

    wrksys : constant Poly_Sys := Polynomial_Equations(planes,xpm);
    res : constant Poly_Sys := Localize(locmap,wrksys);
    tmp : Poly;

  begin
    for i in wrksys'range loop
      tmp := wrksys(i);
      Clear(tmp);
    end loop;
    return res;
  end Create_General_System;

  procedure Evaluate_Roots ( file : in file_type;
                             p : in Poly_sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates the roots at the system, good for overdetermined systems.

    tmp : Solution_List := sols;

  begin
    for i in 1..Length_Of(sols) loop
      put(file,"Solution "); put(file,i,1);
      put_line(file," evaluated at the system : ");
      put_line(file,Eval(p,Head_Of(tmp).v));
      tmp := Tail_Of(tmp);
    end loop;
  end Evaluate_Roots;

  function Square ( n : natural32; p : Poly_Sys ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns a n-by-n system, by adding random linear combinations of 
  --   the polynomials p(i), i > n, to the first n polynomials.

    res : Poly_Sys(1..integer32(n));
    acc : Poly;

  begin
    for i in res'range loop
      Copy(p(i),res(i));
      for j in integer32(n)+1..p'last loop
        acc := Random1*p(j);
        Add(res(i),acc);
        Clear(acc);
      end loop;
    end loop;
    return res;
  end Square;

  procedure Refine_Roots ( file : in file_type;
                           p : in Poly_Sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Calls the root refiner, for square polynomial systems only.

    epsxa : constant double_float := 10.0**(-8);
    epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural32 := 3;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,max,deflate,false);
  end Refine_Roots;

  function General_Homotopy
              ( xpm : Standard_Complex_Poly_Matrices.Matrix;
                locmap : Standard_Natural_Matrices.Matrix;
                start,target : VecMat ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns a general homotopy between start and target planes.
  --   The continuation parameter is inside the determinants.

  -- REQUIRED : dimensions of start and target matrices must correspond.

    res : Link_to_Poly_Sys;
    n : constant natural32 := natural32(xpm'length(1));
    p : constant natural32 := natural32(xpm'length(2));
    nva : constant integer32 := integer32(n*p + 1);

  begin
    for i in start'range loop
      declare
        moving : Standard_Complex_Poly_Matrices.Matrix
                   (start(i)'range(1),start(i)'range(2))
               := Moving_U_Matrix(nva,start(i).all,target(i).all);
        kd : constant integer32 := integer32(p) + start(i)'length(2);
        bm : constant Bracket_Monomial := Maximal_Minors(n,natural32(kd));
        bs : constant Bracket_System(0..integer32(Number_of_Brackets(bm)))
           := Minor_Equations(natural32(kd),natural32(kd)-p,bm); 
        sys : constant Poly_Sys(1..bs'last)
            := Lifted_Expanded_Minors(moving,xpm,bs);
      begin
        Concat(res,sys);
        Standard_Complex_Poly_Matrices.Clear(moving);
      end;
    end loop;
    declare
      locres : constant Poly_Sys := Localize(locmap,res.all);
    begin
      Clear(res);
      return locres;
    end;
  end General_Homotopy;

  procedure Continuation ( file : in file_type; sols : in out Solution_List;
                           report : in boolean ) is

  -- DESCRIPTION :
  --   Calls the path trackers starting at the given solutions.

  -- REQUIRED : Homotopy is properly created.

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                             Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Set_Continuation_Parameter(sols,Create(0.0));
    if report
     then Rep_Cont(file,sols,false,target=>Create(1.0));
     else Sil_Cont(sols,false,target=>Create(1.0));
    end if;
  end Continuation;

  procedure General_Path_Following
               ( file : in file_type; target,homsys : in Poly_Sys;
                 sols : in out Solution_List; report : in boolean ) is

  -- DESCRIPTION :
  --   Given a homotopy with start solutions, the path trackers will
  --   trace the paths defined by this homotopy.

  -- REQUIRED : not Is_Null(sols).

    timer : Timing_Widget;
    dim : constant integer32 := Head_Of(sols).n;   -- actual dimension
    squhom : constant Poly_Sys(1..dim) := Square(natural32(dim),homsys);

  begin
    tstart(timer);
    Standard_Homotopy.Create(squhom,dim+1);
    Continuation(file,sols,report);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Determinantal Cheater's homotopy to target system");
    new_line(file);
    Evaluate_Roots(file,target,sols);
  end General_Path_Following;

  procedure Path_Following_with_Cheater
               ( file : in file_type; start,target : in Poly_Sys;
                 sols : in out Solution_List; report : in boolean ) is

  -- DESCRIPTION :
  --   Calls the standard continuation routines to solve a specific
  --   target system, starting at the solutions of the start system.
  --   This is the usual linear cheater between start and target system,
  --   although it will take nonsquare inputs, it should only be used
  --   for hypersurface intersection conditions.

  -- REQUIRED : not Is_Null(sols).

    timer : Timing_Widget;
    n : constant integer32 := target'last;
    a : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    b : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    dimsol : constant integer32 := Head_Of(sols).n;
    squ_target,squ_start : Poly_Sys(1..dimsol);

  begin
    tstart(timer);
    if dimsol = n then
      Standard_Homotopy.Create(target,start,2,a,b,true);  -- linear cheater
    else
      squ_target := Square(natural32(dimsol),target);
      squ_start  := Square(natural32(dimsol),start);
      Standard_Homotopy.Create(squ_target,squ_start,2,a,b,true);
    end if;
    Continuation(file,sols,report);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Linear Cheater's homotopy to target system");
    if dimsol = n then
      Refine_Roots(file,target,sols);
    else
      Refine_Roots(file,squ_target,sols);
      Evaluate_Roots(file,target,sols);
    end if;
  end Path_Following_with_Cheater;

  procedure Solve_Hypersurface_Target_System
                     ( file : in file_type; m,p : in natural32;
                       start_planes,target_planes : in VecMat;
                       index_poset : in Array_of_Array_of_Nodes;
                       deform_poset : in Array_of_Array_of_VecMats;
                       target_level : in integer32; report : in boolean ) is

  -- DESCRIPTION :
  --   This procedure tests the output of the deformation poset,
  --   creating polynomial representations of the intersection conditions
  --   and solution lists representing the solution planes.

  -- ON ENTRY :
  --   file            to write intermediate output on;
  --   m               dimension of the input planes;
  --   p               dimension of the solution planes;
  --   start_planes    input m-planes in general position;
  --   target_planes   specific input m-planes;
  --   index_poset     indexed localization poset; 
  --   deform_poset    poset with the solution p-planes;
  --   target_level    indicates lowest node in deform_poset that is filled;
  --   report          switch for intermediate output during continuation.

    dim : constant integer32 := integer32(m*p);
    xpm : constant Standard_Complex_Poly_Matrices.Matrix
                     (1..integer32(m+p),1..integer32(p))
        := Localization_Pattern(m+p,index_poset(dim)(1).top,
                                    index_poset(dim)(1).bottom);
    solplanes : constant VecMat := deform_poset(target_level)(1).all;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..integer32(m+p),1..integer32(p))
           := Standard_Coordinate_Frame(xpm,solplanes(1).all);
    locsys : constant Poly_Sys(start_planes'range)
           := Create_Hypersurface_System
                (natural32(dim),locmap,xpm,start_planes);
    target : constant Poly_Sys(target_planes'range)
           := Create_Hypersurface_System
                (natural32(dim),locmap,xpm,target_planes);
    sols : Solution_List := Solution_Planes(locmap,solplanes);

  begin
    Matrix_Indeterminates.Reduce_Symbols(locmap);
    put_line(file,"THE GENERIC SYSTEM : ");
    put_line(file,locsys);
    new_line(file);
    put_line(file,"The localization map :"); put(file,locmap);
    Refine_Roots(file,locsys,sols);
    new_line(file);
    put_line(file,"THE TARGET SYSTEM : ");
    put_line(file,target);
    new_line(file);
    Path_Following_with_Cheater(file,locsys,target,sols,report);
    new_line(file);
    put_line(file,"THE SOLUTION PLANES : ");
    Write_Solution_Planes(file,sols,locmap);
  end Solve_Hypersurface_Target_System;

  procedure Solve_General_Target_System
                     ( file : in file_type; m,p : in natural32;
                       start_planes,target_planes : in VecMat;
                       index_poset : in Array_of_Array_of_Nodes;
                       deform_poset : in Array_of_Array_of_VecMats;
                       target_level : in integer32; report : in boolean ) is

  -- DESCRIPTION :
  --   This procedure tests the output of the deformation poset,
  --   creating polynomial representations of the intersection conditions
  --   and solution lists representing the solution planes.

  -- ON ENTRY :
  --   file            to write intermediate output on;
  --   m               dimension of the input planes;
  --   p               dimension of the solution planes;
  --   start_planes    input (m+1-k(i))-planes in general position;
  --   target_planes   specific input (m+1-k(i))-planes;
  --   index_poset     indexed localization poset; 
  --   deform_poset    poset with the solution p-planes;
  --   target_level    indicates lowest node in deform_poset that is filled;
  --   report          switch for intermediate output during continuation.

    dim : constant integer32 := integer32(m*p);
    xpm : constant Standard_Complex_Poly_Matrices.Matrix
                     (1..integer32(m+p),1..integer32(p))
        := Localization_Pattern(m+p,index_poset(dim)(1).top,
                                    index_poset(dim)(1).bottom);
    solplanes : constant VecMat := deform_poset(target_level)(1).all;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..integer32(m+p),1..integer32(p))
           := Standard_Coordinate_Frame(xpm,solplanes(1).all);
    homtpy : constant Poly_Sys
           := General_Homotopy(xpm,locmap,start_planes,target_planes);
    locsys : constant Poly_Sys
           := Create_General_System(locmap,xpm,start_planes);
    target : constant Poly_Sys
           := Create_General_System(locmap,xpm,target_planes);
    sols : Solution_List := Solution_Planes(locmap,solplanes);

  begin
    Matrix_Indeterminates.Reduce_Symbols(locmap);
    put_line(file,"THE GENERIC SYSTEM : ");
    put_line(file,locsys);
    new_line(file);
    put_line(file,"The localization map :"); put(file,locmap);
    Evaluate_Roots(file,locsys,sols);
    new_line(file);
    put_line(file,"THE TARGET SYSTEM : ");
    put_line(file,target);
    new_line(file);
    General_Path_Following(file,target,homtpy,sols,report);
    new_line(file);
    put_line(file,"THE SOLUTION PLANES : ");
    Write_Solution_Planes(file,sols,locmap);
  end Solve_General_Target_System;

  procedure Write_Poset_Times
               ( file : in file_type; timer : in Timing_Widget;
                 npaths : in Standard_Natural_Vectors.Vector;
                 timings : in Duration_Array ) is

  -- DESCRIPTION :
  --   Writes a overview of #paths and timings spent during deformations
  --   along the poset structure.

  begin
    new_line(file);
    put_line(file,"--------------------------------------");
    put_line(file,"|    TIMING INFORMATION OVERVIEW     |");
    put_line(file,"--------------------------------------");
    put_line(file,"|   n   |  #paths  |  user cpu time  |");
    put_line(file,"--------------------------------------");
    for i in npaths'range loop
      if npaths(i) /= 0 then
        put(file,"|"); put(file,i,4); put(file,"   |");
        put(file,npaths(i),7);        put(file,"   | ");
        print_hms(file,timings(i));   put(file,"  |"); new_line(file);
      end if;
    end loop;
    put_line(file,"--------------------------------------");
    put(file,"| total |");
    put(file,Standard_Natural_Vectors.Sum(npaths),7);
    put(file,"   | ");
    print_hms(file,Elapsed_User_Time(timer));
    put(file,"  |"); new_line(file);
    put_line(file,"--------------------------------------");
  end Write_Poset_Times;

  procedure Solve_Hypersurface_Deformation_Poset
              ( file : in file_type; m,p : in natural32;
                level_poset : in Array_of_Nodes;
                index_poset : in Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates a deformation poset and applies the Solve operator.

    deform_poset : Array_of_Array_of_VecMats(index_poset'range)
                 := Create(index_poset);
    mp : constant integer32 := integer32(m*p);
    planes : constant VecMat(1..mp)
           := Make_Input_Planes.Random_Complex_Planes(m,p);
    target_planes : VecMat(1..mp);
    ans : character;
    report,outlog : boolean;
    timer : Timing_Widget;
    target_level : integer32 := mp;
    nbp : natural32 := 0;
    npaths : Standard_Natural_Vectors.Vector(1..target_level)
           := (1..target_level => 0);
    timings : Duration_Array(1..target_level) := (1..target_level => 0.0);
    nocheater : boolean;

  begin
    put_line("The size of the deformation poset : ");
    put_line(file,"The size of the deformation poset : ");
    put_roco(index_poset);
    put_roco(file,index_poset);
    new_line;
    put("Give target level <= "); put(target_level,1);
    put(" = root level by default : "); get(target_level);
    for i in 1..target_level loop
      nbp := nbp + Row_Root_Count_Sum(level_poset,natural32(i));
    end loop;
    put("The number of paths : "); put(nbp,1); new_line;
    put(file,"The number of paths : "); put(file,nbp,1); new_line(file);
    new_line;
    put("Do you want to have the homotopies on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    outlog := (ans = 'y');
    if target_level < mp then
      planes(mp).all := Last_Standard_Plane(m,p);
      if target_level < mp-1
       then planes(mp-1).all := First_Standard_Plane(m,p);
      end if;
    end if;
    Make_Input_Planes.Main(file,m,p,target_planes,nocheater); 
    Matrix_Indeterminates.Initialize_Symbols(m+p,p);
    Add_t_Symbol;
    Set_Parameters(file,report);
    tstart(timer);
    for i in index_poset(target_level)'range loop
      declare
        root : constant Node := index_poset(target_level)(i).all;
      begin
        Solve(file,m+p,deform_poset,root,planes(1..target_level),
              report,outlog,npaths,timings);
      end;
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Solving along the deformation poset");
    Write_Poset_Times(file,timer,
                      npaths(1..target_level),timings(1..target_level));
    if nocheater then
      new_line(file);
      put_line(file,"THE INPUT PLANES : ");
      put(file,planes);
      new_line(file);
      put_line(file,"THE SOLUTION PLANES : ");
      declare
        solplanes : constant VecMat := deform_poset(target_level)(1).all;
      begin
        put(file,solplanes);
      end;
     -- the statement below seems equivalent to the block above
     -- Write_Solution_Planes(file,sols,locmap);
    else
      if deform_poset(target_level)(1) /= null then
        new_line(file);
        Solve_Hypersurface_Target_System
          (file,m,p,planes,target_planes,index_poset,deform_poset,
           target_level,report);
      end if;
    end if;
  end Solve_Hypersurface_Deformation_Poset;

  procedure Solve_General_Deformation_Poset
              ( file : in file_type; m,p : in natural32; k : in Bracket;
                index_poset : in Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Applies the solver to general intersection conditions.

    deform_poset : Array_of_Array_of_VecMats(index_poset'range)
                 := Create(index_poset);
    planes : constant VecMat(k'range)
           := Make_Input_Planes.Random_Complex_Planes(m,p,k);
    target_planes : VecMat(k'range);
    ans : character;
    report,outlog : boolean;
    timer : Timing_Widget;
    target_level : integer32 := integer32(m*p);
    npaths : Standard_Natural_Vectors.Vector(1..target_level)
           := (1..target_level => 0);
    timings : Duration_Array(1..target_level) := (1..target_level => 0.0);
    nocheater : boolean;

  begin
    put_line("The size of the deformation poset : ");
    put_line(file,"The size of the deformation poset : ");
    put_roco(index_poset);
    put_roco(file,index_poset);
    new_line;
    put("Give target level <= "); put(target_level,1);
    put(" = root level by default : "); get(target_level);
    new_line;
    put("Do you want to have the homotopies on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    outlog := (ans = 'y');
    Make_Input_Planes.Main(file,m,p,k,target_planes,nocheater);
    Matrix_Indeterminates.Initialize_Symbols(m+p,p);
    Add_t_Symbol;
    Set_Parameters(file,report);
    tstart(timer);
    for i in index_poset(target_level)'range loop
      declare
        root : constant Node := index_poset(target_level)(i).all;
      begin
        Solve(file,m+p,k,deform_poset,root,planes,report,outlog,
              npaths,timings);
      end;
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Solving along the deformation poset");
    Write_Poset_Times(file,timer,
                      npaths(1..target_level),timings(1..target_level));
    if nocheater then
      new_line(file);
      put_line(file,"THE INPUT PLANES :");
      put(file,planes);
      new_line(file);
      put_line(file,"THE SOLUTION PLANES : ");
      declare
        solplanes : constant VecMat := deform_poset(target_level)(1).all;
      begin
        put(file,solplanes);
      end;
    else
      if deform_poset(target_level)(1) /= null then
        new_line(file);
        Solve_General_Target_System
          (file,m,p,planes,target_planes,index_poset,deform_poset,
           target_level,report);
      end if;
    end if;
  end Solve_General_Deformation_Poset;

-- HYPERSURFACE SOLVERS :

  procedure Top_Hypersurface_Poset
              ( m,p : in natural32; level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing only top pivots.
  --   The range of the posets on return is 0..m*p.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    Top_Create(lnkroot,m+p);
    put_line("The poset created from the top : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
  end Top_Hypersurface_Poset;

  procedure Top_Hypersurface_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Create the poset by incrementing only top pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    new_line(file);
    tstart(timer);
    Top_Create(lnkroot,m+p);
    put_line("The poset created from the top : ");
    put_line(file,"The poset created from the top : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Hypersurface_Deformation_Poset(file,m,p,level_poset,index_poset);
    tstop(timer);
    new_line(file);
    print_times(file,timer,
                "Total time for Hypersurface Pieri Homotopy Algorithm");
  end Top_Hypersurface_Poset;

  procedure Bottom_Hypersurface_Poset
              ( m,p : in natural32; level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Create the poset by decrementing only bottom pivots.
  --   The range of the posets on return is 0..m*p.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    Bottom_Create(lnkroot);
    put_line("The poset created from the bottom : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
  end Bottom_Hypersurface_Poset;

  procedure Bottom_Hypersurface_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Create the poset by decrementing only bottom pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    new_line(file);
    tstart(timer);
    Bottom_Create(lnkroot);
    put_line("The poset created from the bottom : ");
    put_line(file,"The poset created from the bottom : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Hypersurface_Deformation_Poset(file,m,p,level_poset,index_poset);
    tstop(timer);
    new_line(file);
    print_times(file,timer,
                "Total time for Hypersurface Pieri Homotopy Algorithm");
  end Bottom_Hypersurface_Poset;

  procedure Mixed_Hypersurface_Poset
              ( m,p : in natural32; level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Create the poset by incrementing top and decrementing bottom pivots.
  --   The posets on return are of range 0..m*p.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    Top_Bottom_Create(lnkroot,m+p);
    put_line("The poset created in a mixed fashion : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
  end Mixed_Hypersurface_Poset;

  procedure Mixed_Hypersurface_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Create the poset by incrementing top and decrementing bottom pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    new_line(file);
    tstart(timer);
    Top_Bottom_Create(lnkroot,m+p);
    put_line("The poset created in a mixed fashion : ");
    put_line(file,"The poset created in a mixed fashion : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_Hypersurface_Deformation_Poset(file,m,p,level_poset,index_poset);
    tstop(timer);
    new_line(file);
    print_times(file,timer,
                "Total time for Hypersurface Pieri Homotopy Algorithm");
  end Mixed_Hypersurface_Poset;

-- GENERAL SOLVERS :

  procedure Top_General_Poset
              ( m,p : in natural32; k : out Link_to_Bracket;
                level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

  -- ON ENTRY :
  --   m        dimension of the input planes;
  --   p        dimension of the output planes.
  
  -- ON RETURN :
  --   k        codimension conditions, k(i) = 1 : hypersurface condition;
  --   level_poset is poset of range 0..m*p;
  --   index_poset is poset of range 0..m*p.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);

  begin
    Top_Create(lnkroot,codim,m+p);
    put_line("The poset created from the top : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    k := new Bracket'(codim);
  end Top_General_Poset;

  procedure Top_General_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    put(file,"  k = "); put(file,codim); new_line(file);
    tstart(timer);
    Top_Create(lnkroot,codim,m+p);
    put_line("The poset created from the top : ");
    put_line(file,"The poset created from the top : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
	put(file,index_poset);                           
    Solve_General_Deformation_Poset(file,m,p,codim,index_poset);
    tstop(timer);
    new_line(file);
    print_times(file,timer,
                "Total time for General Pieri Homotopy Algorithm");
  end Top_General_Poset;

  procedure Bottom_General_Poset
              ( m,p : in natural32; k : out Link_to_Bracket;
                level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the bottom pivots.

  -- ON ENTRY :
  --   m        dimension of the input planes;
  --   p        dimension of the output planes.
  
  -- ON RETURN :
  --   k        codimension conditions, k(i) = 1 : hypersurface condition;
  --   level_poset is poset of range 0..m*p;
  --   index_poset is poset of range 0..m*p.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);

  begin
    Bottom_Create(lnkroot,codim);
    put_line("The poset created from the bottom : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    k := new Bracket'(codim);
  end Bottom_General_Poset;

  procedure Bottom_General_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the bottom pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    put(file,"  k = "); put(file,codim); new_line(file);
    tstart(timer);
    Bottom_Create(lnkroot,codim);
    put_line("The poset created from the bottom : ");
    put_line(file,"The poset created from the bottom : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_General_Deformation_Poset(file,m,p,codim,index_poset);
    tstop(timer);
    new_line(file);
    print_times(file,timer,
                "Total time for General Pieri Homotopy Algorithm");
  end Bottom_General_Poset;

  procedure Mixed_General_Poset
              ( m,p : in natural32; k : out Link_to_Bracket;
                level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by incrementing the top and decrementing the bottom pivots.

  -- ON ENTRY :
  --   m        dimension of the input planes;
  --   p        dimension of the output planes.
  
  -- ON RETURN :
  --   k        codimension conditions, k(i) = 1 : hypersurface condition;
  --   level_poset is poset of range 0..m*p;
  --   index_poset is poset of range 0..m*p.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);

  begin
    Top_Bottom_Create(lnkroot,codim,m+p);
    put_line("The poset created in a mixed fashion : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    k := new Bracket'(codim);
  end Mixed_General_Poset;

  procedure Mixed_General_Poset
              ( file : in file_type; m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by incrementing the top and decrementing the bottom pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));

  begin
    put(file,"  k = "); put(file,codim); new_line(file);
    tstart(timer);
    Top_Bottom_Create(lnkroot,codim,m+p);
    put_line("The poset created in a mixed fashion : ");
    put_line(file,"The poset created in a mixed fashion : ");
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    Solve_General_Deformation_Poset(file,m,p,codim,index_poset);
    tstop(timer);
    new_line(file);
    print_times(file,timer,
                "Total time for General Pieri Homotopy Algorithm");
  end Mixed_General_Poset;

  function Menu_Choice return character is

  -- DESCRIPTION :
  --   Displays the menu of available choices and returns the character
  --   the user typed in.
    ans : character;

  begin
    put_line("MENU for deforming p-planes in (m+p)-space, co-dimensions k : ");
    put_line("  1. k_i = 1  consistently incrementing the top pivots.");
    put_line("  2.          consistently decrementing the bottom pivots.");
    put_line("  3.          mixed top-bottom sequence for poset creation.");
    put_line("  4. k_i >= 1 consistently incrementing the top pivots.");
    put_line("  5.          consistently decrementing the bottom pivots.");
    put_line("  6.          mixed top-bottom sequence for poset creation.");
    put("Type 1, 2, 3, 4, 5, or 6 to choose : ");
    Ask_Alternative(ans,"123456");
    return ans;
  end Menu_Choice;

  procedure Run_Pieri_Homotopies
              ( file : in file_type;
                m,p : in natural32; k : in Link_to_Bracket;
                level_poset : in Array_of_Nodes;
                index_poset : in Array_of_Array_of_Nodes ) is

    timer : Timing_Widget;

  begin
    put(file,"Pieri Homotopies for m = "); put(file,m,1);
    put(file," and p = "); put(file,p,1);
    if k = null then
      new_line(file);
      put(file,index_poset);
      tstart(timer);
      Solve_Hypersurface_Deformation_Poset(file,m,p,level_poset,index_poset);
      tstop(timer);
      new_line(file);
      print_times(file,timer,
                 "Total time for Hypersurface Pieri Homotopy Algorithm");
    else
      put(file," and k = "); Brackets_io.put(file,k.all); new_line(file);
      put(file,index_poset);
      tstart(timer);
      Solve_General_Deformation_Poset(file,m,p,k.all,index_poset);
      tstop(timer);
      new_line(file);
      print_times(file,timer,
                  "Total time for General Pieri Homotopy Algorithm");
    end if;
  end Run_Pieri_Homotopies;

  procedure Main ( n,d : in natural32 ) is

    p : constant natural32 := d;
    m : constant natural32 := n-d;
    ans : character := Menu_Choice;
    k : Link_to_Bracket;
    level_poset : Array_of_Nodes(0..integer32(m*p));
    index_poset : Array_of_Array_of_Nodes(0..integer32(m*p));
    file : file_type;

  begin
    new_line;
    case ans is
      when '1' => Top_Hypersurface_Poset(m,p,level_poset,index_poset);
      when '2' => Bottom_Hypersurface_Poset(m,p,level_poset,index_poset);
      when '3' => Mixed_Hypersurface_Poset(m,p,level_poset,index_poset);
      when '4' => Top_General_Poset(m,p,k,level_poset,index_poset);
      when '5' => Bottom_General_Poset(m,p,k,level_poset,index_poset);
      when '6' => Mixed_General_Poset(m,p,k,level_poset,index_poset);
      when others => put_line("Option not recognized.  Please try again.");
    end case;
    new_line;
    put("Do you want to run the Pieri homotopies ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file);
      new_line;
      Run_Pieri_Homotopies(file,m,p,k,level_poset,index_poset);
    end if;
  end Main;

  procedure Main ( file : in file_type; n,d : in natural32 ) is

    p : constant natural32 := d;
    m : constant natural32 := n-d;
    ans : character;

  begin
    put(file,"Pieri Homotopies for m = "); put(file,m,1);
    put(file," and p = "); put(file,p,1); new_line(file);
    ans := Menu_Choice;
    new_line;
    case ans is
      when '1' => Top_Hypersurface_Poset(file,m,p);
      when '2' => Bottom_Hypersurface_Poset(file,m,p);
      when '3' => Mixed_Hypersurface_Poset(file,m,p);
      when '4' => Top_General_Poset(file,m,p);
      when '5' => Bottom_General_Poset(file,m,p);
      when '6' => Mixed_General_Poset(file,m,p);
      when others => put_line("Option not recognized.  Please try again.");
    end case;
  end Main;

end Main_Pieri_Homotopies;
