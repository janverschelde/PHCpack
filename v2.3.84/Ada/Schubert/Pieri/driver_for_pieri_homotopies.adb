with text_io,integer_io;                 use text_io,integer_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Matrix_Inversion;          use Standard_Matrix_Inversion;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Symbol_Table;                       use Symbol_Table;
with Matrix_Indeterminates;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Homotopy;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Increment_and_Fix_Continuation;     use Increment_and_Fix_Continuation;
with Drivers_for_Poly_Continuation;      use Drivers_for_Poly_Continuation;
with Brackets,Brackets_io;               use Brackets,Brackets_io;
with Bracket_Monomials,Bracket_Systems;  use Bracket_Monomials,Bracket_Systems;
with Specialization_of_Planes;           use Specialization_of_Planes;
with Symbolic_Minor_Equations;           use Symbolic_Minor_Equations;
with Numeric_Minor_Equations;            use Numeric_Minor_Equations;
with Determinantal_Systems;              use Determinantal_Systems;
with Plane_Representations;              use Plane_Representations;
with Localization_Posets;                use Localization_Posets;
with Localization_Posets_io;             use Localization_Posets_io;
with Deformation_Posets;                 use Deformation_Posets;
with Drivers_for_Input_Planes;           use Drivers_for_Input_Planes;

procedure Driver_for_Pieri_Homotopies
                ( file : in file_type; n,d : in natural ) is

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

    oc : natural;

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
             ( m,p : natural ) return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the m-plane spanned by the first m standard basis vectors.

    res : Standard_Complex_Matrices.Matrix(1..m+p,1..m);

  begin
    for j in 1..m loop                 -- assign j-th column
      for i in 1..m+p loop
        res(i,j) := Create(0.0);       -- initialize with zeros
      end loop;
      res(j,j) := Create(1.0);         -- j-th column = j-th basis vector
    end loop;
    return res;
  end First_Standard_Plane;

  function Last_Standard_Plane
             ( m,p : natural ) return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the m-plane spanned by the last m standard basis vectors.

    res : Standard_Complex_Matrices.Matrix(1..m+p,1..m);

  begin
    for j in 1..m loop                  -- assign j-th column
      for i in 1..m+p loop
        res(i,j) := Create(0.0);        -- initialize with zeros
      end loop;
      res(p+j,j) := Create(1.0);        -- j-th vector = (p+j)-th basis vector
    end loop;
    return res;
  end Last_Standard_Plane;

  procedure Basis_Change ( n : in natural; vm : in VecMat; transvm : out VecMat;
                           trans : out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Changes basis so that the last two planes in vm are spanned by
  --   the first and last standard basis vectors respectively.

    wrk : Standard_Complex_Matrices.Matrix(1..n,1..n);
   -- penuplane : constant Standard_Complex_Matrices.Matrix
   --   := vm(vm'last-1).all;
    frstplane : constant Standard_Complex_Matrices.Matrix := vm(vm'first).all;
    lastplane : constant Standard_Complex_Matrices.Matrix := vm(vm'last).all;
    colind : natural := 0;
    use Standard_Complex_Matrices;
    ranflt : double_float;

  begin
    for j in frstplane'range(2) loop              -- fill up the work matrix
      colind := colind + 1;
      for i in frstplane'range(1) loop
         wrk(i,colind) := frstplane(i,j);
      end loop;
    end loop;
    for j in colind+1..(n-lastplane'length(2)) loop  -- random spacing
      for i in 1..n loop
        ranflt := Random;
        wrk(i,j) := Create(ranflt);
      end loop;
    end loop;
    colind := n+1;
    for j in reverse lastplane'range(2) loop
      colind := colind - 1;
      for i in lastplane'range(1) loop
         wrk(i,colind) := lastplane(i,j);
      end loop;
    end loop;
    trans := Inverse(wrk);                        -- transformation = inverse
    for i in vm'first+1..vm'last-1 loop
      transvm(i-1) := new Standard_Complex_Matrices.Matrix'
                           (trans*vm(i).all);
    end loop;
    transvm(transvm'last-1)
      := new Standard_Complex_Matrices.Matrix'(trans*vm(vm'first).all);
    transvm(transvm'last)
      := new Standard_Complex_Matrices.Matrix'(trans*vm(vm'last).all);
  end Basis_Change;

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
                ( file : in file_type; n,p : in natural;
                  sols : in Solution_List;
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
             ( n : natural; locmap : Standard_Natural_Matrices.Matrix;
               xpm : Standard_Complex_Poly_Matrices.Matrix; planes : VecMat )
             return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the polynomial system that collects the hypersurface
  --   intersection conditions for meeting the given m-planes.
  --   The system is localized according to the given localization map.

    wrksys : Poly_Sys(1..n) := Polynomial_Equations(planes(1..n),xpm);
    res : Poly_Sys(1..n) := Localize(locmap,wrksys);

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

  function Square ( n : natural; p : Poly_Sys ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns a n-by-n system, by adding random linear combinations of 
  --   the polynomials p(i), i > n, to the first n polynomials.

    res : Poly_Sys(1..n);
    acc : Poly;

  begin
    for i in res'range loop
      Copy(p(i),res(i));
      for j in n+1..p'last loop
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
    max : constant natural := 3;
    numit : natural := 0;

  begin
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,max,false,false);
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
    n : constant natural := xpm'length(1);
    p : constant natural := xpm'length(2);
    m : constant natural := n-p;
    nva : constant natural := n*p + 1;

  begin
    for i in start'range loop
      declare
        moving : Standard_Complex_Poly_Matrices.Matrix
                   (start(i)'range(1),start(i)'range(2))
               := Moving_U_Matrix(nva,start(i).all,target(i).all);
        kd : constant natural := p + start(i)'length(2);
        bm : Bracket_Monomial := Maximal_Minors(n,kd);
        bs : Bracket_System(0..Number_of_Brackets(bm))
           := Minor_Equations(kd,kd-p,bm); 
        sys : Poly_Sys(1..bs'last)
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
     then Rep_Cont(file,sols,false,Create(1.0));
     else Sil_Cont(sols,false,Create(1.0));
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
    neq : constant natural := homsys'last;       -- number of equations 
    dim : constant natural := Head_Of(sols).n;   -- actual dimension
    squhom : Poly_Sys(1..dim) := Square(dim,homsys);

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
    n : constant natural := target'last;
    a : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    b : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    dimsol : constant natural := Head_Of(sols).n;
    squ_target,squ_start : Poly_Sys(1..dimsol);

  begin
    tstart(timer);
    if dimsol = n
     then Standard_Homotopy.Create(target,start,2,a,b,true);  -- linear cheater
     else squ_target := Square(dimsol,target);
          squ_start  := Square(dimsol,start);
          Standard_Homotopy.Create(squ_target,squ_start,2,a,b,true);
    end if;
    Continuation(file,sols,report);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Linear Cheater's homotopy to target system");
    if dimsol = n
     then Refine_Roots(file,target,sols);
     else Refine_Roots(file,squ_target,sols);
          Evaluate_Roots(file,target,sols);
    end if;
  end Path_Following_with_Cheater;

  procedure Solve_Hypersurface_Target_System
                     ( file : in file_type; m,p : in natural;
                       start_planes,target_planes : in VecMat;
                       index_poset : in Array_of_Array_of_Nodes;
                       deform_poset : in Array_of_Array_of_VecMats;
                       target_level : in natural; report : in boolean ) is

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

    dim : constant natural := m*p;
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..m+p,1..p)
        := Localization_Pattern(m+p,index_poset(dim)(1).top,
                                    index_poset(dim)(1).bottom);
    solplanes : constant VecMat := deform_poset(target_level)(1).all;
    locmap : Standard_Natural_Matrices.Matrix(1..m+p,1..p)
           := Standard_Coordinate_Frame(xpm,solplanes(1).all);
    locsys : Poly_Sys(start_planes'range)
           := Create_Hypersurface_System(dim,locmap,xpm,start_planes);
    target : Poly_Sys(target_planes'range)
           := Create_Hypersurface_System(dim,locmap,xpm,target_planes);
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
    Write_Solution_Planes(file,m+p,p,sols,locmap);
  end Solve_Hypersurface_Target_System;

  procedure Solve_General_Target_System
                     ( file : in file_type; m,p : in natural; k : in Bracket;
                       start_planes,target_planes : in VecMat;
                       index_poset : in Array_of_Array_of_Nodes;
                       deform_poset : in Array_of_Array_of_VecMats;
                       target_level : in natural; report : in boolean ) is

  -- DESCRIPTION :
  --   This procedure tests the output of the deformation poset,
  --   creating polynomial representations of the intersection conditions
  --   and solution lists representing the solution planes.

  -- ON ENTRY :
  --   file            to write intermediate output on;
  --   m               dimension of the input planes;
  --   p               dimension of the solution planes;
  --   k               co-dimension conditions;
  --   start_planes    input (m+1-k(i))-planes in general position;
  --   target_planes   specific input (m+1-k(i))-planes;
  --   index_poset     indexed localization poset; 
  --   deform_poset    poset with the solution p-planes;
  --   target_level    indicates lowest node in deform_poset that is filled;
  --   report          switch for intermediate output during continuation.

    dim : constant natural := m*p;
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..m+p,1..p)
        := Localization_Pattern(m+p,index_poset(dim)(1).top,
                                    index_poset(dim)(1).bottom);
    solplanes : constant VecMat := deform_poset(target_level)(1).all;
    locmap : Standard_Natural_Matrices.Matrix(1..m+p,1..p)
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
    Write_Solution_Planes(file,m+p,p,sols,locmap);
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
      if npaths(i) /= 0
       then put(file,"|"); put(file,i,4); put(file,"   |");
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
               ( file : in file_type; m,p : in natural;
                 level_poset : in Array_of_Nodes;
                 index_poset : in Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates a deformation poset and applies the Solve operator.

    n : constant natural := m+p;
    deform_poset : Array_of_Array_of_VecMats(index_poset'range)
                 := Create(index_poset);
    planes : VecMat(1..m*p) := Random_Complex_Planes(m,p);
    target_planes : VecMat(1..m*p);
    ans : character;
    report,outlog : boolean;
    timer : Timing_Widget;
    target_level : natural := m*p;
    nbp : natural := 0;
    npaths : Standard_Natural_Vectors.Vector(1..target_level)
           := (1..target_level => 0);
    timings : Duration_Array(1..target_level) := (1..target_level => 0.0);

  begin
    put_line("The size of the deformation poset : ");
    put_line(file,"The size of the deformation poset : ");
    put_roco(index_poset);
    put_roco(file,index_poset);
    new_line;
    put("Give target level <= "); put(target_level,1);
    put(" = root level by default : "); get(target_level);
    for i in 1..target_level loop
      nbp := nbp + Row_Root_Count_Sum(level_poset,i);
    end loop;
    put("The number of paths : "); put(nbp,1); new_line;
    put(file,"The number of paths : "); put(file,nbp,1); new_line(file);
    new_line;
    put("Do you want to have the homotopies on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    outlog := (ans = 'y');
    if target_level < m*p
     then planes(m*p).all := Last_Standard_Plane(m,p);
          if target_level < m*p-1
           then planes(m*p-1).all := First_Standard_Plane(m,p);
          end if;
    end if;
    Driver_for_Input_Planes(file,m,p,target_planes); 
    Matrix_Indeterminates.Initialize_Symbols(m+p,p);
    Add_t_Symbol;
    Set_Parameters(file,report);
    tstart(timer);
    for i in index_poset(target_level)'range loop
      declare
        root : Node := index_poset(target_level)(i).all;
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
    if deform_poset(target_level)(1) /= null
     then new_line(file);
          Solve_Hypersurface_Target_System
            (file,m,p,planes,target_planes,index_poset,deform_poset,
             target_level,report);
    end if;
  end Solve_Hypersurface_Deformation_Poset;

  procedure Solve_General_Deformation_Poset
               ( file : in file_type; m,p : in natural; k : in Bracket;
                 index_poset : in Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Applies the solver to general intersection conditions.

    deform_poset : Array_of_Array_of_VecMats(index_poset'range)
                 := Create(index_poset);
    planes : VecMat(k'range) := Random_Complex_Planes(m,p,k);
    target_planes : VecMat(k'range);
    ans : character;
    report,outlog : boolean;
    timer : Timing_Widget;
    target_level : natural := m*p;
    npaths : Standard_Natural_Vectors.Vector(1..target_level)
           := (1..target_level => 0);
    timings : Duration_Array(1..target_level) := (1..target_level => 0.0);

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
    Driver_for_Input_Planes(file,m,p,k,target_planes);
    Matrix_Indeterminates.Initialize_Symbols(m+p,p);
    Add_t_Symbol;
    Set_Parameters(file,report);
    tstart(timer);
    for i in index_poset(target_level)'range loop
      declare
        root : Node := index_poset(target_level)(i).all;
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
    if deform_poset(target_level)(1) /= null
     then new_line(file);
          Solve_General_Target_System
            (file,m,p,k,planes,target_planes,index_poset,deform_poset,
             target_level,report);
    end if;
  end Solve_General_Deformation_Poset;

-- HYPERSURFACE SOLVERS :

  procedure Create_Top_Hypersurface_Poset
              ( file : in file_type; m,p : in natural ) is

  -- DESCRIPTION :
  --   Create the poset by incrementing only top pivots.

    timer : Timing_Widget;
    root : Node(p) := Trivial_Root(m,p);
    lnkroot : Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..m*p);
    index_poset : Array_of_Array_of_Nodes(0..m*p);

  begin
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
  end Create_Top_Hypersurface_Poset;

  procedure Create_Bottom_Hypersurface_Poset
              ( file : in file_type; m,p : in natural ) is

  -- DESCRIPTION :
  --   Create the poset by decrementing only bottom pivots.

    timer : Timing_Widget;
    root : Node(p) := Trivial_Root(m,p);
    lnkroot : Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..m*p);
    index_poset : Array_of_Array_of_Nodes(0..m*p);

  begin
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
  end Create_Bottom_Hypersurface_Poset;

  procedure Create_Mixed_Hypersurface_Poset
              ( file : in file_type; m,p : in natural ) is

  -- DESCRIPTION :
  --   Create the poset by incrementing top and decrementing bottom pivots.

    timer : Timing_Widget;
    root : Node(p) := Trivial_Root(m,p);
    lnkroot : Link_to_Node := new Node'(root);
    level_poset : Array_of_Nodes(0..m*p);
    index_poset : Array_of_Array_of_Nodes(0..m*p);

  begin
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
  end Create_Mixed_Hypersurface_Poset;

-- GENERAL SOLVERS :

  procedure Create_Top_General_Poset
              ( file : in file_type; m,p : in natural ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    timer : Timing_Widget;
    root : Node(p) := Trivial_Root(m,p);
    lnkroot : Link_to_Node := new Node'(root);
    codim : constant Bracket := Read_Codimensions(m,p,0);
    level_poset : Array_of_Nodes(0..m*p);
    index_poset : Array_of_Array_of_Nodes(0..m*p);

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
  end Create_Top_General_Poset;

  procedure Create_Bottom_General_Poset
              ( file : in file_type; m,p : in natural ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    timer : Timing_Widget;
    root : Node(p) := Trivial_Root(m,p);
    lnkroot : Link_to_Node := new Node'(root);
    codim : constant Bracket := Read_Codimensions(m,p,0);
    level_poset : Array_of_Nodes(0..m*p);
    index_poset : Array_of_Array_of_Nodes(0..m*p);

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
  end Create_Bottom_General_Poset;

  procedure Create_Mixed_General_Poset
              ( file : in file_type; m,p : in natural ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by incrementing the top and decrementing the bottom pivots.

    timer : Timing_Widget;
    root : Node(p) := Trivial_Root(m,p);
    lnkroot : Link_to_Node := new Node'(root);
    codim : constant Bracket := Read_Codimensions(m,p,0);
    level_poset : Array_of_Nodes(0..m*p);
    index_poset : Array_of_Array_of_Nodes(0..m*p);

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
  end Create_Mixed_General_Poset;

  procedure Main is

    p : constant natural := d;
    m : constant natural := n-d;
    ans : character;

  begin
    new_line;
    put_line("MENU for deforming p-planes in (m+p)-space, co-dimensions k : ");
    put_line("  1. k_i = 1  consistently incrementing the top pivots.");
    put_line("  2.          consistently decrementing the bottom pivots.");
    put_line("  3.          mixed top-bottom sequence for poset creation.");
    put_line("  4. k_i >= 1 consistently incrementing the top pivots.");
    put_line("  5.          consistently decrementing the bottom pivots.");
    put_line("  6.          mixed top-bottom sequence for poset creation.");
    put("Type 1, 2, 3, 4, 5, or 6 to choose : ");
    Ask_Alternative(ans,"123456");
    new_line;
    case ans is
      when '1' => new_line(file); Create_Top_Hypersurface_Poset(file,m,p);
      when '2' => new_line(file); Create_Bottom_Hypersurface_Poset(file,m,p);
      when '3' => new_line(file); Create_Mixed_Hypersurface_Poset(file,m,p);
      when '4' => Create_Top_General_Poset(file,m,p);
      when '5' => Create_Bottom_General_Poset(file,m,p);
      when '6' => Create_Mixed_General_Poset(file,m,p);
      when others => put_line("Option not recognized.  Please try again.");
    end case;
  end Main;

begin
  Main;
end Driver_for_Pieri_Homotopies;
