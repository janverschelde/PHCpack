with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;  
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Homotopy;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Brackets;                           use Brackets;
with Localization_Posets;                use Localization_Posets;
with Localization_Posets_io;             use Localization_Posets_io;
with Curves_into_Grassmannian;           use Curves_into_Grassmannian;
with Curves_into_Grassmannian_io;        use Curves_into_Grassmannian_io;
with Deformation_Posets;                 use Deformation_Posets;
with Determinantal_Systems;              use Determinantal_Systems;
with Make_Input_Planes;

package body Main_Quantum_Pieri is

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

  procedure Solve_Target_System
               ( file : in file_type; start,target : in Poly_Sys;
                 sols : in out Solution_List; report : in boolean ) is

  -- DESCRIPTION :
  --   Calls the standard continuation routines to solve a specific
  --   target system, starting at the solutions of the start system.

  -- REQUIRED : not Is_Null(sols).

    timer : Timing_Widget;
    n : constant integer32 := target'last;
    a : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    b : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);

  begin
    tstart(timer);
    Standard_Homotopy.Create(target,start,1,a,b,true);      -- linear cheater
    Set_Continuation_Parameter(sols,Create(0.0));
    declare
      procedure Sil_Cont is
        new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                            Standard_Homotopy.Diff,Standard_Homotopy.Diff);
      procedure Rep_Cont is
        new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                               Standard_Homotopy.Diff,Standard_Homotopy.Diff);
    begin
      if report
       then Rep_Cont(file,sols,false,target=>Create(1.0));
       else Sil_Cont(sols,false,target=>Create(1.0));
      end if;
    end;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Cheater's homotopy to target system");
    Refine_Roots(file,target,sols);
  end Solve_Target_System;

  procedure Solve_Hypersurface_Target_System
                     ( file : in file_type; m,p,q,nb : in natural32;
                       start_svals,target_svals
                         : in Standard_Complex_Vectors.Vector;
                       start_planes,target_planes : in VecMat;
                       index_poset : in Array_of_Array_of_Nodes;
                       deform_poset : in Array_of_Array_of_VecMats;
                       report : in boolean ) is

  -- DESCRIPTION :
  --   This procedure tests the output of the deformation poset,
  --   creating polynomial representations of the intersection conditions
  --   and solution lists representing the solution planes.

  -- ON ENTRY :
  --   file            to write intermediate output on;
  --   m               dimension of the input planes;
  --   p               dimension of the solution planes;
  --   q               degree of the maps;
  --   nb              number of maps computed;
  --   start_svals     interpolation points at the start;
  --   target_svals    interpolation points at the target;
  --   start_planes    input m-planes in general position;
  --   target_planes   specific input m-planes;
  --   index_poset     indexed localization poset; 
  --   deform_poset    poset with the solution p-planes;
  --   report          switch for intermediate output during continuation.

    dim : constant integer32 := integer32(m*p+q*(m+p));
    top : constant Bracket := index_poset(dim)(1).top;
    bot : constant Bracket := index_poset(dim)(1).bottom; 
    xpm : constant Standard_Complex_Poly_Matrices.Matrix
                     (1..integer32(m+p),1..integer32(p))
        := Symbolic_Create(m,p,q,top,bot);
    solplanes : constant VecMat := deform_poset(dim)(1).all;
    locmap : constant Standard_Natural_Matrices.Matrix
                        (1..integer32((m+p)*(q+1)),1..integer32(p))
           := Standard_Coordinate_Frame(m,p,q,top,bot,solplanes(1).all);
    locsys : constant Poly_Sys(start_planes'range)
           := Create_Polynomial_System
                (top,bot,locmap,xpm,start_svals,start_planes);
    target : constant Poly_Sys(target_planes'range)
           := Create_Polynomial_System
                (top,bot,locmap,xpm,target_svals,target_planes);
    sols : Solution_List
         := Solution_Planes(top,bot,locmap,solplanes(1..integer32(nb)));

  begin
    One_Set_up_Symbol_Table(m,p,q,top,bot);
    new_line(file);
    put(file,"The "); put(file,q,1); put(file,"-map of "); put(file,p,1);
    put_line(file,"-planes representation : ");
    put(file,xpm);
    put_line(file,"with as localization map :"); put(file,locmap);
    new_line(file);
    Reduce_Symbols(top,bot,locmap);
    put_line(file,"THE GENERIC SYSTEM : ");
    put_line(file,locsys);
    new_line(file);
    Refine_Roots(file,locsys,sols);
    new_line(file);
    put_line(file,"THE TARGET SYSTEM : ");
    put_line(file,target);
    new_line(file);
    Solve_Target_System(file,locsys,target,sols,report);
  end Solve_Hypersurface_Target_System;

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

  procedure Solve_Deformation_Poset
              ( file : in file_type; m,p,q,nb : in natural32;
                index_poset : in out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Prepares the input for the hypersurface quantum Pieri homotopies
  --   and solves the problem, writing diagnostics and results to file.

  -- ON ENTRY :
  --   file     for intermediate diagnostics and results;
  --   m        dimension of the input planes;
  --   p        dimension of the output planes;
  --   q        degree of the curves which produce p-planes;
  --   nb       number of solutions wanted;
  --   index_poset is the combinatorial structure to count all
  --            curves of degree q which produce p-planes and meet
  --            m-planes at specific interpolation points.

  -- ON RETURN :
  --   index_poset contains at its root the total number of solutions.
 
    deform_poset : Array_of_Array_of_VecMats(index_poset'range)
                 := Create(index_poset);
    dim : constant integer32 := integer32((m*p)+q*(m+p));
    input : constant VecMat(1..dim)
          := Make_Input_Planes.Random_Complex_Planes(m,p,q);
    svals : constant Standard_Complex_Vectors.Vector(1..dim)
          := Random_Vector(1,dim);
    target_planes : VecMat(1..dim);
    target_svals : Standard_Floating_Vectors.Vector(1..dim);
    comp_target_svals : Standard_Complex_Vectors.Vector(1..dim);
    root : constant Node := index_poset(dim)(1).all;
    ans : character;
    report,outlog : boolean;
    timer : Timing_Widget;
    npaths : Standard_Natural_Vectors.Vector(1..dim) := (1..dim => 0);
    timings : Duration_Array(1..dim) := (1..dim => 0.0);
    nocheater : boolean;

  begin
    new_line;
    put("Do you want to have the homotopies on file ? (y/n) "); 
    Ask_Yes_or_No(ans);
    outlog := (ans = 'y');
    Make_Input_Planes.Main(file,m,p,q,target_svals,target_planes,nocheater);
    Set_Parameters(file,report);
    tstart(timer);
    Solve(file,m+p,q,nb,deform_poset,root,input,svals,report,outlog,
          npaths,timings);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Solving along the deformation poset");
    Write_Poset_Times(file,timer,npaths,timings);
    if not nocheater then
      for i in target_svals'range loop
        comp_target_svals(i) := Create(target_svals(i));
      end loop;
      Solve_Hypersurface_Target_System
        (file,m,p,q,nb,svals,comp_target_svals,input,target_planes,
         index_poset,deform_poset,report);
    end if;
  end Solve_Deformation_Poset;

  procedure Hypersurface_Localization_Poset
              ( lnkroot : in Link_to_Node;
                level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates the posets and outputs them to the screen and on file.
  --   The posets are of range 0..nq, where nq = m*p + q*(m+p).

  begin
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put_line("The size of the poset : "); put_roco(index_poset);
  end Hypersurface_Localization_Poset;

  procedure Hypersurface_Localization_Poset
              ( file : in file_type;
                lnkroot : in Link_to_Node; m,p,q,nq : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the posets and outputs them to the screen and on file.
  --   Calls the solver afterwards.

    level_poset : Array_of_Nodes(0..integer32(nq));
    index_poset : Array_of_Array_of_Nodes(0..integer32(nq));
    nbp,nb : natural32 := 0;

  begin
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    put_line("The size of the poset : "); put_roco(index_poset);
    put_line(file,"The size of the poset : "); put_roco(file,index_poset);
    nbp := Root_Count_Sum(level_poset);
    put("The number of paths : "); put(nbp,1); new_line;
    put(file,"The number of paths : "); put(file,nbp,1); new_line(file);
    put("How many solutions do you want ? ");
    get(nb);
    if nb > 0 then
      put(file,"Computing "); put(file,nb,1);
      put_line(file," solution maps.");
      Solve_Deformation_Poset(file,m,p,q,nb,index_poset);
    else
      put_line(file,"Computing no solution maps.");
    end if;
  end Hypersurface_Localization_Poset;

  procedure General_Localization_Poset
                ( lnkroot : in Link_to_Node;
                  level_poset : out Array_of_Nodes;
                  index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates the posets and outputs them to the screen and on file.

  begin
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put_line("The size of the poset : "); put_roco(index_poset);
  end General_Localization_Poset;

  procedure General_Localization_Poset
                ( file : in file_type; lnkroot : in Link_to_Node;
                  nq : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the posets and outputs them to the screen and on file.

    level_poset : Array_of_Nodes(0..integer32(nq));
    index_poset : Array_of_Array_of_Nodes(0..integer32(nq));
    nbp : natural32;

  begin
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(index_poset);
    put(file,index_poset);
    put_line("The size of the poset : "); put_roco(index_poset);
    put_line(file,"The size of the poset : "); put_roco(file,index_poset);
    nbp := Root_Count_Sum(level_poset);
    put("The number of paths : "); put(nbp,1); new_line;
    put(file,"The number of paths : "); put(file,nbp,1); new_line(file);
   -- Solve_Deformation_Poset(file,m,p,q,index_poset);
  end General_Localization_Poset;

  procedure Top_Hypersurface_Poset
              ( m,p,q : in natural32; level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing only top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : Link_to_Node := new Node'(root);

  begin
    Q_Top_Create(lnkroot,root.bottom(integer32(p)),m+p);
    new_line;
    put_line("The poset created from the top : ");
    Hypersurface_Localization_Poset(lnkroot,level_poset,index_poset);
  end Top_Hypersurface_Poset;

  procedure Top_Hypersurface_Poset
              ( file : in file_type; m,p,q,nq : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing only top pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : Link_to_Node := new Node'(root);

  begin
    tstart(timer);
    Q_Top_Create(lnkroot,root.bottom(integer32(p)),m+p);
    put_line("The poset created from the top : ");
    put_line(file,"The poset created from the top : ");
    Hypersurface_Localization_Poset(file,lnkroot,m,p,q,nq);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Total time for Quantum Pieri Homotopy Algorithm");
  end Top_Hypersurface_Poset;

  procedure Bottom_Hypersurface_Poset
              ( m,p,q : in natural32; level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates the poset by decrementing only bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    Q_Bottom_Create(lnkroot,m+p);
    new_line;
    put_line("The poset created from the bottom : ");
    Hypersurface_Localization_Poset(lnkroot,level_poset,index_poset);
  end Bottom_Hypersurface_Poset;

  procedure Bottom_Hypersurface_Poset
              ( file : in file_type; m,p,q,nq : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by decrementing only bottom pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    tstart(timer);
    Q_Bottom_Create(lnkroot,m+p);
    put_line("The poset created from the bottom : ");
    put_line(file,"The poset created from the bottom : ");
    Hypersurface_Localization_Poset(file,lnkroot,m,p,q,nq);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Total time for Quantum Pieri Homotopy Algorithm");
  end Bottom_Hypersurface_Poset;

  procedure Mixed_Hypersurface_Poset
              ( m,p,q : in natural32; level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing top and decrementing bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    Q_Top_Bottom_Create(lnkroot,root.bottom(integer32(p)),m+p);
    new_line;
    put_line("The poset created in a mixed fashion : ");
    Hypersurface_Localization_Poset(lnkroot,level_poset,index_poset);
  end Mixed_Hypersurface_Poset;

  procedure Mixed_Hypersurface_Poset
              ( file : in file_type; m,p,q,nq : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing top and decrementing bottom pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    tstart(timer);
    Q_Top_Bottom_Create(lnkroot,root.bottom(integer32(p)),m+p);
    put_line("The poset created in a mixed fashion : ");
    put_line(file,"The poset created in a mixed fashion :");
    Hypersurface_Localization_Poset(file,lnkroot,m,p,q,nq);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Total time for Quantum Pieri Homotopy Algorithm");
  end Mixed_Hypersurface_Poset;

  procedure Top_General_Poset
               ( m,p,q : in natural32; k : out Link_to_Bracket;
                 level_poset : out Array_of_Nodes;
                 index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);

  begin
    Q_Top_Create(lnkroot,codim,root.bottom(integer32(p)),m+p);
    put_line("The poset created from the top : ");
    General_Localization_Poset(lnkroot,level_poset,index_poset);
    k := new Bracket'(codim);
  end Top_General_Poset;

  procedure Top_General_Poset
               ( file : in file_type; m,p,q,nq : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing top pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);

  begin
    tstart(timer);
    Q_Top_Create(lnkroot,codim,root.bottom(integer32(p)),m+p);
    put_line("The poset created from the top : ");
    put_line(file,"The poset created from the top :");
    General_Localization_Poset(file,lnkroot,nq);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Total time for Quantum Pieri Homotopy Algorithm");
  end Top_General_Poset;

  procedure Bottom_General_Poset
               ( m,p,q : in natural32; k : out Link_to_Bracket;
                 level_poset : out Array_of_Nodes;
                 index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates the poset by decrementing bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);

  begin
    Q_Bottom_Create(lnkroot,codim,m+p);
    put_line("The poset created from the bottom : ");
    General_Localization_Poset(lnkroot,level_poset,index_poset);
    k := new Bracket'(codim);
  end Bottom_General_Poset;

  procedure Bottom_General_Poset
               ( file : in file_type; m,p,q,nq : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by decrementing bottom pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);

  begin
    tstart(timer);
    Q_Bottom_Create(lnkroot,codim,m+p);
    put_line("The poset created from the bottom : ");
    put_line(file,"The poset created from the bottom :");
    General_Localization_Poset(file,lnkroot,nq);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Total time for Quantum Pieri Homotopy Algorithm");
  end Bottom_General_Poset;

  procedure Mixed_General_Poset
              ( m,p,q : in natural32; k : out Link_to_Bracket;
                level_poset : out Array_of_Nodes;
                index_poset : out Array_of_Array_of_Nodes ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing top and decrementing bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);

  begin
    Q_Top_Bottom_Create(lnkroot,codim,root.bottom(integer32(p)),m+p);
    put_line("The poset created from the bottom : ");
    General_Localization_Poset(lnkroot,level_poset,index_poset);
    k := new Bracket'(codim);
  end Mixed_General_Poset;

  procedure Mixed_General_Poset
               ( file : in file_type; m,p,q,nq : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the poset by incrementing top and decrementing bottom pivots.

    timer : Timing_Widget;
    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);

  begin
    tstart(timer);
    Q_Top_Bottom_Create(lnkroot,codim,root.bottom(integer32(p)),m+p);
    put_line("The poset created from the bottom : ");
    put_line(file,"The poset created from the bottom :");
    General_Localization_Poset(file,lnkroot,nq);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Total time for Quantum Pieri Homotopy Algorithm");
  end Mixed_General_Poset;

  function Menu_Choice return character is

  -- DESCRIPTION :
  --   Returns the character the user types in response to the menu display.

    ans : character;

  begin
    put_line("MENU for interpolating maps of fixed degree in Grassmannian.");
    put_line("  1. k_i = 1  consistently incrementing top pivots.");
    put_line("  2.          consistently decrementing bottom pivots.");
    put_line("  3.          mixed top-bottom sequence for poset creation.");
    put_line("  4. k_i >= 1 consistently incrementing top pivots.");
    put_line("  5.          consistently incrementing bottom pivots.");
    put_line("  6.          mixed top-bottom sequence for poset creation.");
    put("Type 1, 2, 3, 4, 5, or 6 to choose : ");
    Ask_Alternative(ans,"123456");
    return ans;
  end Menu_Choice;

  procedure Run_Quantum_Pieri_Homotopies
              ( file : in file_type; m,p,q : in natural32;
                level_poset : in Array_of_Nodes;
                index_poset : in out Array_of_Array_of_Nodes ) is

    timer : Timing_Widget;
    nbp,nb : natural32 := 0;

  begin
    put(file,"Quantum Pieri for m = "); put(file,m,1);
    put(file,", p = "); put(file,p,1);
    put(file," and q = "); put(file,q,1); new_line(file);
    put(file,index_poset);
    put_line(file,"The size of the poset : "); put_roco(file,index_poset);
    nbp := Root_Count_Sum(level_poset);
    put("The number of paths : "); put(nbp,1); new_line;
    put(file,"The number of paths : "); put(file,nbp,1); new_line(file);
    put("How many solutions do you want ? "); get(nb);
    if nb > 0 then
      put(file,"Computing "); put(file,nb,1);
      put_line(file," solution maps.");
      tstart(timer);
      Solve_Deformation_Poset(file,m,p,q,nb,index_poset);
      tstop(timer);
      new_line(file);
      print_times(file,timer,"Quantum Pieri Homotopy Algorithm");
    else
      put_line(file,"Computing no solution maps.");
    end if;
  end Run_Quantum_Pieri_Homotopies;

  procedure Main ( n,d,q : in natural32 ) is

    ans : character := Menu_Choice;
    p : constant natural32 := d;
    m : constant natural32 := n-d;
    nq : constant natural32 := m*p + q*(m+p);
    k : Link_to_Bracket;
    level_poset : Array_of_Nodes(0..integer32(nq));
    index_poset : Array_of_Array_of_Nodes(0..integer32(nq));
    file : file_type;

  begin
    case ans is
      when '1' => Top_Hypersurface_Poset(m,p,q,level_poset,index_poset);
      when '2' => Bottom_Hypersurface_Poset(m,p,q,level_poset,index_poset);
      when '3' => Mixed_Hypersurface_Poset(m,p,q,level_poset,index_poset);
      when '4' => Top_General_Poset(m,p,q,k,level_poset,index_poset);
      when '5' => Bottom_General_Poset(m,p,q,k,level_poset,index_poset);
      when '6' => Mixed_General_Poset(m,p,q,k,level_poset,index_poset);
      when others => put_line("Option not recognized.  Please try again.");
    end case;
    if k = null then
      new_line;
      put("Do you want to run the quantum Pieri homotopies ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(file);
        Run_Quantum_Pieri_Homotopies(file,m,p,q,level_poset,index_poset);
      end if;
    end if;
  end Main;

  procedure Main ( file : in file_type; n,d,q : in natural32 ) is

    p : constant natural32 := d;
    m : constant natural32 := n-d;
    nq : constant natural32 := m*p + q*(m+p);
    ans : character;

  begin
    put(file,"Quantum Pieri for m = "); put(file,m,1);
    put(file,", p = "); put(file,p,1);
    put(file," and q = "); put(file,q,1); new_line(file);
    new_line;
    ans := Menu_Choice;
    new_line;
    case ans is
      when '1' => Top_Hypersurface_Poset(file,m,p,q,nq);
      when '2' => Bottom_Hypersurface_Poset(file,m,p,q,nq);
      when '3' => Mixed_Hypersurface_Poset(file,m,p,q,nq);
      when '4' => Top_General_Poset(file,m,p,q,nq);
      when '5' => Bottom_General_Poset(file,m,p,q,nq);
      when '6' => Mixed_General_Poset(file,m,p,q,nq);
      when others => put_line("Option not recognized.  Please try again.");
    end case;
  end Main;

end Main_Quantum_Pieri;
