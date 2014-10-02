with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Natural_Matrices;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Symbol_Table;
with Matrix_Indeterminates;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io; use Standard_Complex_Poly_Matrices_io;
with Standard_Homotopy;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with Standard_IncFix_Continuation;
with Drivers_for_Poly_Continuation;     use Drivers_for_Poly_Continuation;
with Brackets;
with Checker_Boards_io;
with Checker_Moves;
with Checker_Posets_io;
with Checker_Localization_Patterns;
with Checker_Homotopies;
with Moving_Flag_Homotopies;            use Moving_Flag_Homotopies;

package body Moving_Flag_Continuation is

-- AUXILIARIES :

  function Create ( x : Standard_Complex_Vectors.Vector ) return Solution is

  -- DESCRIPTION :
  --   Returns the solution representation of the vector x.

    res : Solution(x'last-1);

  begin
    res.t := x(x'last);
    res.m := 1;
    res.v := x(x'first..x'last-1);
    res.err := 0.0;
    res.rco := 1.0;
    res.res := 0.0;
    return res;
  end Create;

  function Create ( x : Standard_Complex_Vectors.Vector ) 
                  return Solution_List is

  -- DESCRIPTION :
  --   Returns the solution list representation of the vector x.

    sol : constant Solution(x'last-1) := Create(x);
    res : Solution_List;

  begin
    Add(res,sol);
    return res;
  end Create;

-- TARGET ROUTINES :

  procedure Set_Parameters ( file : in file_type; report : out boolean ) is

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

  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32; h : in Poly_Sys;
                xt : in out Standard_Complex_Vectors.Vector;
                sol : out Link_to_Solution ) is

    use Standard_IncFix_Continuation;

    sols : Solution_List := Create(xt);

    procedure Track is
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                             Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(h,n+1);
    Track(file,sols,false,Create(1.0));
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
  end Call_Path_Trackers;

  procedure Track_Paths
              ( file : in file_type; n : in integer32; h : in Poly_Sys ) is

    x : Standard_Complex_Vectors.Vector(1..n);
    xt : Standard_Complex_Vectors.Vector(1..n+1);
    y : Standard_Complex_Vectors.Vector(h'range);
    fail : boolean;
    res : double_float;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : Standard_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := 1.0E-8;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    ls : Link_to_Solution;
    deflate : boolean := false;

  begin
    Start_Solution(h,fail,x,res);
    if fail then
      put_line(file,"no start solution found...");
    else
      new_line(file);
      put(file,"The residual of the start solution : ");
      put(file,res,3); new_line(file);
      xt(x'range) := x; xt(xt'last) := Create(0.0);
      yh := Eval(sh,xt);
      put_line(file,"Value of the start solution at the squared homotopy :");
      put_line(file,yh);
      sols := Create(xt);
      sh0 := Eval(sh,Create(0.0),n+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,deflate,false);
      Clear(sh0); --Clear(sols);
      Call_Path_Trackers(file,n,sh,xt,ls);
      put(file,"The residual of the end solution : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file);
      new_line(file);
    end if;
    Standard_Complex_Poly_Systems.Clear(sh);
  end Track_Paths;

  procedure Track_First_Move
              ( file : in file_type; n : in integer32; h : in Poly_Sys;
                sol : out Link_to_Solution; fail : out boolean ) is

    x : Standard_Complex_Vectors.Vector(1..n);
    xt : Standard_Complex_Vectors.Vector(1..n+1);
    y : Standard_Complex_Vectors.Vector(h'range);
    res : double_float;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : Standard_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := 1.0E-8;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    Start_Solution(h,fail,x,res);
    if fail then
      put_line(file,"no start solution found...");
    else
      new_line(file);
      put(file,"The residual of the start solution : ");
      put(file,res,3); new_line(file);
      xt(x'range) := x; xt(xt'last) := Create(0.0);
      yh := Eval(sh,xt);
      put_line(file,"Value of the start solution at the squared homotopy :");
      put_line(file,yh);
      sols := Create(xt);
      sh0 := Eval(sh,Create(0.0),n+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,deflate,false);
      Clear(sh0); --Clear(sols);
      Call_Path_Trackers(file,n,sh,xt,sol);
      put(file,"The residual of the end solution : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > epsfa);
    end if;
    Standard_Complex_Poly_Systems.Clear(sh);
  exception
    when others => put_line("exception in Track_First_Move"); raise;
  end Track_First_Move;

  procedure Track_Next_Move
              ( file : in file_type; n : in integer32; h : in Poly_Sys;
                sol : in out Link_to_Solution; fail : out boolean ) is

    xt : Standard_Complex_Vectors.Vector(1..n+1);
    y : Standard_Complex_Vectors.Vector(h'range);
    res : double_float;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : Standard_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := 1.0E-8;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    xt(sol.v'range) := sol.v; xt(xt'last) := Create(0.0);
    y := Eval(h,xt);
    new_line(file);
    put_line(file,"Value of the start solution at the original homotopy :");
    put_line(file,y);
    res := Max_Norm(y);
    put(file,"The residual : "); put(file,res,3); new_line(file);
    fail := (res > tolsing);
    yh := Eval(sh,xt);
    put_line(file,"Value of the start solution at the squared homotopy :");
    put_line(file,yh);
    res := Max_Norm(y);
    put(file,"The residual : "); put(file,res,3); new_line(file);
    fail := fail and (res > tolsing);
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      sols := Create(xt);
      sh0 := Eval(sh,Create(0.0),n+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,deflate,false);
      Call_Path_Trackers(file,n,sh,xt,sol);
      put(file,"The residual of the end solution at original homotopy : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > tolsing);
    end if;
    Clear(sh); Clear(sh0); --Clear(sols);
  end Track_Next_Move;

  procedure Initialize_Symbol_Table
              ( n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                dim : out integer32 ) is

  -- DESCRIPTION :
  --   Uses the localization pattern to initialize the symbol table.
  --   On return in dim is the number of free variables.

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);

  begin
    dim := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    if not Symbol_Table.Empty
     then Symbol_Table.Clear;
    end if;
  end Initialize_Symbol_Table;

  procedure Generalizing_Homotopy
              ( file : in file_type; n,k : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,nf : in Standard_Complex_Matrices.Matrix;
                h : out Link_to_Poly_Sys; dim : out integer32 ) is

   -- f : constant natural := Checker_Moves.Falling_Checker(p);
   -- a : constant natural := Checker_Moves.Ascending_Checker(p,f);
   -- t : Standard_Natural_Matrices.Matrix(1..n,1..n)
   --   := Checker_Localization_Patterns.Transformation(n,q(f));
   -- m : Standard_Natural_Matrices.Matrix(1..n,1..n)
   --   := Checker_Localization_Patterns.Moving_Flag(p);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
        -- rows and cols must be rows and cols with p, not with q!
        --   := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
        -- must use q (current) not p (previous)

  begin
   -- Checker_Boards_io.Write_Permutation(file,p,rows,cols,m,t);
    dim := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
   -- put(file,"The degree of freedom of the localization map : ");
   -- put(file,dim,1); new_line(file);
    if not Symbol_Table.Empty
     then Symbol_Table.Clear;
    end if;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
    if cond'last = 1 then
      declare
        c : constant Brackets.Bracket(1..k)
          := Brackets.Bracket(cond(cond'first).all);
      begin
        One_Flag_Homotopy(file,n,k,q,p,rows,cols,c,vf(vf'first).all,mf,nf,h);
      end;
    else
      Moving_Flag_Homotopy(file,n,k,q,p,rows,cols,cond,vf,mf,nf,h);
    end if;
    put(file,"The moving flag homotopy has ");
    put(file,h'last,1); put(file," equations in ");
    put(file,dim,1); put_line(file,"+1 unknowns ..."); -- put(file,h.all);
  end Generalizing_Homotopy;

  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                x : in Standard_Complex_Vectors.Vector ) is

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : natural32;
    f : Link_to_Poly_Sys;

  begin
    dim := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    if not Symbol_Table.Empty
     then Symbol_Table.Clear;
    end if;
    Matrix_Indeterminates.Initialize_Symbols(dim,locmap);
   -- note: p and parameters mf and nf needed for this call ...
   -- Flag_Conditions(file,n,k,q,p,rows,cols,cond,vf,mf,nf,f);
    Flag_Conditions(n,k,q,rows,cols,cond,vf,f);
    put(file,"At q = "); put(file,q);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols); new_line(file);
   -- put_line(file,"THE FIXED FLAGS :");
   -- for i in vf'range loop
   --   put(file,vf(i).all,3);
   -- end loop;
   -- put_line(file,"THE POLYNOMIAL SYSTEM : "); put_line(file,f.all);
    put_line(file,"Verification of intersection conditions :");
    declare
      z : Standard_Complex_Vectors.Vector(x'range);
      fail : boolean;
      res : double_float;
      y : constant Standard_Complex_Vectors.Vector(f'range) := Eval(f.all,x);
    begin
      put_line(file,"The given solution :"); put_line(file,x);
      put_line(file,"The value of the given solution :"); put_line(file,y);
      First_Solution(f.all,fail,z,res);
      if fail then
        put_line(file,"failed to recompute the solution");
      else
        put_line(file,"The recomputed solution :"); put_line(file,z);
        put(file,"with residual :"); put(file,res,3); new_line(file);
      end if;
    end;
  exception
    when others => put_line("exception in verify_intersection conditions");
                   raise;
  end Verify_Intersection_Conditions;

  procedure Track_Game 
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf : in Standard_Complex_Matrices.Matrix ) is

    m : constant integer32
      := integer32(Checker_Moves.Number_of_Moves(natural32(n))) - 1;
    p,q : Standard_Natural_Vectors.Vector(1..n);
    r,c : Standard_Natural_Vectors.Vector(1..k);
    mp : Standard_Natural_Matrices.Matrix(1..n,1..n);
    nf : constant Standard_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    lnd : Checker_Posets.Link_to_Node;
    ans : character;
    h : Link_to_Poly_Sys;
    dim : integer32;

  begin
    new_line(file);
    put(file,"Number of moves in the game : ");
    put(file,m,1); new_line(file);
    for i in reverse ps.black'first+1..ps.black'last loop
      lnd := ps.white(i);
      loop
        r := lnd.rows; c := lnd.cols;
        put("make next move at current node (");
        Checker_Boards_io.Write_Bracket(r); put(",");
        Checker_Boards_io.Write_Bracket(c); put(") ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans = 'y');
        if ans /= 'y' then
          put("change current node to next sibling ? (y/n) ");
          Ask_Yes_or_No(ans);
          if ans = 'y' then
            lnd := ps.white(i).next_sibling;
          end if;
        end if;
      end loop;
      put(file,"Move "); put(file,m-i+2,1);
      p := ps.black(i).all; q := ps.black(i-1).all;
      put(file," from"); put(file,p);
      put(file," to"); put(file,q); put_line(file," :");
      r := lnd.rows; c := lnd.cols;
      Generalizing_Homotopy(file,n,k,q,p,r,c,cond,vf,mf,nf,h,dim);
      Track_Paths(file,dim,h.all);
    end loop;
    put_line("Final configuration :");
    p := ps.black(ps.black'first).all;
    r := ps.white(ps.white'first).rows;
    c := ps.white(ps.white'first).cols;
    mp := Checker_Localization_Patterns.Moving_Flag(p);
    Checker_Boards_io.Write_Permutation(file,p,r,c,mp);
    put_line("The final flag : "); put(mf,3);
    put_line("The accumulated flag : "); put(nf,3);
  end Track_Game;

  procedure Copy_Flags ( src : in Standard_Complex_VecMats.VecMat;
                         dst : in out Standard_Complex_VecMats.VecMat ) is
  begin
    for i in src'range loop
      declare
        A : constant Standard_Complex_Matrices.Link_to_Matrix := src(i);
      begin
        dst(i) := new Standard_Complex_Matrices.Matrix'(A.all);
      end;
    end loop;
  end Copy_Flags;

  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in out Standard_Complex_VecMats.VecMat;
                ls : in out Link_to_Solution; fail : out boolean ) is

    gh : Link_to_Poly_Sys;
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(q,qr,qc));
    x : Standard_Complex_Vectors.Vector(1..dim);
    res : double_float;

  begin
    fail := false;
    if ind = 0 then
      Flag_Conditions(n,k,q,qr,qc,cond,vf,gh);
      First_Solution(gh.all,fail,x,res);
      put_line(file,"The first solution :"); put_line(file,x);
      put(file,"Residual of first solution : "); put(file,res,3);
      if fail then
        put_line(file," failed first solution.");
      else
        put_line(file," found first solution.");
        declare
          sol : Solution(dim);
        begin
          sol.t := Create(0.0); sol.m := 1; sol.v := x;
          sol.err := 0.0; sol.res := res; sol.rco := 1.0;
          ls := new Solution'(sol);
        end;
      end if;
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
     -- Checker_Homotopies.Inverse_Coordinate_Transformation(ctr,vf);
      Checker_Homotopies.Trivial_Stay_Coordinates
        (file,n,k,ctr,q,p,qr,qc,pr,pc,ls.v);
      put_line(file,"Verifying after coordinate changes ...");
      Verify_Intersection_Conditions(file,n,k,q,qr,qc,cond,vf,ls.v);
    end if;
    Clear(gh);
  end Trivial_Stay;

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in out Standard_Complex_VecMats.VecMat;
                mf,nf : in Standard_Complex_Matrices.Matrix;
                ls : in out Link_to_Solution; fail : out boolean ) is

    gh : Link_to_Poly_Sys;
    dim : integer32;

  begin
    fail := true;
    if ind = 0 then
      Generalizing_Homotopy(file,n,k,q,p,pr,pc,cond,vf,mf,nf,gh,dim);
      Track_First_Move(file,dim,gh.all,ls,fail);
    else
      declare
        xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
           := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
        locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
               := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
        dim : constant integer32
            := integer32(Checker_Localization_Patterns.Degree_of_Freedom
                           (locmap));
      begin
        Initialize_Homotopy_Symbols(natural32(dim),locmap);
        put_line(file,"The moving coordinates : "); put(file,xp);
        Flag_Conditions(n,k,xp,cond,vf,gh);
        Track_Next_Move(file,dim,gh.all,ls,fail);
        Standard_Complex_Poly_Matrices.Clear(xp);
      end;
    end if;
    if not fail then
      put(file,"Transforming input planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Inverse_Coordinate_Transformation(ctr,vf);
      Checker_Homotopies.Homotopy_Stay_Coordinates(file,n,k,ctr,q,qr,qc,ls.v);
      put_line(file,"Verifying after coordinate changes ...");
      Verify_Intersection_Conditions(file,n,k,q,qr,qc,cond,vf,ls.v);
    end if;
    Clear(gh);
  exception
    when others => put_line("exception in Stay_Homotopy"); raise;
  end Stay_Homotopy;

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in out Standard_Complex_VecMats.VecMat;
                ls : in out Link_to_Solution; fail : out boolean ) is

   -- s : constant natural := Checker_Homotopies.Swap_Column(ctr,qr);
    big_r : constant integer32 := Checker_Homotopies.Swap_Checker(q,qr,qc);
    dc : constant integer32 := Checker_Moves.Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    s : constant integer32 := Checker_Homotopies.Swap_Column(ctr,locmap);
    xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Swap_Moving_Plane(file,n,k,ctr,big_r,s,q,p,pr,pc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    gh : Link_to_Poly_Sys;

  begin
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
    put_line(file,"The moving coordinates : "); put(file,xp); 
    fail := true;
    Flag_Conditions(n,k,xp,cond,vf,gh);
    if ind = 0
     then Track_First_Move(file,dim,gh.all,ls,fail);
     else Track_Next_Move(file,dim,gh.all,ls,fail);
    end if;
    if not fail then
      put(file,"Transforming input planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Inverse_Coordinate_Transformation(ctr,vf);
      if big_r > ctr + 1
       then Checker_Homotopies.First_Swap_Coordinates
             -- (file,n,k,ctr,big_r,s,q,qr,qc,ls.v);
              (file,n,k,ctr,big_r,dc,s,q,p,qr,qc,pr,pc,ls.v);
       else Checker_Homotopies.Second_Swap_Coordinates
              (file,n,k,ctr,s,q,qr,qc,ls.v);
      end if;
      put_line(file,"Verifying after coordinate changes ...");
      Verify_Intersection_Conditions(file,n,k,q,qr,qc,cond,vf,ls.v);
    end if;
    Standard_Complex_Poly_Systems.Clear(gh);
    Standard_Complex_Poly_Matrices.Clear(xp);
  exception
    when others => put_line("exception in Swap_Homotopy"); raise;
  end Swap_Homotopy;

  procedure Track_Path_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf : in Standard_Complex_Matrices.Matrix;
                tvf : out Standard_Complex_VecMats.VecMat;
                ls : in out Link_to_Solution; unhappy : out boolean ) is

    leaf : constant Link_to_Node := path(path'first);
    ip : constant Standard_Natural_Vectors.Vector(1..n) 
       := Checker_Moves.Identity_Permutation(natural32(n));
    p,q : Standard_Natural_Vectors.Vector(1..n);
    pr,pc,qr,qc : Standard_Natural_Vectors.Vector(1..k);
    cnd : constant Standard_Natural_Vectors.Vector(1..k)
        := cond(cond'first).all;
    nf : constant Standard_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    dim,ptr,homtp,ctr,ind : integer32;
    stay_child : boolean;
    fail : boolean := false;
    work_vf : Standard_Complex_VecMats.VecMat(vf'range);

  begin
    new_line(file);
    if not Checker_Moves.Happy_Checkers(ip,leaf.cols,cnd) then
      put(file,"No tracking for path "); put(file,count,1);
      put(file," because "); Checker_Posets_io.Write(file,leaf.cols,cnd);
      put_line(file," is not happy.");
      unhappy := true;
    else
      Copy_Flags(vf,work_vf);
      unhappy := false;
      put(file,"Tracking path "); put(file,count,1);
      put(file," in poset starting at a happy ");
      Checker_Posets_io.Write(file,leaf.cols,cnd); put_line(file," ...");
      p := ps.black(ps.black'last).all;
      pr := leaf.rows; pc := leaf.cols;
      for i in path'first+1..path'last loop
        ptr := ps.black'last - i + 1;
        p := ps.black(ptr+1).all; pr := path(i-1).rows; pc := path(i-1).cols;
        q := ps.black(ptr).all; qr := path(i).rows; qc := path(i).cols;
        Checker_Posets_io.Write_Node_in_Path(file,n,k,ps,path,i);
        stay_child := Checker_Posets.Is_Stay_Child(path(i).all,path(i-1).all);
        Checker_Homotopies.Define_Generalizing_Homotopy
           (file,n,q,qr,qc,stay_child,homtp,ctr);
       -- Checker_Homotopies.Define_Generalizing_Homotopy
       --    (file,n,p,pr,pc,stay_child,homtp,ctr);
       -- Checker_Homotopies.Define_Specializing_Homotopy
       --    (file,n,q,qr,qc,homtp,ctr);
        Initialize_Symbol_Table(n,k,q,qr,qc,dim);
        ind := i-path'first-1; -- ind = 0 signals start solution
        if homtp = 0 then
          Trivial_Stay(file,n,k,ctr,ind,q,p,qr,qc,pr,pc,cond,work_vf,ls,fail);
        elsif homtp = 1 then
          Stay_Homotopy(file,n,k,ctr,ind,q,p,qr,qc,pr,pc,cond,
                        work_vf,mf,nf,ls,fail);
        else -- homtp = 2
          Moving_Flag_Homotopies.Add_t_Symbol;
          Swap_Homotopy(file,n,k,ctr,ind,q,p,qr,qc,pr,pc,cond,work_vf,ls,fail);
        end if;
        if fail then
          put_line(file,"no longer a valid solution, abort tracking");
          new_line(file);
          exit;
        end if;
      end loop;
      Copy_Flags(work_vf,tvf);
    end if;
  end Track_Path_in_Poset;

  procedure Track_All_Paths_in_Poset
              ( file : in file_type; n,k : in integer32; ps : in Poset;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf : in Standard_Complex_Matrices.Matrix;
                tvf : out Standard_Complex_VecMats.VecMat;
                sols : out Solution_List ) is

    cnt : integer32 := 0;
    sols_last : Solution_List := sols;

    procedure Track_Path ( nds : in Array_of_Nodes; ct : out boolean ) is

      ls : Link_to_Solution;
      fail : boolean;

    begin
      cnt := cnt + 1;
      Track_Path_in_Poset(file,n,k,ps,nds,cnt,cond,vf,mf,tvf,ls,fail);
      if not fail
       then Append(sols,sols_last,ls.all);
      end if;
      ct := true;
    end Track_Path;
    procedure Enumerate_Paths is new Enumerate_Paths_in_Poset(Track_Path);

  begin
    Enumerate_Paths(ps);
  end Track_All_Paths_in_Poset;

end Moving_Flag_Continuation;
