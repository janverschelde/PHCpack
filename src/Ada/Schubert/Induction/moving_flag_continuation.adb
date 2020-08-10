with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with Symbol_Table;
with Matrix_Indeterminates;
with Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Poly_Matrices_io;
with DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_Matrices_io;
with QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_Matrices_io;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with DoblDobl_Root_Refiners;            use DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;            use QuadDobl_Root_Refiners;
with Brackets;
with Checker_Moves;
with Checker_Localization_Patterns;
with Checker_Homotopies;
with Wrapped_Solution_Vectors;
with Wrapped_Path_Trackers;
with Wrapped_Pade_Trackers;
with Start_Flag_Homotopies;             use Start_Flag_Homotopies;
with Setup_Flag_Homotopies;             use Setup_Flag_Homotopies;
with Moving_Flag_Homotopies;            use Moving_Flag_Homotopies;
with Recondition_Swap_Homotopies;

package body Moving_Flag_Continuation is

  procedure Track_First_Move
              ( file : in file_type; n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out Standard_Complex_Solutions.Link_to_Solution;
                fail : out boolean; rpt : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Solutions;

    x : Standard_Complex_Vectors.Vector(1..n);
    xt : Standard_Complex_Vectors.Vector(1..n+1);
    y : Standard_Complex_Vectors.Vector(h'range);
    res : double_float;
    sh : Poly_Sys(1..n);
    yh : Standard_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_First_Move 1 ...");
    end if;
    Start_Solution(h,fail,x,res);
    if fail then
      put_line(file,"no start solution found...");
    else
      new_line(file);
      put(file,"Residual of the start solution : ");
      put(file,res,3); new_line(file);
      xt(x'range) := x;
      xt(xt'last) := Standard_Complex_Numbers.Create(0.0);
      if tosqr then
        sh := Square(n,h);
        yh := Eval(sh,xt);
        put_line(file,"Value of the start solution at the squared homotopy :");
        put_line(file,yh);
      end if;
      y := Eval(h,xt);
      put_line(file,"Value of the start solution at the original homotopy :");
      put_line(file,y);
      sols := Wrapped_Solution_Vectors.Create(xt);
      if sol = null then
        put_line(file,"In Track_First_Move, sol is null.");
      else
        put_line(file,"In Track_First_Move, the solution on input :");
        put(file,sol.all); new_line(file);
      end if;
      put_line(file,"The start solution in Track_First_Move : ");
      put(file,Head_Of(sols).all); new_line(file);
      sh0 := Eval(sh,Standard_Complex_Numbers.Create(0.0),n+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,deflate,false);
      Clear(sh0); --Clear(sols);
      if rpt
       then Wrapped_Pade_Trackers.Run(file,n,sh,xt,sol,false,vrblvl-1);
       else Wrapped_Path_Trackers.Run(file,n,sh,xt,sol,vrblvl-1);
      end if;
      put(file,"Residual of the end solution : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > epsfa);
    end if;
    if tosqr
     then Standard_Complex_Poly_Systems.Clear(sh);
    end if;
  exception
    when others => put_line("exception in Track_First_Move 1"); raise;
  end Track_First_Move;

  procedure Track_First_Move
              ( file : in file_type; n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Solutions;

    x : DoblDobl_Complex_Vectors.Vector(1..n);
    xt : DoblDobl_Complex_Vectors.Vector(1..n+1);
    y : DoblDobl_Complex_Vectors.Vector(h'range);
    res : double_double;
    sh : Poly_Sys(1..n);
    yh : DoblDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_First_Move 2 ...");
    end if;
    Start_Solution(h,fail,x,res);
    if fail then
      put_line(file,"no start solution found...");
    else
      new_line(file);
      put(file,"Residual of the start solution : ");
      put(file,res,3); new_line(file);
      xt(x'range) := x;
      xt(xt'last) := DoblDobl_Complex_Numbers.Create(integer(0));
      if tosqr then
        sh := Square(n,h);
        yh := Eval(sh,xt);
        put_line(file,"Value of the start solution at the squared homotopy :");
        put_line(file,yh);
      end if;
      y := Eval(h,xt);
      put_line(file,"Value of the start solution at the original homotopy :");
      put_line(file,y);
      sols := Wrapped_Solution_Vectors.Create(xt);
      if sol = null then
        put_line(file,"In Track_First_Move, sol is null.");
      else
        put_line(file,"In Track_First_Move, the solution on input :");
        put(file,sol.all); new_line(file);
      end if;
      put_line(file,"The start solution in Track_First_Move : ");
      put(file,Head_Of(sols).all); new_line(file);
      sh0 := Eval(sh,DoblDobl_Complex_Numbers.Create(integer(0)),n+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,false);
      Clear(sh0); --Clear(sols);
      Wrapped_Path_Trackers.Run(file,n,sh,xt,sol,vrblvl-1);
      put(file,"Residual of the end solution : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > epsfa);
    end if;
    if tosqr
     then DoblDobl_Complex_Poly_Systems.Clear(sh);
    end if;
  exception
    when others => put_line("exception in Track_First_Move 2"); raise;
  end Track_First_Move;

  procedure Track_First_Move
              ( file : in file_type; n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Solutions;

    x : QuadDobl_Complex_Vectors.Vector(1..n);
    xt : QuadDobl_Complex_Vectors.Vector(1..n+1);
    y : QuadDobl_Complex_Vectors.Vector(h'range);
    res : quad_double;
    sh : Poly_Sys(1..n);
    yh : QuadDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_First_Move 3 ...");
    end if;
    Start_Solution(h,fail,x,res);
    if fail then
      put_line(file,"no start solution found...");
    else
      new_line(file);
      put(file,"Residual of the start solution : ");
      put(file,res,3); new_line(file);
      xt(x'range) := x;
      xt(xt'last) := QuadDobl_Complex_Numbers.Create(integer(0));
      if tosqr then
        sh := Square(n,h);
        yh := Eval(sh,xt);
        put_line(file,"Value of the start solution at the squared homotopy :");
        put_line(file,yh);
      end if;
      y := Eval(h,xt);
      put_line(file,"Value of the start solution at the original homotopy :");
      put_line(file,y);
      sols := Wrapped_Solution_Vectors.Create(xt);
      if sol = null then
        put_line(file,"In Track_First_Move, sol is null.");
      else
        put_line(file,"In Track_First_Move, the solution on input :");
        put(file,sol.all); new_line(file);
      end if;
      put_line(file,"The start solution in Track_First_Move : ");
      put(file,Head_Of(sols).all); new_line(file);
      sh0 := Eval(sh,QuadDobl_Complex_Numbers.Create(integer(0)),n+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,false);
      Clear(sh0); --Clear(sols);
      Wrapped_Path_Trackers.Run(file,n,sh,xt,sol,vrblvl-1);
      put(file,"Residual of the end solution : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > epsfa);
    end if;
    if tosqr
     then QuadDobl_Complex_Poly_Systems.Clear(sh);
    end if;
  exception
    when others => put_line("exception in Track_First_Move 3"); raise;
  end Track_First_Move;

  procedure Track_Next_Move
              ( file : in file_type; nv : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out Standard_Complex_Solutions.Link_to_Solution;
                fail : out boolean; rpt : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Solutions;

    xt : Standard_Complex_Vectors.Vector(1..nv+1);
    y : Standard_Complex_Vectors.Vector(h'range);
    res : double_float;
    sh : Poly_Sys(1..nv);
    yh : Standard_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_Next_Move 1 ...");
    end if;
    xt(sol.v'range) := sol.v;
    xt(xt'last) := Standard_Complex_Numbers.Create(0.0);
    y := Eval(h,xt);
    new_line(file);
   -- put_line(file,"Value of the start solution at the original homotopy :");
   -- put_line(file,y);
    res := Max_Norm(y);
    put(file,"Residual of the start solution at original homotopy : ");
    put(file,res,3); new_line(file);
    fail := (res > tolsing);
    if tosqr then
      sh := Square(nv,h);
      yh := Eval(sh,xt);
     -- put_line(file,"Value of the start solution at the squared homotopy :");
     -- put_line(file,yh);
      res := Max_Norm(yh);
      put(file,"Residual of the start solution at the squared homotopy : ");
      put(file,res,3); new_line(file);
    end if;
    fail := fail and (res > tolsing);
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      sols := Wrapped_Solution_Vectors.Create(xt);
      sh0 := Eval(sh,Standard_Complex_Numbers.Create(0.0),nv+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,deflate,false);
      if rpt 
       then Wrapped_Pade_Trackers.Run(file,nv,sh,xt,sol,false,vrblvl-1);
       else Wrapped_Path_Trackers.Run(file,nv,sh,xt,sol,vrblvl-1);
      end if;
      put(file,"Residual of the end solution at original homotopy : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > tolsing);
    end if;
    if tosqr
     then Clear(sh); Clear(sh0); --Clear(sols);
    end if;
  exception
    when others => put_line("exception in Track_Next_Move 1 ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( file : in file_type; nv : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Solutions;

    xt : DoblDobl_Complex_Vectors.Vector(1..nv+1);
    y : DoblDobl_Complex_Vectors.Vector(h'range);
    res : double_double;
    sh : Poly_Sys(1..nv);
    yh : DoblDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_Next_Move 2 ...");
    end if;
    xt(sol.v'range) := sol.v;
    xt(xt'last) := DoblDobl_Complex_Numbers.Create(integer(0));
    y := Eval(h,xt);
    new_line(file);
   -- put_line(file,"Value of the start solution at the original homotopy :");
   -- put_line(file,y);
    res := Max_Norm(y);
    put(file,"Residual of the start solution at the original homotopy : ");
    put(file,res,3); new_line(file);
    fail := (res > tolsing);
    if tosqr then
      sh := Square(nv,h);
      yh := Eval(sh,xt);
     -- put_line(file,"Value of the start solution at the squared homotopy :");
     -- put_line(file,yh);
      res := Max_Norm(yh);
      put(file,"Residual of the start solution at the squared homotopy : ");
      put(file,res,3); new_line(file);
    end if;
    fail := fail and (res > tolsing);
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      sols := Wrapped_Solution_Vectors.Create(xt);
      sh0 := Eval(sh,DoblDobl_Complex_Numbers.Create(integer(0)),nv+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,false);
      Wrapped_Path_Trackers.Run(file,nv,sh,xt,sol,vrblvl-1);
      put(file,"Residual of the end solution at original homotopy : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > tolsing);
    end if;
    if tosqr
     then Clear(sh); Clear(sh0); --Clear(sols);
    end if;
  exception
    when others => put_line("exception in Track_Next_Move 2 ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( file : in file_type; nv : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sol : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Solutions;

    xt : QuadDobl_Complex_Vectors.Vector(1..nv+1);
    y : QuadDobl_Complex_Vectors.Vector(h'range);
    res : quad_double;
    sh : Poly_Sys(1..nv);
    yh : QuadDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_Next_Move 3 ...");
    end if;
    xt(sol.v'range) := sol.v;
    xt(xt'last) := QuadDobl_Complex_Numbers.Create(integer(0));
    y := Eval(h,xt);
    new_line(file);
   -- put_line(file,"Value of the start solution at the original homotopy :");
   -- put_line(file,y);
    res := Max_Norm(y);
    put(file,"Residual of the start solution at the original homtopy : ");
    put(file,res,3); new_line(file);
    fail := (res > tolsing);
    if tosqr then
      sh := Square(nv,h);
      yh := Eval(sh,xt);
     -- put_line(file,"Value of the start solution at the squared homotopy :");
     -- put_line(file,yh);
      res := Max_Norm(yh);
      put(file,"Residual of the start solution at the squared homotopy : ");
      put(file,res,3); new_line(file);
    end if;
    fail := fail and (res > tolsing);
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      sols := Wrapped_Solution_Vectors.Create(xt);
      sh0 := Eval(sh,QuadDobl_Complex_Numbers.Create(integer(0)),nv+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,false);
      Wrapped_Path_Trackers.Run(file,nv,sh,xt,sol,vrblvl-1);
      put(file,"Residual of the end solution at original homotopy : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > tolsing);
    end if;
    if tosqr
     then Clear(sh); Clear(sh0); --Clear(sols);
    end if;
  exception
    when others => put_line("exception in Track_Next_Move 3 ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( file : in file_type; nv,nt : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; rpt : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Solutions;

    xt : Standard_Complex_Vectors.Vector(1..nv+1);
    y : Standard_Complex_Vectors.Vector(h'range);
    res : double_float;
    sh : Poly_Sys(1..nv);
    yh : Standard_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    deflate : boolean := false;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;
    zero : constant Standard_Complex_Numbers.Complex_Number 
         := Standard_Complex_Numbers.Create(0.0);
    one : constant Standard_Complex_Numbers.Complex_Number 
        := Standard_Complex_Numbers.Create(1.0);

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_Next_Move 4 ...");
    end if;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := zero;
      y := Eval(h,xt);
      put(file,"The original homotopy has "); put(file,h'last,1);
      put(file," equations in "); put(file,nv,1); put(file," variables.");
      new_line(file);
     -- put_line(file,h);
     -- put_line(file,"Value of the start solution at the original homotopy :");
     -- put_line(file,y);
      res := Max_Norm(y);
      put(file,"Residual of the start solution at the original homotopy : ");
      put(file,res,3); new_line(file);
      fail := (res > tolsing);
      if tosqr then
        sh := Square(nv,h);
        yh := Eval(sh,xt);
      -- put_line(file,"Value of the start solution at the squared homotopy :");
      -- put_line(file,yh);
        res := Max_Norm(yh);
        put(file,"Residual of the start solution at the squared homotopy : ");
        put(file,res,3); new_line(file);
      end if;
      fail := fail and (res > tolsing);
      Append(xtsols,xt_sols_last,Wrapped_Solution_Vectors.Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      if tosqr then
        sh0 := Eval(sh,zero,nv+1);
        Reporting_Root_Refiner
          (file,sh0,xtsols,epsxa,epsfa,tolsing,numit,3,deflate,false);
        put(file,"Number of solutions in xtsols : ");
        put(file,Length_Of(xtsols),1); new_line(file);
        put(file,"Number of solutions in sols   : ");
        put(file,Length_Of(sols),1); new_line(file);
      end if;
      put_line(file,"Sharpening the roots on the original system ...");
      declare
        h0 : Poly_Sys(h'range) := Eval(h,zero,nv+1);
        h0fz : constant Poly_Sys := Filter_Zero_Equations(h0);
      begin
        numit := 0;
        Reporting_Root_Sharpener
          (file,h0fz,xtsols,epsxa,epsfa,tolsing,numit,3,deflate,false);
        Standard_Complex_Poly_Systems.Clear(h0);
      end;
      if tosqr then
        if nt > 0 then
          Wrapped_Path_Trackers.Multitasked_Run(file,nv,nt,sh,xtsols,sols);
        else
          if rpt then
            Wrapped_Pade_Trackers.Run(file,nv+1,sh,xtsols,false,vrblvl-1);
           -- put_line(file,"update is just a copy ...");
            Standard_Complex_Solutions.Clear(sols); sols := xtsols;
           -- put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
           -- Wrapped_Solution_Vectors.Update(sols,xtsols);
          else
            Wrapped_Path_Trackers.Run(file,nv,sh,xtsols,sols,vrblvl-1);
          end if;
        end if;
        Clear(sh0);
        sh0 := Eval(sh,one,nv+1);
      else
        if nt > 0
         then Wrapped_Path_Trackers.Multitasked_Run(file,nv,nt,h,xtsols,sols);
         else Wrapped_Path_Trackers.Run(file,nv,h,xtsols,sols,vrblvl-1);
        end if;
      end if;
      tmp := xtsols;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        if tosqr then
          put(file,"Residual of the end solution at squared homotopy : ");
          yh := Eval(sh0,ls.v); res := Max_Norm(yh);
          put(file,res,3); new_line(file);
        end if;
     -- put_line(file,"Evaluating the end solution at the original homotopy :");
        declare
          xt : Standard_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
        begin
          xt(ls.v'range) := ls.v;
          xt(xt'last) := ls.t;
          y := Eval(h,xt);
        end;
       -- put_line(file,y);
        put(file,"Residual of the end solution at original homotopy : ");
        res := Max_Norm(y);
        put(file,res,3); new_line(file);
        fail := (res > tolsing);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    if tosqr
     then Clear(sh); Clear(sh0); --Clear(sols);
    end if;
  exception
    when others => put_line("exception in Track_Next_Move 4 ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( file : in file_type; nv,nt : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Solutions;

    xt : DoblDobl_Complex_Vectors.Vector(1..nv+1);
    y : DoblDobl_Complex_Vectors.Vector(h'range);
    res : double_double;
    sh : Poly_Sys(1..nv);
    yh : DoblDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));
    one : constant DoblDobl_Complex_Numbers.Complex_Number
        := DoblDobl_Complex_Numbers.Create(integer(1));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_Next_Move 5 ...");
    end if;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := zero;
      y := Eval(h,xt);
      new_line(file);
     -- put_line(file,"Value of the start solution at the original homotopy :");
     -- put_line(file,y);
      put(file,"The original homotopy has "); put(file,h'last,1);
      put(file," equations in "); put(file,nv,1); put(file," variables.");
      new_line(file);
      res := Max_Norm(y);
      put(file,"Residual of the start solution at the original homotopy : ");
      put(file,res,3); new_line(file);
      fail := (res > tolsing);
      if tosqr then
        sh := Square(nv,h);
        yh := Eval(sh,xt);
     -- put_line(file,"Value of the start solution at the squared homotopy :");
     -- put_line(file,yh);
        res := Max_Norm(yh);
        put(file,"Residual of the start solution at the squared homotopy : ");
        put(file,res,3); new_line(file);
      end if;
      fail := fail and (res > tolsing);
      Append(xtsols,xt_sols_last,Wrapped_Solution_Vectors.Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      if tosqr then
        sh0 := Eval(sh,zero,nv+1);
        Reporting_Root_Refiner
          (file,sh0,xtsols,epsxa,epsfa,tolsing,numit,3,false);
        put(file,"Number of solutions in xtsols : ");
        put(file,Length_Of(xtsols),1); new_line(file);
        put(file,"Number of solutions in sols   : ");
        put(file,Length_Of(sols),1); new_line(file);
      end if;
      put_line(file,"Sharpening the roots on the original system ...");
      declare
        h0 : Poly_Sys(h'range) := Eval(h,zero,nv+1);
        h0fz : constant Poly_Sys := Filter_Zero_Equations(h0);
      begin
        numit := 0;
        Reporting_Root_Refiner
          (file,h0fz,xtsols,epsxa,epsfa,tolsing,numit,3,false);
        DoblDobl_Complex_Poly_Systems.Clear(h0);
      end;
      if tosqr then
        if nt > 0
         then Wrapped_Path_Trackers.Multitasked_Run(file,nv,nt,sh,xtsols,sols);
         else Wrapped_Path_Trackers.Run(file,nv,sh,xtsols,sols,vrblvl-1);
        end if;
        Clear(sh0);
        sh0 := Eval(sh,one,nv+1);
      else
        if nt > 0
         then Wrapped_Path_Trackers.Multitasked_Run(file,nv,nt,h,xtsols,sols);
         else Wrapped_Path_Trackers.Run(file,nv,h,xtsols,sols,vrblvl-1);
        end if;
      end if;
      tmp := xtsols;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        if tosqr then
          put(file,"Residual of the end solution at squared homotopy : ");
          yh := Eval(sh0,ls.v); res := Max_Norm(yh);
          put(file,res,3); new_line(file);
        end if;
     -- put_line(file,"Evaluating the end solution at the original homotopy :");
        declare
          xt : DoblDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
        begin
          xt(ls.v'range) := ls.v;
          xt(xt'last) := ls.t;
          y := Eval(h,xt);
        end;
       -- put_line(file,y);
        put(file,"Residual of the end solution at original homotopy : ");
        res := Max_Norm(y);
        put(file,res,3); new_line(file);
        fail := (res > tolsing);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    if tosqr
     then Clear(sh); Clear(sh0); --Clear(sols);
    end if;
  exception
    when others => put_line("exception in Track_Next_Move 5 ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( file : in file_type; nv,nt : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Solutions;

    xt : QuadDobl_Complex_Vectors.Vector(1..nv+1);
    y : QuadDobl_Complex_Vectors.Vector(h'range);
    res : quad_double;
    sh : Poly_Sys(1..nv);
    yh : QuadDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));
    one : constant QuadDobl_Complex_Numbers.Complex_Number
        := QuadDobl_Complex_Numbers.Create(integer(1));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_Next_Move 6 ...");
    end if;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := zero;
      y := Eval(h,xt);
      new_line(file);
     -- put_line(file,"Value of the start solution at the original homotopy :");
     -- put_line(file,y);
      put(file,"The original homotopy has "); put(file,h'last,1);
      put(file," equations in "); put(file,nv,1); put(file," variables.");
      new_line(file);
      res := Max_Norm(y);
      put(file,"Residual of the start solution at the original homotopy : ");
      put(file,res,3); new_line(file);
      fail := (res > tolsing);
      if tosqr then
        sh := Square(nv,h);
        yh := Eval(sh,xt);
     -- put_line(file,"Value of the start solution at the squared homotopy :");
     -- put_line(file,yh);
        res := Max_Norm(yh);
        put(file,"Residual of the start solution at the squared homotopy : ");
        put(file,res,3); new_line(file);
      end if;
      fail := fail and (res > tolsing);
      Append(xtsols,xt_sols_last,Wrapped_Solution_Vectors.Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      if tosqr then
        sh0 := Eval(sh,zero,nv+1);
        Reporting_Root_Refiner
          (file,sh0,xtsols,epsxa,epsfa,tolsing,numit,3,false);
        put(file,"Number of solutions in xtsols : ");
        put(file,Length_Of(xtsols),1); new_line(file);
        put(file,"Number of solutions in sols   : ");
        put(file,Length_Of(sols),1); new_line(file);
      end if;
      put_line(file,"Sharpening the roots on the original system ...");
      declare
        h0 : Poly_Sys(h'range) := Eval(h,zero,nv+1);
        h0fz : constant Poly_Sys := Filter_Zero_Equations(h0);
      begin
        numit := 0;
        Reporting_Root_Refiner
          (file,h0fz,xtsols,epsxa,epsfa,tolsing,numit,3,false);
        QuadDobl_Complex_Poly_Systems.Clear(h0);
      end;
      if tosqr then
        if nt > 0
         then Wrapped_Path_Trackers.Multitasked_Run(file,nv,nt,sh,xtsols,sols);
         else Wrapped_Path_Trackers.Run(file,nv,sh,xtsols,sols,vrblvl-1);
        end if;
        Clear(sh0);
        sh0 := Eval(sh,one,nv+1);
      else
        if nt > 0
         then Wrapped_Path_Trackers.Multitasked_Run(file,nv,nt,h,xtsols,sols);
         else Wrapped_Path_Trackers.Run(file,nv,h,xtsols,sols,vrblvl-1);
        end if;
      end if;
      tmp := xtsols;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        if tosqr then
          put(file,"Residual of the end solution at squared homotopy : ");
          yh := Eval(sh0,ls.v); res := Max_Norm(yh);
          put(file,res,3); new_line(file);
        end if;
     -- put_line(file,"Evaluating the end solution at the original homotopy :");
        declare
          xt : QuadDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
        begin
          xt(ls.v'range) := ls.v;
          xt(xt'last) := ls.t;
          y := Eval(h,xt);
        end;
       -- put_line(file,y);
        put(file,"Residual of the end solution at original homotopy : ");
        res := Max_Norm(y);
        put(file,res,3); new_line(file);
        fail := (res > tolsing);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    if tosqr
     then Clear(sh); Clear(sh0); --Clear(sols);
    end if;
  exception
    when others => put_line("exception in Track_Next_Move 6 ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( nv,nt : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; rpt : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    xt : Standard_Complex_Vectors.Vector(1..nv+1);
    sh : Poly_Sys(1..nv);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_Next_Move 7 ...");
    end if;
    fail := false;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := Standard_Complex_Numbers.Create(0.0);
      fail := fail and (ls.res > tol);
      Append(xtsols,xt_sols_last,Wrapped_Solution_Vectors.Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if not fail then
      if tosqr then
        sh := Square(nv,h);
        if nt > 0 then
          Wrapped_Path_Trackers.Multitasked_Run(nv,nt,sh,xtsols,sols);
        else
          if rpt then
            Wrapped_Pade_Trackers.Run(nv+1,sh,xtsols,vrblvl-1);
           -- just a copy ...
           -- Wrapped_Solution_Vectors.Update(sols,xtsols);
            Standard_Complex_Solutions.Clear(sols); sols := xtsols;
          else
            Wrapped_Path_Trackers.Run(nv,sh,xtsols,sols,vrblvl-1);
          end if;
        end if;
      else
        if nt > 0
         then Wrapped_Path_Trackers.Multitasked_Run(nv,nt,h,xtsols,sols);
         else Wrapped_Path_Trackers.Run(nv,h,xtsols,sols,vrblvl-1);
        end if;
      end if;
      tmp := xtsols;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        declare
          xt : Standard_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
        begin
          xt(ls.v'range) := ls.v;
          xt(xt'last) := ls.t;
        end;
        fail := fail and (ls.res > tol);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    if tosqr
     then Clear(sh);
    end if;
  exception
    when others => put_line("exception in Track_Next_Move 7 ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( nv,nt : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    xt : DoblDobl_Complex_Vectors.Vector(1..nv+1);
    sh : Poly_Sys(1..nv);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_Next_Move 8 ...");
    end if;
    fail := false;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := DoblDobl_Complex_Numbers.Create(integer(0));
      fail := fail and (ls.res > tol);
      Append(xtsols,xt_sols_last,Wrapped_Solution_Vectors.Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if not fail then
      if tosqr then
        sh := Square(nv,h);
        if nt > 0
         then Wrapped_Path_Trackers.Multitasked_Run(nv,nt,sh,xtsols,sols);
         else Wrapped_Path_Trackers.Run(nv,sh,xtsols,sols,vrblvl-1);
        end if;
      else
        if nt > 0
         then Wrapped_Path_Trackers.Multitasked_Run(nv,nt,h,xtsols,sols);
         else Wrapped_Path_Trackers.Run(nv,h,xtsols,sols,vrblvl-1);
        end if;
      end if;
      tmp := xtsols;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        declare
          xt : DoblDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
        begin
          xt(ls.v'range) := ls.v;
          xt(xt'last) := ls.t;
        end;
        fail := fail and (ls.res > tol);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    if tosqr
     then Clear(sh);
    end if;
  exception
    when others => put_line("exception in Track_Next_Move 8 ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( nv,nt : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tosqr : in boolean; tol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    xt : QuadDobl_Complex_Vectors.Vector(1..nv+1);
    sh : Poly_Sys(1..nv);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Track_Next_Move 9 ...");
    end if;
    fail := false;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := QuadDobl_Complex_Numbers.Create(integer(0));
      fail := fail and (ls.res > tol);
      Append(xtsols,xt_sols_last,Wrapped_Solution_Vectors.Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if not fail then
      if tosqr then
        sh := Square(nv,h);
        if nt > 0
         then Wrapped_Path_Trackers.Multitasked_Run(nv,nt,sh,xtsols,sols);
         else Wrapped_Path_Trackers.Run(nv,sh,xtsols,sols,vrblvl-1);
        end if;
      else
        if nt > 0
         then Wrapped_Path_Trackers.Multitasked_Run(nv,nt,h,xtsols,sols);
         else Wrapped_Path_Trackers.Run(nv,h,xtsols,sols,vrblvl-1);
        end if;
      end if;
      tmp := xtsols;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        declare
          xt : QuadDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
        begin
          xt(ls.v'range) := ls.v;
          xt(xt'last) := ls.t;
        end;
        fail := fail and (ls.res > tol);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    if tosqr
     then Clear(sh);
    end if;
  exception
    when others => put_line("exception in Track_Next_Move 9 ..."); raise;
  end Track_Next_Move;

  procedure Initialize_Symbol_Table
              ( n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                dim : out integer32 ) is

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
                h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                dim : out integer32 ) is

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

  procedure Generalizing_Homotopy
              ( n,k : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,nf : in Standard_Complex_Matrices.Matrix;
                h : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                dim : out integer32 ) is

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
    dim := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    if cond'last = 1 then
      declare
        c : constant Brackets.Bracket(1..k)
          := Brackets.Bracket(cond(cond'first).all);
      begin
        One_Flag_Homotopy(n,k,q,p,rows,cols,c,vf(vf'first).all,mf,nf,h);
      end;
    else
      Moving_Flag_Homotopy(n,k,q,p,rows,cols,cond,vf,mf,nf,h);
    end if;
  end Generalizing_Homotopy;

  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                x : in Standard_Complex_Vectors.Vector ) is

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : natural32;
    f : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    dim := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    if not Symbol_Table.Empty
     then Symbol_Table.Clear;
    end if;
    Matrix_Indeterminates.Initialize_Symbols(dim,locmap);
   -- note: p and parameters mf and nf needed for this call ...
   -- Flag_Conditions(file,n,k,q,p,rows,cols,cond,vf,mf,nf,f);
    if minrep
     then Minimal_Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
     else Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
    end if;
    put(file,"At q = "); put(file,q);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols); new_line(file);
    put_line(file,"Verification of intersection conditions :");
   -- put_line(file,"The moving flag : ");
   -- Setup_Flag_Homotopies.Write_Standard_Moving_Flag(file,mf);
    declare
      z : Standard_Complex_Vectors.Vector(x'range);
      fail : boolean;
      res : double_float;
      y : constant Standard_Complex_Vectors.Vector(f'range)
        := Standard_Complex_Poly_SysFun.Eval(f.all,x);
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
   -- Clear(f); -- executing this Clear(f) results in a crash ...
  exception
    when others => put_line("exception in verify_intersection conditions");
                   raise;
  end Verify_Intersection_Conditions;

  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                x : in DoblDobl_Complex_Vectors.Vector ) is

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : natural32;
    f : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    dim := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    if not Symbol_Table.Empty
     then Symbol_Table.Clear;
    end if;
    Matrix_Indeterminates.Initialize_Symbols(dim,locmap);
   -- note: p and parameters mf and nf needed for this call ...
   -- Flag_Conditions(file,n,k,q,p,rows,cols,cond,vf,mf,nf,f);
    if minrep
     then Minimal_Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
     else Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
    end if;
    put(file,"At q = "); put(file,q);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols); new_line(file);
    put_line(file,"Verification of intersection conditions :");
   -- put_line(file,"The moving flag : ");
   -- Setup_Flag_Homotopies.Write_DoblDobl_Moving_Flag(file,mf);
    declare
      z : DoblDobl_Complex_Vectors.Vector(x'range);
      fail : boolean;
      res : double_double;
      y : constant DoblDobl_Complex_Vectors.Vector(f'range)
        := DoblDobl_Complex_Poly_SysFun.Eval(f.all,x);
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
   -- Clear(f); -- executing this Clear(f) results in a crash ...
  exception
    when others => put_line("exception in verify_intersection conditions");
                   raise;
  end Verify_Intersection_Conditions;

  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                x : in QuadDobl_Complex_Vectors.Vector ) is

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : natural32;
    f : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    dim := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    if not Symbol_Table.Empty
     then Symbol_Table.Clear;
    end if;
    Matrix_Indeterminates.Initialize_Symbols(dim,locmap);
   -- note: p and parameters mf and nf needed for this call ...
   -- Flag_Conditions(file,n,k,q,p,rows,cols,cond,vf,mf,nf,f);
    if minrep
     then Minimal_Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
     else Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
    end if;
    put(file,"At q = "); put(file,q);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols); new_line(file);
    put_line(file,"Verification of intersection conditions :");
   -- put_line(file,"The moving flag : ");
   -- Setup_Flag_Homotopies.Write_QuadDobl_Moving_Flag(file,mf);
    declare
      z : QuadDobl_Complex_Vectors.Vector(x'range);
      fail : boolean;
      res : quad_double;
      y : constant QuadDobl_Complex_Vectors.Vector(f'range)
        := QuadDobl_Complex_Poly_SysFun.Eval(f.all,x);
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
   -- Clear(f); -- executing this Clear(f) results in a crash ...
  exception
    when others => put_line("exception in verify_intersection conditions");
                   raise;
  end Verify_Intersection_Conditions;

  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : natural32;
    f : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    dim := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    if not Symbol_Table.Empty
     then Symbol_Table.Clear;
    end if;
    Matrix_Indeterminates.Initialize_Symbols(dim,locmap);
   -- note: p and parameters mf and nf needed for this call ...
   -- Flag_Conditions(file,n,k,q,p,rows,cols,cond,vf,mf,nf,f);
    if minrep
     then Minimal_Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
     else Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
    end if;
    put(file,"At q = "); put(file,q);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols); new_line(file);
    put_line(file,"Verification of intersection conditions :");
   -- put_line(file,"The moving flag : ");
   -- Setup_Flag_Homotopies.Write_Standard_Moving_Flag(file,mf);
    fail := true; -- assume all solutions are failures
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        res : double_float;
        y : constant Standard_Complex_Vectors.Vector(f'range)
          := Eval(f.all,ls.v);
      begin
        put_line(file,"The given solution :"); put_line(file,ls.v);
       -- put_line(file,"The value of the given solution :"); put_line(file,y);
        res := Max_Norm(y);
        put(file,"Residual of the given solution : ");
        put(file,res,3); new_line(file);
        if fail
         then fail := (res > tol); -- no fail as soon as one succeeds
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
   -- Clear(f); -- executing this Clear(f) results in a crash ...
  exception
    when others => put_line("exception in verify_intersection conditions");
                   raise;
  end Verify_Intersection_Conditions;

  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : natural32;
    f : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    dim := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    if not Symbol_Table.Empty
     then Symbol_Table.Clear;
    end if;
    Matrix_Indeterminates.Initialize_Symbols(dim,locmap);
   -- note: p and parameters mf and nf needed for this call ...
   -- Flag_Conditions(file,n,k,q,p,rows,cols,cond,vf,mf,nf,f);
    if minrep
     then Minimal_Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
     else Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
    end if;
    put(file,"At q = "); put(file,q);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols); new_line(file);
    put_line(file,"Verification of intersection conditions :");
   -- put_line(file,"The moving flag : ");
   -- Setup_Flag_Homotopies.Write_DoblDobl_Moving_Flag(file,mf);
    fail := true; -- assume all solutions are failures
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        res : double_double;
        y : constant DoblDobl_Complex_Vectors.Vector(f'range)
          := Eval(f.all,ls.v);
      begin
        put_line(file,"The given solution :"); put_line(file,ls.v);
       -- put_line(file,"The value of the given solution :"); put_line(file,y);
        res := Max_Norm(y);
        put(file,"Residual of the given solution : ");
        put(file,res,3); new_line(file);
        if fail
         then fail := (res > tol); -- no fail as soon as one succeeds
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
   -- Clear(f); -- executing this Clear(f) results in a crash ...
  exception
    when others => put_line("exception in verify_intersection conditions");
                   raise;
  end Verify_Intersection_Conditions;

  procedure Verify_Intersection_Conditions
              ( file : in file_type; n,k : in integer32;
                q,rows,cols : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : natural32;
    f : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    dim := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    if not Symbol_Table.Empty
     then Symbol_Table.Clear;
    end if;
    Matrix_Indeterminates.Initialize_Symbols(dim,locmap);
   -- note: p and parameters mf and nf needed for this call ...
   -- Flag_Conditions(file,n,k,q,p,rows,cols,cond,vf,mf,nf,f);
    if minrep
     then Minimal_Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
     else Flag_Conditions(n,k,q,rows,cols,cond,mf,vf,f);
    end if;
    put(file,"At q = "); put(file,q);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols); new_line(file);
    put_line(file,"Verification of intersection conditions :");
   -- put_line(file,"The moving flag : ");
   -- Setup_Flag_Homotopies.Write_QuadDobl_Moving_Flag(file,mf);
    fail := true; -- assume all solutions are failures
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        res : quad_double;
        y : constant QuadDobl_Complex_Vectors.Vector(f'range)
          := Eval(f.all,ls.v);
      begin
        put_line(file,"The given solution :"); put_line(file,ls.v);
       -- put_line(file,"The value of the given solution :"); put_line(file,y);
        res := Max_Norm(y);
        put(file,"Residual of the given solution : ");
        put(file,res,3); new_line(file);
        if fail
         then fail := (res > tol); -- no fail as soon as one succeeds
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
   -- Clear(f); -- executing this Clear(f) results in a crash ...
  exception
    when others => put_line("exception in verify_intersection conditions");
                   raise;
  end Verify_Intersection_Conditions;

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
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                ls : in out Standard_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(q,qr,qc));
    x : Standard_Complex_Vectors.Vector(1..dim);
    res : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Trivial_Stay 1 ...");
    end if;
    fail := false;
    if ind = 0 then
      if minrep
       then Minimal_Flag_Conditions(n,k,q,qr,qc,cond,vf,gh);
       else Flag_Conditions(n,k,q,qr,qc,cond,vf,gh);
      end if;
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
          sol.t := Standard_Complex_Numbers.Create(0.0);
          sol.m := 1; sol.v := x;
          sol.err := 0.0; sol.res := res; sol.rco := 1.0;
          ls := new Solution'(sol);
        end;
      end if;
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Trivial_Stay_Coordinates
        (file,n,k,ctr,q,p,qr,qc,pr,pc,ls.v);
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions(file,n,k,q,qr,qc,minrep,cond,mf,vf,ls.v);
      end if;
    end if;
    Standard_Complex_Poly_Systems.Clear(gh);
  end Trivial_Stay;

  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(q,qr,qc));
    x : DoblDobl_Complex_Vectors.Vector(1..dim);
    res : double_double;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Trivial_Stay 2 ...");
    end if;
    fail := false;
    if ind = 0 then
      if minrep
       then Minimal_Flag_Conditions(n,k,q,qr,qc,cond,vf,gh);
       else Flag_Conditions(n,k,q,qr,qc,cond,vf,gh);
      end if;
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
          sol.t := DoblDobl_Complex_Numbers.Create(integer(0));
          sol.m := 1; sol.v := x;
          sol.err := create(0.0);
          sol.res := res;
          sol.rco := create(1.0);
          ls := new Solution'(sol);
        end;
      end if;
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Trivial_Stay_Coordinates
        (file,n,k,ctr,q,p,qr,qc,pr,pc,ls.v);
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions(file,n,k,q,qr,qc,minrep,cond,mf,vf,ls.v);
      end if;
    end if;
    DoblDobl_Complex_Poly_Systems.Clear(gh);
  end Trivial_Stay;

  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(q,qr,qc));
    x : QuadDobl_Complex_Vectors.Vector(1..dim);
    res : quad_double;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Trivial_Stay 3 ...");
    end if;
    fail := false;
    if ind = 0 then
      if minrep
       then Minimal_Flag_Conditions(n,k,q,qr,qc,cond,vf,gh);
       else Flag_Conditions(n,k,q,qr,qc,cond,vf,gh);
      end if;
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
          sol.t := QuadDobl_Complex_Numbers.Create(integer(0));
          sol.m := 1; sol.v := x;
          sol.err := create(0.0);
          sol.res := res;
          sol.rco := create(1.0);
          ls := new Solution'(sol);
        end;
      end if;
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Trivial_Stay_Coordinates
        (file,n,k,ctr,q,p,qr,qc,pr,pc,ls.v);
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions(file,n,k,q,qr,qc,minrep,cond,mf,vf,ls.v);
      end if;
    end if;
    QuadDobl_Complex_Poly_Systems.Clear(gh);
  end Trivial_Stay;

  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Trivial_Stay 4 ...");
    end if;
    put(file,"Transforming solution planes with critical row = ");
    put(file,ctr,1); put_line(file,".");
    put_line(file,"The solution given to the Trivial_Stay_Coordinates : ");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    Checker_Homotopies.Trivial_Stay_Coordinates
      (file,n,k,ctr,q,p,qr,qc,pr,pc,sols);
    if verify then
      put_line(file,"Verifying after coordinate changes ...");
      Verify_Intersection_Conditions
        (file,n,k,q,qr,qc,minrep,cond,mf,vf,sols,tol,fail);
    end if;
  end Trivial_Stay;

  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Trivial_Stay 5 ...");
    end if;
    put(file,"Transforming solution planes with critical row = ");
    put(file,ctr,1); put_line(file,".");
    put_line(file,"The solution given to the Trivial_Stay_Coordinates : ");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    Checker_Homotopies.Trivial_Stay_Coordinates
      (file,n,k,ctr,q,p,qr,qc,pr,pc,sols);
    if verify then
      put_line(file,"Verifying after coordinate changes ...");
      Verify_Intersection_Conditions
        (file,n,k,q,qr,qc,minrep,cond,mf,vf,sols,tol,fail);
    end if;
  end Trivial_Stay;

  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Trivial_Stay 6 ...");
    end if;
    put(file,"Transforming solution planes with critical row = ");
    put(file,ctr,1); put_line(file,".");
    put_line(file,"The solution given to the Trivial_Stay_Coordinates : ");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    Checker_Homotopies.Trivial_Stay_Coordinates
      (file,n,k,ctr,q,p,qr,qc,pr,pc,sols);
    if verify then
      put_line(file,"Verifying after coordinate changes ...");
      Verify_Intersection_Conditions
        (file,n,k,q,qr,qc,minrep,cond,mf,vf,sols,tol,fail);
    end if;
  end Trivial_Stay;

  procedure Trivial_Stay
              ( n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Trivial_Stay 7 ...");
    end if;
    fail := false; -- no checks anymore ...
    Checker_Homotopies.Trivial_Stay_Coordinates
      (n,k,ctr,q,p,qr,qc,pr,pc,sols);
  end Trivial_Stay;

  procedure Trivial_Stay
              ( n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Trivial_Stay 8 ...");
    end if;
    fail := false; -- no checks anymore ...
    Checker_Homotopies.Trivial_Stay_Coordinates
      (n,k,ctr,q,p,qr,qc,pr,pc,sols);
  end Trivial_Stay;

  procedure Trivial_Stay
              ( n,k,ctr : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Trivial_Stay 9 ...");
    end if;
    fail := false; -- no checks anymore ...
    Checker_Homotopies.Trivial_Stay_Coordinates
      (n,k,ctr,q,p,qr,qc,pr,pc,sols);
  end Trivial_Stay;

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                ls : in out Standard_Complex_Solutions.Link_to_Solution; 
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Stay_Homotopy 1 ...");
    end if;
    fail := true;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp);
   -- put_line(file,"the new moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_Standard_Moving_Flag(file,mf);
   -- put_line(file,"the old moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_Standard_Moving_Flag(file,start_mf);
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if ind = 0
     then Track_First_Move(file,dim,gh.all,tosqr,tol,ls,fail,rpt,vrblvl-1);
     else Track_Next_Move(file,dim,gh.all,tosqr,tol,ls,fail,rpt,vrblvl-1);
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Homotopy_Stay_Coordinates
        (file,n,k,ctr,q,qr,qc,mf,xpm,ls.v);
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions(file,n,k,q,qr,qc,minrep,cond,mf,vf,ls.v);
      end if;
    end if;
    Standard_Complex_Poly_Matrices.Clear(xp);
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Standard_Complex_Poly_Systems.Clear(gh);
  exception
    when others => put_line("exception in Stay_Homotopy"); raise;
  end Stay_Homotopy;

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution; 
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Stay_Homotopy 2 ...");
    end if;
    fail := true;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp);
   -- put_line(file,"the new moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_DoblDobl_Moving_Flag(file,mf);
   -- put_line(file,"the old moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_DoblDobl_Moving_Flag(file,start_mf);
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if ind = 0
     then Track_First_Move(file,dim,gh.all,tosqr,tol,ls,fail,vrblvl-1);
     else Track_Next_Move(file,dim,gh.all,tosqr,tol,ls,fail,vrblvl-1);
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Homotopy_Stay_Coordinates
        (file,n,k,ctr,q,qr,qc,mf,xpm,ls.v);
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions(file,n,k,q,qr,qc,minrep,cond,mf,vf,ls.v);
      end if;
    end if;
    DoblDobl_Complex_Poly_Matrices.Clear(xp);
    DoblDobl_Complex_Poly_Matrices.Clear(xpm);
    DoblDobl_Complex_Poly_Systems.Clear(gh);
  exception
    when others => put_line("exception in Stay_Homotopy"); raise;
  end Stay_Homotopy;

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution; 
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Stay_Homotopy 3 ...");
    end if;
    fail := true;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp);
   -- put_line(file,"the new moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_QuadDobl_Moving_Flag(file,mf);
   -- put_line(file,"the old moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_QuadDobl_Moving_Flag(file,start_mf);
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if ind = 0
     then Track_First_Move(file,dim,gh.all,tosqr,tol,ls,fail,vrblvl-1);
     else Track_Next_Move(file,dim,gh.all,tosqr,tol,ls,fail,vrblvl-1);
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Homotopy_Stay_Coordinates
        (file,n,k,ctr,q,qr,qc,mf,xpm,ls.v);
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions(file,n,k,q,qr,qc,minrep,cond,mf,vf,ls.v);
      end if;
    end if;
    QuadDobl_Complex_Poly_Matrices.Clear(xp);
    QuadDobl_Complex_Poly_Matrices.Clear(xpm);
    QuadDobl_Complex_Poly_Systems.Clear(gh);
  exception
    when others => put_line("exception in Stay_Homotopy"); raise;
  end Stay_Homotopy;

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Stay_Homotopy 4 ...");
    end if;
    fail := true;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp);
   -- put_line(file,"the new moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_Standard_Moving_Flag(file,mf);
   -- put_line(file,"the old moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_Standard_Moving_Flag(file,start_mf);
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(file,dim,nt,gh.all,tosqr,tol,sols,fail,rpt,vrblvl-1);
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Homotopy_Stay_Coordinates
        (file,n,k,ctr,q,qr,qc,mf,xpm,sols);
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions
          (file,n,k,q,qr,qc,minrep,cond,mf,vf,sols,tol,fail);
      end if;
    end if;
    Standard_Complex_Poly_Matrices.Clear(xp);
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Standard_Complex_Poly_Systems.Clear(gh);
  exception
    when others => put_line("exception in Stay_Homotopy"); raise;
  end Stay_Homotopy;

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Stay_Homotopy 5 ...");
    end if;
    fail := true;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp);
   -- put_line(file,"the new moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_DoblDobl_Moving_Flag(file,mf);
   -- put_line(file,"the old moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_DoblDobl_Moving_Flag(file,start_mf);
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(file,dim,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Homotopy_Stay_Coordinates
        (file,n,k,ctr,q,qr,qc,mf,xpm,sols);
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions
          (file,n,k,q,qr,qc,minrep,cond,mf,vf,sols,tol,fail);
      end if;
    end if;
    DoblDobl_Complex_Poly_Matrices.Clear(xp);
    DoblDobl_Complex_Poly_Matrices.Clear(xpm);
    DoblDobl_Complex_Poly_Systems.Clear(gh);
  exception
    when others => put_line("exception in Stay_Homotopy"); raise;
  end Stay_Homotopy;

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Stay_Homotopy 6 ...");
    end if;
    fail := true;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp);
   -- put_line(file,"the new moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_QuadDobl_Moving_Flag(file,mf);
   -- put_line(file,"the old moving flag when making the stay homotopy :");
   -- Setup_Flag_Homotopies.Write_QuadDobl_Moving_Flag(file,start_mf);
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(file,dim,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      Checker_Homotopies.Homotopy_Stay_Coordinates
        (file,n,k,ctr,q,qr,qc,mf,xpm,sols);
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions
          (file,n,k,q,qr,qc,minrep,cond,mf,vf,sols,tol,fail);
      end if;
    end if;
    QuadDobl_Complex_Poly_Matrices.Clear(xp);
    QuadDobl_Complex_Poly_Matrices.Clear(xpm);
    QuadDobl_Complex_Poly_Systems.Clear(gh);
  exception
    when others => put_line("exception in Stay_Homotopy"); raise;
  end Stay_Homotopy;

  procedure Stay_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Stay_Homotopy 7 ...");
    end if;
    fail := true;
    xpm := Moving_Flag(start_mf,xp);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(dim,nt,gh.all,tosqr,tol,sols,fail,rpt,vrblvl-1);
    if not fail then
      Checker_Homotopies.Homotopy_Stay_Coordinates
        (n,k,ctr,q,qr,qc,mf,xpm,sols);
    end if;
    Standard_Complex_Poly_Matrices.Clear(xp);
    Standard_Complex_Poly_Matrices.Clear(xpm);
    Standard_Complex_Poly_Systems.Clear(gh);
  exception
    when others => put_line("exception in Stay_Homotopy"); raise;
  end Stay_Homotopy;

  procedure Stay_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Stay_Homotopy 8 ...");
    end if;
    fail := true;
    xpm := Moving_Flag(start_mf,xp);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(dim,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
    if not fail then
      Checker_Homotopies.Homotopy_Stay_Coordinates
        (n,k,ctr,q,qr,qc,mf,xpm,sols);
    end if;
    DoblDobl_Complex_Poly_Matrices.Clear(xp);
    DoblDobl_Complex_Poly_Matrices.Clear(xpm);
    DoblDobl_Complex_Poly_Systems.Clear(gh);
  exception
    when others => put_line("exception in Stay_Homotopy"); raise;
  end Stay_Homotopy;

  procedure Stay_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Stay_Homotopy 9 ...");
    end if;
    fail := true;
    xpm := Moving_Flag(start_mf,xp);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(dim,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
    if not fail then
      Checker_Homotopies.Homotopy_Stay_Coordinates
        (n,k,ctr,q,qr,qc,mf,xpm,sols);
    end if;
    QuadDobl_Complex_Poly_Matrices.Clear(xp);
    QuadDobl_Complex_Poly_Matrices.Clear(xpm);
    QuadDobl_Complex_Poly_Systems.Clear(gh);
  exception
    when others => put_line("exception in Stay_Homotopy"); raise;
  end Stay_Homotopy;

  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Poly_Matrices.Matrix;
                ls : in out Standard_Complex_Solutions.Link_to_Solution;
                rlq : out Standard_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 ) is

    use Recondition_Swap_Homotopies;

    rowpiv : constant integer32
           := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);
    sol : Standard_Complex_Solutions.Solution(ls.n+1);

  begin
    if vrblvl > 0 then
      put("-> in moving_flag_continuation.");
      put_line("Recondition_Swap_Homotopy 1 ...");
    end if;
    put_line(file,"reconditioning the swap homotopy ...");
    pividx := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    put(file,"the index of variable x(r+1,s+1) : ");
    put(file,pividx,1); new_line(file);
    put_line(file,"The polynomial matrix on input :");
    Standard_Complex_Poly_Matrices_io.put(file,x);
    if pividx /= 0 then
      Recondition(x,locmap,dim,s);
      Insert_Scaling_Symbol(natural32(rowpiv),natural32(s+1));
      put_line(file,"the polynomial matrix for reconditioning :");
      Standard_Complex_Poly_Matrices_io.put(file,x);
      rlq := Recondition_Equation(x,s,dim+2,pividx);
      put_line(file,"the linear recondition equation :");
      Standard_Complex_Polynomials_io.put(file,rlq); new_line(file);
      sol := Recondition_Solution(ls.all,pividx,s,locmap,x);
      put_line(file,"the reconditioned solution :");
      put_vector(file,sol);
      Standard_Complex_Solutions.Clear(ls);
      ls := new Standard_Complex_Solutions.Solution'(sol);
    end if;
  end Recondition_Swap_Homotopy;

  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                rlq : out DoblDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 ) is

    use Recondition_Swap_Homotopies;

    rowpiv : constant integer32
           := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);
    sol : DoblDobl_Complex_Solutions.Solution(ls.n+1);

  begin
    if vrblvl > 0 then
      put("-> in moving_flag_continuation.");
      put_line("Recondition_Swap_Homotopy 2 ...");
    end if;
    put_line(file,"reconditioning the swap homotopy ...");
    pividx := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    put(file,"the index of variable x(r+1,s+1) : ");
    put(file,pividx,1); new_line(file);
    put_line(file,"The polynomial matrix on input :");
    DoblDobl_Complex_Poly_Matrices_io.put(file,x);
    if pividx /= 0 then
      Recondition(x,locmap,dim,s);
      Insert_Scaling_Symbol(natural32(rowpiv),natural32(s+1));
      put_line(file,"the polynomial matrix for reconditioning :");
      DoblDobl_Complex_Poly_Matrices_io.put(file,x);
      rlq := Recondition_Equation(x,s,dim+2,pividx);
      put_line(file,"the linear recondition equation :");
      DoblDobl_Complex_Polynomials_io.put(file,rlq); new_line(file);
      sol := Recondition_Solution(ls.all,pividx,s,locmap,x);
      put_line(file,"the reconditioned solution :");
      put_vector(file,sol);
      DoblDobl_Complex_Solutions.Clear(ls);
      ls := new DoblDobl_Complex_Solutions.Solution'(sol);
    end if;
  end Recondition_Swap_Homotopy;

  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                rlq : out QuadDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 ) is

    use Recondition_Swap_Homotopies;

    rowpiv : constant integer32
           := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);
    sol : QuadDobl_Complex_Solutions.Solution(ls.n+1);

  begin
    if vrblvl > 0 then
      put("-> in moving_flag_continuation.");
      put_line("Recondition_Swap_Homotopy 3 ...");
    end if;
    put_line(file,"reconditioning the swap homotopy ...");
    pividx := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    put(file,"the index of variable x(r+1,s+1) : ");
    put(file,pividx,1); new_line(file);
    put_line(file,"The polynomial matrix on input :");
    QuadDobl_Complex_Poly_Matrices_io.put(file,x);
    if pividx /= 0 then
      Recondition(x,locmap,dim,s);
      Insert_Scaling_Symbol(natural32(rowpiv),natural32(s+1));
      put_line(file,"the polynomial matrix for reconditioning :");
      QuadDobl_Complex_Poly_Matrices_io.put(file,x);
      rlq := Recondition_Equation(x,s,dim+2,pividx);
      put_line(file,"the linear recondition equation :");
      QuadDobl_Complex_Polynomials_io.put(file,rlq); new_line(file);
      sol := Recondition_Solution(ls.all,pividx,s,locmap,x);
      put_line(file,"the reconditioned solution :");
      put_vector(file,sol);
      QuadDobl_Complex_Solutions.Clear(ls);
      ls := new QuadDobl_Complex_Solutions.Solution'(sol);
    end if;
  end Recondition_Swap_Homotopy;

  procedure Recondition_Swap_Homotopy
              ( dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List;
                rlq : out Standard_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
    use Recondition_Swap_Homotopies;

   -- rowpiv : constant integer32
   --        := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);
    rcndsols : Standard_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in moving_flag_continuation.");
      put_line("Recondition_Swap_Homotopy 4 ...");
    end if;
    pividx := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    if pividx /= 0 then
      Recondition(x,locmap,dim,s);
      rlq := Recondition_Equation(x,s,dim+2,pividx);
      rcndsols := Recondition_Solutions(sols,pividx,s,locmap,x);
      Standard_Complex_Solutions.Clear(sols);
      sols := rcndsols;
    end if;
  end Recondition_Swap_Homotopy;

  procedure Recondition_Swap_Homotopy
              ( dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                rlq : out DoblDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
    use Recondition_Swap_Homotopies;

   -- rowpiv : constant integer32
   --        := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);
    rcndsols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in moving_flag_continuation.");
      put_line("Recondition_Swap_Homotopy 5 ...");
    end if;
    pividx := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    if pividx /= 0 then
      Recondition(x,locmap,dim,s);
      rlq := Recondition_Equation(x,s,dim+2,pividx);
      rcndsols := Recondition_Solutions(sols,pividx,s,locmap,x);
      DoblDobl_Complex_Solutions.Clear(sols);
      sols := rcndsols;
    end if;
  end Recondition_Swap_Homotopy;

  procedure Recondition_Swap_Homotopy
              ( dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                rlq : out QuadDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
    use Recondition_Swap_Homotopies;

   -- rowpiv : constant integer32
   --        := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);
    rcndsols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in moving_flag_continuation.");
      put_line("Recondition_Swap_Homotopy 6 ...");
    end if;
    pividx := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    if pividx /= 0 then
      Recondition(x,locmap,dim,s);
      rlq := Recondition_Equation(x,s,dim+2,pividx);
      rcndsols := Recondition_Solutions(sols,pividx,s,locmap,x);
      QuadDobl_Complex_Solutions.Clear(sols);
      sols := rcndsols;
    end if;
  end Recondition_Swap_Homotopy;

  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List;
                rlq : out Standard_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
    use Recondition_Swap_Homotopies;

    rowpiv : constant integer32
           := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);
    rcndsols : Standard_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in moving_flag_continuation.");
      put_line("Recondition_Swap_Homotopy 7 ...");
    end if;
    put_line(file,"reconditioning the swap homotopy ...");
    pividx := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    put(file,"the index of variable x(r+1,s+1) : ");
    put(file,pividx,1); new_line(file);
    put_line(file,"The polynomial matrix on input :");
    Standard_Complex_Poly_Matrices_io.put(file,x);
    if pividx /= 0 then
      Recondition(x,locmap,dim,s);
      Insert_Scaling_Symbol(natural32(rowpiv),natural32(s+1));
      put_line(file,"the polynomial matrix for reconditioning :");
      Standard_Complex_Poly_Matrices_io.put(file,x);
      rlq := Recondition_Equation(x,s,dim+2,pividx);
      put_line(file,"the linear recondition equation :");
      Standard_Complex_Polynomials_io.put(file,rlq); new_line(file);
      rcndsols := Recondition_Solutions(sols,pividx,s,locmap,x);
      put_line(file,"the reconditioned solution list :");
      put(file,Length_Of(rcndsols),natural32(Head_Of(rcndsols).n),rcndsols);
      Standard_Complex_Solutions.Clear(sols);
      sols := rcndsols;
    end if;
  end Recondition_Swap_Homotopy;

  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                rlq : out DoblDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
    use Recondition_Swap_Homotopies;

    rowpiv : constant integer32
           := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);
    rcndsols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in moving_flag_continuation.");
      put_line("Recondition_Swap_Homotopy 8 ...");
    end if;
    put_line(file,"reconditioning the swap homotopy ...");
    pividx := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    put(file,"the index of variable x(r+1,s+1) : ");
    put(file,pividx,1); new_line(file);
    put_line(file,"The polynomial matrix on input :");
    DoblDobl_Complex_Poly_Matrices_io.put(file,x);
    if pividx /= 0 then
      Recondition(x,locmap,dim,s);
      Insert_Scaling_Symbol(natural32(rowpiv),natural32(s+1));
      put_line(file,"the polynomial matrix for reconditioning :");
      DoblDobl_Complex_Poly_Matrices_io.put(file,x);
      rlq := Recondition_Equation(x,s,dim+2,pividx);
      put_line(file,"the linear recondition equation :");
      DoblDobl_Complex_Polynomials_io.put(file,rlq); new_line(file);
      rcndsols := Recondition_Solutions(sols,pividx,s,locmap,x);
      put_line(file,"the reconditioned solution list :");
      put(file,Length_Of(rcndsols),natural32(Head_Of(rcndsols).n),rcndsols);
      DoblDobl_Complex_Solutions.Clear(sols);
      sols := rcndsols;
    end if;
  end Recondition_Swap_Homotopy;

  procedure Recondition_Swap_Homotopy
              ( file : in file_type; dim,r,s : in integer32;
                locmap : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                rlq : out QuadDobl_Complex_Polynomials.Poly;
                pividx : out integer32; vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
    use Recondition_Swap_Homotopies;

    rowpiv : constant integer32
           := Checker_Localization_Patterns.Row_of_Pivot(locmap,s+1);
    rcndsols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in moving_flag_continuation.");
      put_line("Recondition_Swap_Homotopy 9 ...");
    end if;
    put_line(file,"reconditioning the swap homotopy ...");
    pividx := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    put(file,"the index of variable x(r+1,s+1) : ");
    put(file,pividx,1); new_line(file);
    put_line(file,"The polynomial matrix on input :");
    QuadDobl_Complex_Poly_Matrices_io.put(file,x);
    if pividx /= 0 then
      Recondition(x,locmap,dim,s);
      Insert_Scaling_Symbol(natural32(rowpiv),natural32(s+1));
      put_line(file,"the polynomial matrix for reconditioning :");
      QuadDobl_Complex_Poly_Matrices_io.put(file,x);
      rlq := Recondition_Equation(x,s,dim+2,pividx);
      put_line(file,"the linear recondition equation :");
      QuadDobl_Complex_Polynomials_io.put(file,rlq); new_line(file);
      rcndsols := Recondition_Solutions(sols,pividx,s,locmap,x);
      put_line(file,"the reconditioned solution list :");
      put(file,Length_Of(rcndsols),natural32(Head_Of(rcndsols).n),rcndsols);
      QuadDobl_Complex_Solutions.Clear(sols);
      sols := rcndsols;
    end if;
  end Recondition_Swap_Homotopy;

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                ls : in out Standard_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    big_r : constant integer32 := Checker_Homotopies.Swap_Checker(q,qr,qc);
    dc : constant integer32 := Checker_Moves.Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    s : constant integer32 := Checker_Homotopies.Swap_Column(ctr,locmap);
    xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Swap_Moving_Plane(file,n,k,ctr,big_r,s,q,p,pr,pc);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlq : Standard_Complex_Polynomials.Poly;
    pivot : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Swap_Homotopy 1 ...");
    end if;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    if Checker_Homotopies.Is_Zone_A_Empty(locmap,p,ctr,s,dc)
     then pivot := 0; -- no reconditioning needed
     else Recondition_Swap_Homotopy(file,dim,ctr,s,locmap,xp,ls,rlq,pivot);
    end if;
    xpm := Moving_Flag(start_mf,xp);
    if pivot /= 0 then
      put_line(file,"The moving coordinates after multiplication by M :");
      Standard_Complex_Poly_Matrices_io.put(file,xpm);
    end if;
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if pivot = 0 then
      if ind = 0
       then Track_First_Move(file,dim,gh.all,tosqr,tol,ls,fail,rpt,vrblvl-1);
       else Track_Next_Move(file,dim,gh.all,tosqr,tol,ls,fail,rpt,vrblvl-1);
      end if;
    else -- pivot /= 0, reconditioned homotopy has one extra variable
      Setup_Flag_Homotopies.Append(gh,rlq);
      put_line(file,"The reconditioned swap homotopy :");
      put_line(file,gh.all);
      if ind = 0
       then Track_First_Move(file,dim+1,gh.all,tosqr,tol,ls,fail,rpt,vrblvl-1);
       else Track_Next_Move(file,dim+1,gh.all,tosqr,tol,ls,fail,rpt,vrblvl-1);
      end if;
      declare
        nls : constant Standard_Complex_Solutions.Solution
            := Recondition_Swap_Homotopies.Rescale_Solution
                 (ls.all,s,locmap,xp,pivot);
      begin
        Standard_Complex_Solutions.Clear(ls);
        ls := new Standard_Complex_Solutions.Solution'(nls);
        Recondition_Swap_Homotopies.Remove_One_Variable(xpm,dim+1);
      end;
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      if big_r > ctr + 1
       then Checker_Homotopies.First_Swap_Coordinates
              (file,n,k,ctr,big_r,dc,s,q,p,qr,qc,pr,pc,mf,xpm,ls.v);
       else Checker_Homotopies.Second_Swap_Coordinates
              (file,n,k,ctr,s,q,qr,qc,mf,xpm,ls.v);
      end if;
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions(file,n,k,q,qr,qc,minrep,cond,mf,vf,ls.v);
      end if;
    end if;
    Standard_Complex_Poly_Systems.Clear(gh);
    Standard_Complex_Poly_Matrices.Clear(xp);
    Standard_Complex_Poly_Matrices.Clear(xpm);
  exception
    when others => put_line("exception in Swap_Homotopy"); raise;
  end Swap_Homotopy;

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    big_r : constant integer32 := Checker_Homotopies.Swap_Checker(q,qr,qc);
    dc : constant integer32 := Checker_Moves.Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    s : constant integer32 := Checker_Homotopies.Swap_Column(ctr,locmap);
    xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Swap_Moving_Plane(file,n,k,ctr,big_r,s,q,p,pr,pc);
    xpm : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlq : DoblDobl_Complex_Polynomials.Poly;
    pivot : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Swap_Homotopy 2 ...");
    end if;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    if Checker_Homotopies.Is_Zone_A_Empty(locmap,p,ctr,s,dc)
     then pivot := 0; -- no reconditioning needed
     else Recondition_Swap_Homotopy(file,dim,ctr,s,locmap,xp,ls,rlq,pivot);
    end if;
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if pivot = 0 then
      if ind = 0
       then Track_First_Move(file,dim,gh.all,tosqr,tol,ls,fail,vrblvl-1);
       else Track_Next_Move(file,dim,gh.all,tosqr,tol,ls,fail,vrblvl-1);
      end if;
    else -- pivot /= 0, reconditioned homotopy has one extra variable
      Setup_Flag_Homotopies.Append(gh,rlq);
      put_line(file,"The reconditioned swap homotopy :");
      put_line(file,gh.all);
      if ind = 0
       then Track_First_Move(file,dim+1,gh.all,tosqr,tol,ls,fail,vrblvl-1);
       else Track_Next_Move(file,dim+1,gh.all,tosqr,tol,ls,fail,vrblvl-1);
      end if;
      declare
        nls : constant DoblDobl_Complex_Solutions.Solution
            := Recondition_Swap_Homotopies.Rescale_Solution
                 (ls.all,s,locmap,xp,pivot);
      begin
        DoblDobl_Complex_Solutions.Clear(ls);
        ls := new DoblDobl_Complex_Solutions.Solution'(nls);
        Recondition_Swap_Homotopies.Remove_One_Variable(xpm,dim+1);
      end;
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      if big_r > ctr + 1
       then Checker_Homotopies.First_Swap_Coordinates
              (file,n,k,ctr,big_r,dc,s,q,p,qr,qc,pr,pc,mf,xpm,ls.v);
       else Checker_Homotopies.Second_Swap_Coordinates
              (file,n,k,ctr,s,q,qr,qc,mf,xpm,ls.v);
      end if;
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions(file,n,k,q,qr,qc,minrep,cond,mf,vf,ls.v);
      end if;
    end if;
    DoblDobl_Complex_Poly_Systems.Clear(gh);
    DoblDobl_Complex_Poly_Matrices.Clear(xp);
    DoblDobl_Complex_Poly_Matrices.Clear(xpm);
  exception
    when others => put_line("exception in Swap_Homotopy"); raise;
  end Swap_Homotopy;

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    big_r : constant integer32 := Checker_Homotopies.Swap_Checker(q,qr,qc);
    dc : constant integer32 := Checker_Moves.Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    s : constant integer32 := Checker_Homotopies.Swap_Column(ctr,locmap);
    xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Swap_Moving_Plane(file,n,k,ctr,big_r,s,q,p,pr,pc);
    xpm : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlq : QuadDobl_Complex_Polynomials.Poly;
    pivot : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Swap_Homotopy 3 ...");
    end if;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    if Checker_Homotopies.Is_Zone_A_Empty(locmap,p,ctr,s,dc)
     then pivot := 0; -- no reconditioning needed
     else Recondition_Swap_Homotopy(file,dim,ctr,s,locmap,xp,ls,rlq,pivot);
    end if;
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if pivot = 0 then
      if ind = 0
       then Track_First_Move(file,dim,gh.all,tosqr,tol,ls,fail,vrblvl-1);
       else Track_Next_Move(file,dim,gh.all,tosqr,tol,ls,fail,vrblvl-1);
      end if;
    else -- pivot /= 0, reconditioned homotopy has one extra variable
      Setup_Flag_Homotopies.Append(gh,rlq);
      put_line(file,"The reconditioned swap homotopy :");
      put_line(file,gh.all);
      if ind = 0
       then Track_First_Move(file,dim+1,gh.all,tosqr,tol,ls,fail,vrblvl-1);
       else Track_Next_Move(file,dim+1,gh.all,tosqr,tol,ls,fail,vrblvl-1);
      end if;
      declare
        nls : constant QuadDobl_Complex_Solutions.Solution
            := Recondition_Swap_Homotopies.Rescale_Solution
                 (ls.all,s,locmap,xp,pivot);
      begin
        QuadDobl_Complex_Solutions.Clear(ls);
        ls := new QuadDobl_Complex_Solutions.Solution'(nls);
        Recondition_Swap_Homotopies.Remove_One_Variable(xpm,dim+1);
      end;
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      if big_r > ctr + 1
       then Checker_Homotopies.First_Swap_Coordinates
              (file,n,k,ctr,big_r,dc,s,q,p,qr,qc,pr,pc,mf,xpm,ls.v);
       else Checker_Homotopies.Second_Swap_Coordinates
              (file,n,k,ctr,s,q,qr,qc,mf,xpm,ls.v);
      end if;
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions(file,n,k,q,qr,qc,minrep,cond,mf,vf,ls.v);
      end if;
    end if;
    QuadDobl_Complex_Poly_Systems.Clear(gh);
    QuadDobl_Complex_Poly_Matrices.Clear(xp);
    QuadDobl_Complex_Poly_Matrices.Clear(xpm);
  exception
    when others => put_line("exception in Swap_Homotopy"); raise;
  end Swap_Homotopy;

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    big_r : constant integer32 := Checker_Homotopies.Swap_Checker(q,qr,qc);
    dc : constant integer32 := Checker_Moves.Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    s : constant integer32 := Checker_Homotopies.Swap_Column(ctr,locmap);
    xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Swap_Moving_Plane(file,n,k,ctr,big_r,s,q,p,pr,pc);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlq : Standard_Complex_Polynomials.Poly;
    pivot : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Swap_Homotopy 4 ...");
    end if;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    if Checker_Homotopies.Is_Zone_A_Empty(locmap,p,ctr,s,dc)
     then pivot := 0; -- no reconditioning needed
     else Recondition_Swap_Homotopy(file,dim,ctr,s,locmap,xp,sols,rlq,pivot);
    end if;
    xpm := Moving_Flag(start_mf,xp);
    if pivot /= 0 then
      put_line(file,"The moving coordinates after multiplication by M :");
      Standard_Complex_Poly_Matrices_io.put(file,xpm);
    end if;
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if pivot = 0 then
      Track_Next_Move(file,dim,nt,gh.all,tosqr,tol,sols,fail,rpt,vrblvl-1);
    else -- pivot /= 0, reconditioned homotopy has one extra variable
      Setup_Flag_Homotopies.Append(gh,rlq);
      put_line(file,"The reconditioned swap homotopy :");
      put_line(file,gh.all);
      Track_Next_Move(file,dim+1,nt,gh.all,tosqr,tol,sols,fail,rpt,vrblvl-1);
      declare
        rsols : constant Standard_Complex_Solutions.Solution_List
              := Recondition_Swap_Homotopies.Rescale_Solutions
                   (sols,s,locmap,xp,pivot);
        use Standard_Complex_Solutions;
      begin
        Standard_Complex_Solutions.Clear(sols);
        sols := rsols;
        put_line(file,"The rescaled solutions :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
        Recondition_Swap_Homotopies.Remove_One_Variable(xpm,dim+1);
      end;
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      if big_r > ctr + 1
       then Checker_Homotopies.First_Swap_Coordinates
              (file,n,k,ctr,big_r,dc,s,q,p,qr,qc,pr,pc,mf,xpm,sols);
       else Checker_Homotopies.Second_Swap_Coordinates
              (file,n,k,ctr,s,q,qr,qc,mf,xpm,sols);
      end if;
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions
          (file,n,k,q,qr,qc,minrep,cond,mf,vf,sols,tol,fail);
      end if;
    end if;
    Standard_Complex_Poly_Systems.Clear(gh);
    Standard_Complex_Poly_Matrices.Clear(xp);
    Standard_Complex_Poly_Matrices.Clear(xpm);
  exception
    when others => put_line("exception in Swap_Homotopy"); raise;
  end Swap_Homotopy;

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    big_r : constant integer32 := Checker_Homotopies.Swap_Checker(q,qr,qc);
    dc : constant integer32 := Checker_Moves.Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    s : constant integer32 := Checker_Homotopies.Swap_Column(ctr,locmap);
    xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Swap_Moving_Plane(file,n,k,ctr,big_r,s,q,p,pr,pc);
    xpm : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlq : DoblDobl_Complex_Polynomials.Poly;
    pivot : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Swap_Homotopy 5 ...");
    end if;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    if Checker_Homotopies.Is_Zone_A_Empty(locmap,p,ctr,s,dc)
     then pivot := 0; -- no reconditioning needed
     else Recondition_Swap_Homotopy(file,dim,ctr,s,locmap,xp,sols,rlq,pivot);
    end if;
    xpm := Moving_Flag(start_mf,xp);
    if pivot /= 0 then
      put_line(file,"The moving coordinates after multiplication by M :");
      DoblDobl_Complex_Poly_Matrices_io.put(file,xpm);
    end if;
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if pivot = 0 then
      Track_Next_Move(file,dim,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
    else -- pivot /= 0, reconditioned homotopy has one extra variable
      Setup_Flag_Homotopies.Append(gh,rlq);
      put_line(file,"The reconditioned swap homotopy :");
      put_line(file,gh.all);
      Track_Next_Move(file,dim+1,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
      declare
        rsols : constant DoblDobl_Complex_Solutions.Solution_List
              := Recondition_Swap_Homotopies.Rescale_Solutions
                   (sols,s,locmap,xp,pivot);
        use DoblDobl_Complex_Solutions;
      begin
        DoblDobl_Complex_Solutions.Clear(sols);
        sols := rsols;
        put_line(file,"The rescaled solutions :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
        Recondition_Swap_Homotopies.Remove_One_Variable(xpm,dim+1);
      end;
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      if big_r > ctr + 1
       then Checker_Homotopies.First_Swap_Coordinates
              (file,n,k,ctr,big_r,dc,s,q,p,qr,qc,pr,pc,mf,xpm,sols);
       else Checker_Homotopies.Second_Swap_Coordinates
              (file,n,k,ctr,s,q,qr,qc,mf,xpm,sols);
      end if;
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions
          (file,n,k,q,qr,qc,minrep,cond,mf,vf,sols,tol,fail);
      end if;
    end if;
    DoblDobl_Complex_Poly_Systems.Clear(gh);
    DoblDobl_Complex_Poly_Matrices.Clear(xp);
    DoblDobl_Complex_Poly_Matrices.Clear(xpm);
  exception
    when others => put_line("exception in Swap_Homotopy"); raise;
  end Swap_Homotopy;

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    big_r : constant integer32 := Checker_Homotopies.Swap_Checker(q,qr,qc);
    dc : constant integer32 := Checker_Moves.Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    s : constant integer32 := Checker_Homotopies.Swap_Column(ctr,locmap);
    xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Swap_Moving_Plane(file,n,k,ctr,big_r,s,q,p,pr,pc);
    xpm : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlq : QuadDobl_Complex_Polynomials.Poly;
    pivot : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Swap_Homotopy 6 ...");
    end if;
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    if Checker_Homotopies.Is_Zone_A_Empty(locmap,p,ctr,s,dc)
     then pivot := 0; -- no reconditioning needed
     else Recondition_Swap_Homotopy(file,dim,ctr,s,locmap,xp,sols,rlq,pivot);
    end if;
    xpm := Moving_Flag(start_mf,xp);
    if pivot /= 0 then
      put_line(file,"The moving coordinates after multiplication by M :");
      QuadDobl_Complex_Poly_Matrices_io.put(file,xpm);
    end if;
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if pivot = 0 then
      Track_Next_Move(file,dim,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
    else -- pivot /= 0, reconditioned homotopy has one extra variable
      Setup_Flag_Homotopies.Append(gh,rlq);
      put_line(file,"The reconditioned swap homotopy :");
      put_line(file,gh.all);
      Track_Next_Move(file,dim+1,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
      declare
        rsols : constant QuadDobl_Complex_Solutions.Solution_List
              := Recondition_Swap_Homotopies.Rescale_Solutions
                   (sols,s,locmap,xp,pivot);
        use QuadDobl_Complex_Solutions;
      begin
        QuadDobl_Complex_Solutions.Clear(sols);
        sols := rsols;
        put_line(file,"The rescaled solutions :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
        Recondition_Swap_Homotopies.Remove_One_Variable(xpm,dim+1);
      end;
    end if;
    if not fail then
      put(file,"Transforming solution planes with critical row = ");
      put(file,ctr,1); put_line(file,".");
      if big_r > ctr + 1
       then Checker_Homotopies.First_Swap_Coordinates
              (file,n,k,ctr,big_r,dc,s,q,p,qr,qc,pr,pc,mf,xpm,sols);
       else Checker_Homotopies.Second_Swap_Coordinates
              (file,n,k,ctr,s,q,qr,qc,mf,xpm,sols);
      end if;
      if verify then
        put_line(file,"Verifying after coordinate changes ...");
        Verify_Intersection_Conditions
          (file,n,k,q,qr,qc,minrep,cond,mf,vf,sols,tol,fail);
      end if;
    end if;
    QuadDobl_Complex_Poly_Systems.Clear(gh);
    QuadDobl_Complex_Poly_Matrices.Clear(xp);
    QuadDobl_Complex_Poly_Matrices.Clear(xpm);
  exception
    when others => put_line("exception in Swap_Homotopy"); raise;
  end Swap_Homotopy;

  procedure Swap_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    big_r : constant integer32 := Checker_Homotopies.Swap_Checker(q,qr,qc);
    dc : constant integer32 := Checker_Moves.Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    s : constant integer32 := Checker_Homotopies.Swap_Column(ctr,locmap);
    xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Swap_Moving_Plane(n,k,ctr,big_r,s,q,p,pr,pc);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlq : Standard_Complex_Polynomials.Poly;
    pivot : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Swap_Homotopy 7 ...");
    end if;
    if Checker_Homotopies.Is_Zone_A_Empty(locmap,p,ctr,s,dc)
     then pivot := 0; -- no reconditioning needed
     else Recondition_Swap_Homotopy(dim,ctr,s,locmap,xp,sols,rlq,pivot);
    end if;
    xpm := Moving_Flag(start_mf,xp);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if pivot = 0 then
      Track_Next_Move(dim,nt,gh.all,tosqr,tol,sols,fail,rpt,vrblvl-1);
    else -- pivot /= 0, reconditioned homotopy has one extra variable
      Setup_Flag_Homotopies.Append(gh,rlq);
      Track_Next_Move(dim+1,nt,gh.all,tosqr,tol,sols,fail,rpt,vrblvl-1);
      declare
        rsols : constant Standard_Complex_Solutions.Solution_List
              := Recondition_Swap_Homotopies.Rescale_Solutions
                   (sols,s,locmap,xp,pivot);
      begin
        Standard_Complex_Solutions.Clear(sols);
        sols := rsols;
        Recondition_Swap_Homotopies.Remove_One_Variable(xpm,dim+1);
      end;
    end if;
    if not fail then
      if big_r > ctr + 1
       then Checker_Homotopies.First_Swap_Coordinates
              (n,k,ctr,big_r,dc,s,q,p,qr,qc,pr,pc,mf,xpm,sols);
       else Checker_Homotopies.Second_Swap_Coordinates
              (n,k,ctr,s,q,qr,qc,mf,xpm,sols);
      end if;
    end if;
    Standard_Complex_Poly_Systems.Clear(gh);
    Standard_Complex_Poly_Matrices.Clear(xp);
    Standard_Complex_Poly_Matrices.Clear(xpm);
  exception
    when others => put_line("exception in Swap_Homotopy"); raise;
  end Swap_Homotopy;

  procedure Swap_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    big_r : constant integer32 := Checker_Homotopies.Swap_Checker(q,qr,qc);
    dc : constant integer32 := Checker_Moves.Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    s : constant integer32 := Checker_Homotopies.Swap_Column(ctr,locmap);
    xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Swap_Moving_Plane(n,k,ctr,big_r,s,q,p,pr,pc);
    xpm : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlq : DoblDobl_Complex_Polynomials.Poly;
    pivot : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Swap_Homotopy 8 ...");
    end if;
    if Checker_Homotopies.Is_Zone_A_Empty(locmap,p,ctr,s,dc)
     then pivot := 0; -- no reconditioning needed
     else Recondition_Swap_Homotopy(dim,ctr,s,locmap,xp,sols,rlq,pivot);
    end if;
    xpm := Moving_Flag(start_mf,xp);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if pivot = 0 then
      Track_Next_Move(dim,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
    else -- pivot /= 0, reconditioned homotopy has one extra variable
      Setup_Flag_Homotopies.Append(gh,rlq);
      Track_Next_Move(dim+1,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
      declare
        rsols : constant DoblDobl_Complex_Solutions.Solution_List
              := Recondition_Swap_Homotopies.Rescale_Solutions
                   (sols,s,locmap,xp,pivot);
      begin
        DoblDobl_Complex_Solutions.Clear(sols);
        sols := rsols;
        Recondition_Swap_Homotopies.Remove_One_Variable(xpm,dim+1);
      end;
    end if;
    if not fail then
      if big_r > ctr + 1
       then Checker_Homotopies.First_Swap_Coordinates
              (n,k,ctr,big_r,dc,s,q,p,qr,qc,pr,pc,mf,xpm,sols);
       else Checker_Homotopies.Second_Swap_Coordinates
              (n,k,ctr,s,q,qr,qc,mf,xpm,sols);
      end if;
    end if;
    DoblDobl_Complex_Poly_Systems.Clear(gh);
    DoblDobl_Complex_Poly_Matrices.Clear(xp);
    DoblDobl_Complex_Poly_Matrices.Clear(xpm);
  exception
    when others => put_line("exception in Swap_Homotopy"); raise;
  end Swap_Homotopy;

  procedure Swap_Homotopy
              ( n,k,ctr,nt : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep,tosqr : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean;
                vrblvl : in integer32 := 0 ) is

    big_r : constant integer32 := Checker_Homotopies.Swap_Checker(q,qr,qc);
    dc : constant integer32 := Checker_Moves.Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    s : constant integer32 := Checker_Homotopies.Swap_Column(ctr,locmap);
    xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Swap_Moving_Plane(n,k,ctr,big_r,s,q,p,pr,pc);
    xpm : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    rlq : QuadDobl_Complex_Polynomials.Poly;
    pivot : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in moving_flag_continuation.Swap_Homotopy 9 ...");
    end if;
    if Checker_Homotopies.Is_Zone_A_Empty(locmap,p,ctr,s,dc)
     then pivot := 0; -- no reconditioning needed
     else Recondition_Swap_Homotopy(dim,ctr,s,locmap,xp,sols,rlq,pivot);
    end if;
    xpm := Moving_Flag(start_mf,xp);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if pivot = 0 then
      Track_Next_Move(dim,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
    else -- pivot /= 0, reconditioned homotopy has one extra variable
      Setup_Flag_Homotopies.Append(gh,rlq);
      Track_Next_Move(dim+1,nt,gh.all,tosqr,tol,sols,fail,vrblvl-1);
      declare
        rsols : constant QuadDobl_Complex_Solutions.Solution_List
              := Recondition_Swap_Homotopies.Rescale_Solutions
                   (sols,s,locmap,xp,pivot);
      begin
        QuadDobl_Complex_Solutions.Clear(sols);
        sols := rsols;
        Recondition_Swap_Homotopies.Remove_One_Variable(xpm,dim+1);
      end;
    end if;
    if not fail then
      if big_r > ctr + 1
       then Checker_Homotopies.First_Swap_Coordinates
              (n,k,ctr,big_r,dc,s,q,p,qr,qc,pr,pc,mf,xpm,sols);
       else Checker_Homotopies.Second_Swap_Coordinates
              (n,k,ctr,s,q,qr,qc,mf,xpm,sols);
      end if;
    end if;
    QuadDobl_Complex_Poly_Systems.Clear(gh);
    QuadDobl_Complex_Poly_Matrices.Clear(xp);
    QuadDobl_Complex_Poly_Matrices.Clear(xpm);
  exception
    when others => put_line("exception in Swap_Homotopy"); raise;
  end Swap_Homotopy;

end Moving_Flag_Continuation;
