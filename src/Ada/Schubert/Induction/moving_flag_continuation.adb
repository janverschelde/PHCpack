with Communications_with_User;          use Communications_with_User;
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
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;      use Standard_Natural_Matrices_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with Symbol_Table;
with Matrix_Indeterminates;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io; use Standard_Complex_Poly_Matrices_io;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_Matrices;
with DoblDobl_Complex_Poly_Matrices_io; use DoblDobl_Complex_Poly_Matrices_io;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Poly_Matrices_io; use QuadDobl_Complex_Poly_Matrices_io;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with DoblDobl_Root_Refiners;            use DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;            use QuadDobl_Root_Refiners;
with Brackets;
with Checker_Boards_io;
with Checker_Moves;
with Checker_Localization_Patterns;
with Checker_Homotopies;
with Wrapped_Path_Trackers;             use Wrapped_Path_Trackers;
with Start_Flag_Homotopies;             use Start_Flag_Homotopies;
with Setup_Flag_Homotopies;             use Setup_Flag_Homotopies;
with Moving_Flag_Homotopies;            use Moving_Flag_Homotopies;

package body Moving_Flag_Continuation is

  procedure Track_First_Move
              ( file : in file_type; n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sol : in out Standard_Complex_Solutions.Link_to_Solution;
                fail : out boolean ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Solutions;

    x : Standard_Complex_Vectors.Vector(1..n);
    xt : Standard_Complex_Vectors.Vector(1..n+1);
    y : Standard_Complex_Vectors.Vector(h'range);
    res : double_float;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : Standard_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    Start_Solution(h,fail,x,res);
    if fail then
      put_line(file,"no start solution found...");
    else
      new_line(file);
      put(file,"Residual of the start solution : ");
      put(file,res,3); new_line(file);
      xt(x'range) := x;
      xt(xt'last) := Standard_Complex_Numbers.Create(0.0);
      yh := Eval(sh,xt);
      put_line(file,"Value of the start solution at the squared homotopy :");
      put_line(file,yh);
      sols := Create(xt);
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
      Call_Path_Trackers(file,n,sh,xt,sol);
      put(file,"Residual of the end solution : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > epsfa);
    end if;
    Standard_Complex_Poly_Systems.Clear(sh);
  exception
    when others => put_line("exception in Track_First_Move"); raise;
  end Track_First_Move;

  procedure Track_First_Move
              ( file : in file_type; n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sol : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Solutions;

    x : DoblDobl_Complex_Vectors.Vector(1..n);
    xt : DoblDobl_Complex_Vectors.Vector(1..n+1);
    y : DoblDobl_Complex_Vectors.Vector(h'range);
    res : double_double;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : DoblDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;

  begin
    Start_Solution(h,fail,x,res);
    if fail then
      put_line(file,"no start solution found...");
    else
      new_line(file);
      put(file,"Residual of the start solution : ");
      put(file,res,3); new_line(file);
      xt(x'range) := x;
      xt(xt'last) := DoblDobl_Complex_Numbers.Create(integer(0));
      yh := Eval(sh,xt);
      put_line(file,"Value of the start solution at the squared homotopy :");
      put_line(file,yh);
      sols := Create(xt);
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
      Call_Path_Trackers(file,n,sh,xt,sol);
      put(file,"Residual of the end solution : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > epsfa);
    end if;
    DoblDobl_Complex_Poly_Systems.Clear(sh);
  exception
    when others => put_line("exception in Track_First_Move"); raise;
  end Track_First_Move;

  procedure Track_First_Move
              ( file : in file_type; n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sol : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Solutions;

    x : QuadDobl_Complex_Vectors.Vector(1..n);
    xt : QuadDobl_Complex_Vectors.Vector(1..n+1);
    y : QuadDobl_Complex_Vectors.Vector(h'range);
    res : quad_double;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : QuadDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;

  begin
    Start_Solution(h,fail,x,res);
    if fail then
      put_line(file,"no start solution found...");
    else
      new_line(file);
      put(file,"Residual of the start solution : ");
      put(file,res,3); new_line(file);
      xt(x'range) := x;
      xt(xt'last) := QuadDobl_Complex_Numbers.Create(integer(0));
      yh := Eval(sh,xt);
      put_line(file,"Value of the start solution at the squared homotopy :");
      put_line(file,yh);
      sols := Create(xt);
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
      Call_Path_Trackers(file,n,sh,xt,sol);
      put(file,"Residual of the end solution : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > epsfa);
    end if;
    QuadDobl_Complex_Poly_Systems.Clear(sh);
  exception
    when others => put_line("exception in Track_First_Move"); raise;
  end Track_First_Move;

  procedure Track_Next_Move
              ( file : in file_type; n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sol : in out Standard_Complex_Solutions.Link_to_Solution;
                fail : out boolean ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Solutions;

    xt : Standard_Complex_Vectors.Vector(1..n+1);
    y : Standard_Complex_Vectors.Vector(h'range);
    res : double_float;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : Standard_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
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
    yh := Eval(sh,xt);
   -- put_line(file,"Value of the start solution at the squared homotopy :");
   -- put_line(file,yh);
    res := Max_Norm(y);
    put(file,"Residual of the start solution at the squared homotopy : ");
    put(file,res,3); new_line(file);
    fail := fail and (res > tolsing);
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      sols := Create(xt);
      sh0 := Eval(sh,Standard_Complex_Numbers.Create(0.0),n+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,deflate,false);
      Call_Path_Trackers(file,n,sh,xt,sol);
      put(file,"Residual of the end solution at original homotopy : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > tolsing);
    end if;
    Clear(sh); Clear(sh0); --Clear(sols);
  exception
    when others => put_line("exception in Track_Next_Move ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( file : in file_type; n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sol : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Solutions;

    xt : DoblDobl_Complex_Vectors.Vector(1..n+1);
    y : DoblDobl_Complex_Vectors.Vector(h'range);
    res : double_double;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : DoblDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;

  begin
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
    yh := Eval(sh,xt);
   -- put_line(file,"Value of the start solution at the squared homotopy :");
   -- put_line(file,yh);
    res := Max_Norm(y);
    put(file,"Residual of the start solution at the squared homotopy : ");
    put(file,res,3); new_line(file);
    fail := fail and (res > tolsing);
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      sols := Create(xt);
      sh0 := Eval(sh,DoblDobl_Complex_Numbers.Create(integer(0)),n+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,false);
      Call_Path_Trackers(file,n,sh,xt,sol);
      put(file,"Residual of the end solution at original homotopy : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > tolsing);
    end if;
    Clear(sh); Clear(sh0); --Clear(sols);
  exception
    when others => put_line("exception in Track_Next_Move ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( file : in file_type; n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sol : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                fail : out boolean ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Solutions;

    xt : QuadDobl_Complex_Vectors.Vector(1..n+1);
    y : QuadDobl_Complex_Vectors.Vector(h'range);
    res : quad_double;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : QuadDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    sols : Solution_List;
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;

  begin
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
    yh := Eval(sh,xt);
   -- put_line(file,"Value of the start solution at the squared homotopy :");
   -- put_line(file,yh);
    res := Max_Norm(y);
    put(file,"Residual of the start solution at the squared homotopy : ");
    put(file,res,3); new_line(file);
    fail := fail and (res > tolsing);
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      sols := Create(xt);
      sh0 := Eval(sh,QuadDobl_Complex_Numbers.Create(integer(0)),n+1);
      Reporting_Root_Refiner
        (file,sh0,sols,epsxa,epsfa,tolsing,numit,3,false);
      Call_Path_Trackers(file,n,sh,xt,sol);
      put(file,"Residual of the end solution at original homotopy : ");
      y := Eval(h,xt); res := Max_Norm(y);
      put(file,res,3); new_line(file); new_line(file);
      fail := (res > tolsing);
    end if;
    Clear(sh); Clear(sh0); --Clear(sols);
  exception
    when others => put_line("exception in Track_Next_Move ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( file : in file_type; n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Solutions;

    xt : Standard_Complex_Vectors.Vector(1..n+1);
    y : Standard_Complex_Vectors.Vector(h'range);
    res : double_float;
    sh : Poly_Sys(1..n) := Square(n,h);
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

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := Standard_Complex_Numbers.Create(0.0);
      y := Eval(h,xt);
      new_line(file);
     -- put_line(file,"Value of the start solution at the original homotopy :");
     -- put_line(file,y);
      res := Max_Norm(y);
      put(file,"Residual of the start solution at the original homotopy : ");
      put(file,res,3); new_line(file);
      fail := (res > tolsing);
      yh := Eval(sh,xt);
     -- put_line(file,"Value of the start solution at the squared homotopy :");
     -- put_line(file,yh);
      res := Max_Norm(y);
      put(file,"Residual of the start solution at the squared hootopy : ");
      put(file,res,3); new_line(file);
      fail := fail and (res > tolsing);
      Append(xtsols,xt_sols_last,Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      sh0 := Eval(sh,Standard_Complex_Numbers.Create(0.0),n+1);
      Reporting_Root_Refiner
        (file,sh0,xtsols,epsxa,epsfa,tolsing,numit,3,deflate,false);
      put(file,"Number of solutions in xtsols : ");
      put(file,Length_Of(xtsols),1); new_line(file);
      put(file,"Number of solutions in sols   : ");
      put(file,Length_Of(sols),1); new_line(file);
      Call_Path_Trackers(file,n,sh,xtsols,sols);
      tmp := xtsols;
      Clear(sh0);
      sh0 := Eval(sh,Standard_Complex_Numbers.Create(1.0),n+1);
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        put(file,"Residual of the end solution at squared homotopy : ");
        yh := Eval(sh0,ls.v); res := Max_Norm(yh);
        put(file,res,3); new_line(file);
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
    Clear(sh); Clear(sh0); --Clear(sols);
  exception
    when others => put_line("exception in Track_Next_Move ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( file : in file_type; n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Solutions;

    xt : DoblDobl_Complex_Vectors.Vector(1..n+1);
    y : DoblDobl_Complex_Vectors.Vector(h'range);
    res : double_double;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : DoblDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := DoblDobl_Complex_Numbers.Create(integer(0));
      y := Eval(h,xt);
      new_line(file);
     -- put_line(file,"Value of the start solution at the original homotopy :");
     -- put_line(file,y);
      res := Max_Norm(y);
      put(file,"Residual of the start solution at the original homotopy : ");
      put(file,res,3); new_line(file);
      fail := (res > tolsing);
      yh := Eval(sh,xt);
     -- put_line(file,"Value of the start solution at the squared homotopy :");
     -- put_line(file,yh);
      res := Max_Norm(y);
      put(file,"Residual of the start solution at the squared homotopy : ");
      put(file,res,3); new_line(file);
      fail := fail and (res > tolsing);
      Append(xtsols,xt_sols_last,Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      sh0 := Eval(sh,DoblDobl_Complex_Numbers.Create(integer(0)),n+1);
      Reporting_Root_Refiner
        (file,sh0,xtsols,epsxa,epsfa,tolsing,numit,3,false);
      put(file,"Number of solutions in xtsols : ");
      put(file,Length_Of(xtsols),1); new_line(file);
      put(file,"Number of solutions in sols   : ");
      put(file,Length_Of(sols),1); new_line(file);
      Call_Path_Trackers(file,n,sh,xtsols,sols);
      tmp := xtsols;
      Clear(sh0);
      sh0 := Eval(sh,DoblDobl_Complex_Numbers.Create(integer(1)),n+1);
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        put(file,"Residual of the end solution at squared homotopy : ");
        yh := Eval(sh0,ls.v); res := Max_Norm(yh);
        put(file,res,3); new_line(file);
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
    Clear(sh); Clear(sh0); --Clear(sols);
  exception
    when others => put_line("exception in Track_Next_Move ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( file : in file_type; n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Solutions;

    xt : QuadDobl_Complex_Vectors.Vector(1..n+1);
    y : QuadDobl_Complex_Vectors.Vector(h'range);
    res : quad_double;
    sh : Poly_Sys(1..n) := Square(n,h);
    yh : QuadDobl_Complex_Vectors.Vector(sh'range);
    sh0 : Poly_Sys(sh'range);
    epsxa : constant double_float := 1.0E-12;
    tolsing : constant double_float := tol;
    epsfa : constant double_float := 1.0E-12;
    numit : natural32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := QuadDobl_Complex_Numbers.Create(integer(0));
      y := Eval(h,xt);
      new_line(file);
     -- put_line(file,"Value of the start solution at the original homotopy :");
     -- put_line(file,y);
      res := Max_Norm(y);
      put(file,"Residual of the start solution at the original homotopy : ");
      put(file,res,3); new_line(file);
      fail := (res > tolsing);
      yh := Eval(sh,xt);
     -- put_line(file,"Value of the start solution at the squared homotopy :");
     -- put_line(file,yh);
      res := Max_Norm(y);
      put(file,"Residual of the start solution at the squared homotopy : ");
      put(file,res,3); new_line(file);
      fail := fail and (res > tolsing);
      Append(xtsols,xt_sols_last,Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if fail then
      put_line(file,"-> residual too high, abort path tracking");
    else
      sh0 := Eval(sh,QuadDobl_Complex_Numbers.Create(integer(0)),n+1);
      Reporting_Root_Refiner
        (file,sh0,xtsols,epsxa,epsfa,tolsing,numit,3,false);
      put(file,"Number of solutions in xtsols : ");
      put(file,Length_Of(xtsols),1); new_line(file);
      put(file,"Number of solutions in sols   : ");
      put(file,Length_Of(sols),1); new_line(file);
      Call_Path_Trackers(file,n,sh,xtsols,sols);
      tmp := xtsols;
      Clear(sh0);
      sh0 := Eval(sh,QuadDobl_Complex_Numbers.Create(integer(1)),n+1);
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        put(file,"Residual of the end solution at squared homotopy : ");
        yh := Eval(sh0,ls.v); res := Max_Norm(yh);
        put(file,res,3); new_line(file);
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
    Clear(sh); Clear(sh0); --Clear(sols);
  exception
    when others => put_line("exception in Track_Next_Move ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    xt : Standard_Complex_Vectors.Vector(1..n+1);
    sh : Poly_Sys(1..n) := Square(n,h);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;

  begin
    fail := false;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := Standard_Complex_Numbers.Create(0.0);
      fail := fail and (ls.res > tol);
      Append(xtsols,xt_sols_last,Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if not fail then
      Call_Path_Trackers(n,sh,xtsols,sols);
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
    Clear(sh);
  exception
    when others => put_line("exception in Track_Next_Move ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    xt : DoblDobl_Complex_Vectors.Vector(1..n+1);
    sh : Poly_Sys(1..n) := Square(n,h);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;

  begin
    fail := false;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := DoblDobl_Complex_Numbers.Create(integer(0));
      fail := fail and (ls.res > tol);
      Append(xtsols,xt_sols_last,Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if not fail then
      Call_Path_Trackers(n,sh,xtsols,sols);
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
    Clear(sh);
  exception
    when others => put_line("exception in Track_Next_Move ..."); raise;
  end Track_Next_Move;

  procedure Track_Next_Move
              ( n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                tol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    xt : QuadDobl_Complex_Vectors.Vector(1..n+1);
    sh : Poly_Sys(1..n) := Square(n,h);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    xtsols,xt_sols_last : Solution_List;

  begin
    fail := false;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      xt(ls.v'range) := ls.v;
      xt(xt'last) := QuadDobl_Complex_Numbers.Create(integer(0));
      fail := fail and (ls.res > tol);
      Append(xtsols,xt_sols_last,Create(xt));
      tmp := Tail_Of(tmp);
    end loop;
    if not fail then
      Call_Path_Trackers(n,sh,xtsols,sols);
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
    Clear(sh);
  exception
    when others => put_line("exception in Track_Next_Move ..."); raise;
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
                fail : out boolean ) is

    use Standard_Complex_Solutions;

    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(q,qr,qc));
    x : Standard_Complex_Vectors.Vector(1..dim);
    res : double_float;

  begin
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
                fail : out boolean ) is

    use DoblDobl_Complex_Solutions;

    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(q,qr,qc));
    x : DoblDobl_Complex_Vectors.Vector(1..dim);
    res : double_double;

  begin
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
                fail : out boolean ) is

    use QuadDobl_Complex_Solutions;

    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(q,qr,qc));
    x : QuadDobl_Complex_Vectors.Vector(1..dim);
    res : quad_double;

  begin
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
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    use Standard_Complex_Solutions;

  begin
    put(file,"Transforming solution planes with critical row = ");
    put(file,ctr,1); put_line(file,".");
    put_line(file,"The solution given to the Trivial_Stay_Coordinates : ");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    Checker_Homotopies.Trivial_Stay_Coordinates
      (file,n,k,ctr,q,p,qr,qc,pr,pc,sols);
    put_line(file,"Verifying after coordinate changes ...");
    Verify_Intersection_Conditions
      (file,n,k,q,qr,qc,minrep,cond,mf,vf,sols,tol,fail);
  end Trivial_Stay;

  procedure Trivial_Stay
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    use DoblDobl_Complex_Solutions;

  begin
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
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    use QuadDobl_Complex_Solutions;

  begin
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
              ( n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                fail : out boolean ) is
  begin
    fail := false; -- no checks anymore ...
    Checker_Homotopies.Trivial_Stay_Coordinates
      (n,k,ctr,q,p,qr,qc,pr,pc,sols);
  end Trivial_Stay;

  procedure Trivial_Stay
              ( n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is
  begin
    fail := false; -- no checks anymore ...
    Checker_Homotopies.Trivial_Stay_Coordinates
      (n,k,ctr,q,p,qr,qc,pr,pc,sols);
  end Trivial_Stay;

  procedure Trivial_Stay
              ( n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is
  begin
    fail := false; -- no checks anymore ...
    Checker_Homotopies.Trivial_Stay_Coordinates
      (n,k,ctr,q,p,qr,qc,pr,pc,sols);
  end Trivial_Stay;

  procedure Stay_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                ls : in out Standard_Complex_Solutions.Link_to_Solution; 
                tol : in double_float; fail : out boolean ) is

    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
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
     then Track_First_Move(file,dim,gh.all,tol,ls,fail);
     else Track_Next_Move(file,dim,gh.all,tol,ls,fail);
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
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution; 
                tol : in double_float; fail : out boolean ) is

    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
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
     then Track_First_Move(file,dim,gh.all,tol,ls,fail);
     else Track_Next_Move(file,dim,gh.all,tol,ls,fail);
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
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution; 
                tol : in double_float; fail : out boolean ) is

    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
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
     then Track_First_Move(file,dim,gh.all,tol,ls,fail);
     else Track_Next_Move(file,dim,gh.all,tol,ls,fail);
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
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
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
    Track_Next_Move(file,dim,gh.all,tol,sols,fail);
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
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
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
    Track_Next_Move(file,dim,gh.all,tol,sols,fail);
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
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
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
    Track_Next_Move(file,dim,gh.all,tol,sols,fail);
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
              ( n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in Standard_Complex_VecMats.VecMat;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    gh : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    fail := true;
    xpm := Moving_Flag(start_mf,xp);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(dim,gh.all,tol,sols,fail);
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
              ( n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    gh : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    fail := true;
    xpm := Moving_Flag(start_mf,xp);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(dim,gh.all,tol,sols,fail);
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
              ( n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

    gh : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    xp : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k)
       := Checker_Homotopies.Stay_Moving_Plane(n,k,ctr,p,pr,pc);
    xpm : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));

  begin
    fail := true;
    xpm := Moving_Flag(start_mf,xp);
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(dim,gh.all,tol,sols,fail);
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

  procedure Swap_Homotopy
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                ls : in out Standard_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean ) is

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

  begin
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if ind = 0
     then Track_First_Move(file,dim,gh.all,tol,ls,fail);
     else Track_Next_Move(file,dim,gh.all,tol,ls,fail);
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
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean ) is

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

  begin
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if ind = 0
     then Track_First_Move(file,dim,gh.all,tol,ls,fail);
     else Track_Next_Move(file,dim,gh.all,tol,ls,fail);
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
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                ls : in out QuadDobl_Complex_Solutions.Link_to_Solution;
                tol : in double_float; fail : out boolean ) is

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

  begin
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    if ind = 0
     then Track_First_Move(file,dim,gh.all,tol,ls,fail);
     else Track_Next_Move(file,dim,gh.all,tol,ls,fail);
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
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

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

  begin
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(file,dim,gh.all,tol,sols,fail);
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
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

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

  begin
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(file,dim,gh.all,tol,sols,fail);
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
              ( file : in file_type; n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                verify,minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

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

  begin
    Initialize_Homotopy_Symbols(natural32(dim),locmap);
   -- put_line(file,"The moving coordinates : "); put(file,xp); 
    xpm := Moving_Flag(start_mf,xp);
   -- put_line(file,"The moving coordinates after multiplication by M :");
   -- put(file,xpm);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(file,dim,gh.all,tol,sols,fail);
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
              ( n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in Standard_Complex_Matrices.Matrix;
                vf : in Standard_Complex_VecMats.VecMat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

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

  begin
    xpm := Moving_Flag(start_mf,xp);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(dim,gh.all,tol,sols,fail);
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
              ( n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in DoblDobl_Complex_Matrices.Matrix;
                vf : in DoblDobl_Complex_VecMats.VecMat;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

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

  begin
    xpm := Moving_Flag(start_mf,xp);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(dim,gh.all,tol,sols,fail);
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
              ( n,k,ctr,ind : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                minrep : in boolean;
                cond : in Standard_Natural_VecVecs.VecVec;
                mf,start_mf : in QuadDobl_Complex_Matrices.Matrix;
                vf : in QuadDobl_Complex_VecMats.VecMat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                tol : in double_float; fail : out boolean ) is

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

  begin
    xpm := Moving_Flag(start_mf,xp);
    fail := true;
    if minrep
     then Minimal_Flag_Conditions(n,k,xpm,cond,vf,gh);
     else Flag_Conditions(n,k,xpm,cond,vf,gh);
    end if;
    Track_Next_Move(dim,gh.all,tol,sols,fail);
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
