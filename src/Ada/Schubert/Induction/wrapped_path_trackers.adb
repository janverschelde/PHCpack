with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_IncFix_Continuation;
with DoblDobl_IncFix_Continuation;
with QuadDobl_IncFix_Continuation;
with Drivers_for_Poly_Continuation;
with Multitasking_Continuation;         use Multitasking_Continuation;

package body Wrapped_Path_Trackers is

  function Create ( x : Standard_Complex_Vectors.Vector )
                  return Standard_Complex_Solutions.Solution is

    res : Standard_Complex_Solutions.Solution(x'last-1);

  begin
    res.t := x(x'last);
    res.m := 1;
    res.v := x(x'first..x'last-1);
    res.err := 0.0;
    res.rco := 1.0;
    res.res := 0.0;
    return res;
  end Create;

  function Create ( x : DoblDobl_Complex_Vectors.Vector )
                  return DoblDobl_Complex_Solutions.Solution is

    res : DoblDobl_Complex_Solutions.Solution(x'last-1);

  begin
    res.t := x(x'last);
    res.m := 1;
    res.v := x(x'first..x'last-1);
    res.err := create(0.0);
    res.rco := create(1.0);
    res.res := create(0.0);
    return res;
  end Create;

  function Create ( x : QuadDobl_Complex_Vectors.Vector )
                  return QuadDobl_Complex_Solutions.Solution is

    res : QuadDobl_Complex_Solutions.Solution(x'last-1);

  begin
    res.t := x(x'last);
    res.m := 1;
    res.v := x(x'first..x'last-1);
    res.err := create(0.0);
    res.rco := create(1.0);
    res.res := create(0.0);
    return res;
  end Create;

  function Create ( x : Standard_Complex_Vectors.Vector ) 
                  return Standard_Complex_Solutions.Solution_List is

    sol : constant Standard_Complex_Solutions.Solution(x'last-1) := Create(x);
    res : Standard_Complex_Solutions.Solution_List;

  begin
    Standard_Complex_Solutions.Add(res,sol);
    return res;
  end Create;

  function Create ( x : DoblDobl_Complex_Vectors.Vector ) 
                  return DoblDobl_Complex_Solutions.Solution_List is

    sol : constant DoblDobl_Complex_Solutions.Solution(x'last-1) := Create(x);
    res : DoblDobl_Complex_Solutions.Solution_List;

  begin
    DoblDobl_Complex_Solutions.Add(res,sol);
    return res;
  end Create;

  function Create ( x : QuadDobl_Complex_Vectors.Vector ) 
                  return QuadDobl_Complex_Solutions.Solution_List is

    sol : constant QuadDobl_Complex_Solutions.Solution(x'last-1) := Create(x);
    res : QuadDobl_Complex_Solutions.Solution_List;

  begin
    QuadDobl_Complex_Solutions.Add(res,sol);
    return res;
  end Create;

  procedure Set_Parameters ( file : in file_type; report : out boolean ) is

    oc : natural32;
    use Drivers_for_Poly_Continuation;

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

-- TRACKING ONE PATH WITHOUT OUTPUT :

  procedure Call_Path_Trackers
              ( n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                xt : in out Standard_Complex_Vectors.Vector;
                sol : out Standard_Complex_Solutions.Link_to_Solution ) is

    use Standard_Complex_Solutions;
    use Standard_IncFix_Continuation;

    sols : Solution_List := Create(xt);
    nbequ : constant integer32 := h'last;

    procedure Track is
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(sols,false,target=>Standard_Complex_Numbers.Create(1.0));
    else      -- overdetermined homotopy
      Track(sols,false,nbq=>nbequ,
            target=>Standard_Complex_Numbers.Create(1.0));
    end if;
    sol := Head_Of(sols);
    xt(xt'first..xt'last-1) := sol.v;
    xt(xt'last) := sol.t;
    Standard_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 1"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                xt : in out DoblDobl_Complex_Vectors.Vector;
                sol : out DoblDobl_Complex_Solutions.Link_to_Solution ) is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_IncFix_Continuation;

    sols : Solution_List := Create(xt);
    one : constant double_double := create(1.0);
    nbequ : constant integer32 := h'last;

    procedure Track is
      new Silent_Continue(Max_Norm,DoblDobl_Homotopy.Eval,
                          DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

  begin
    DoblDobl_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(sols,target=>DoblDobl_Complex_Numbers.Create(one));
    else      -- overdetermined homotopy
      Track(sols,nbq=>nbequ,target=>DoblDobl_Complex_Numbers.Create(one));
    end if;
    sol := Head_Of(sols);
    xt(xt'first..xt'last-1) := sol.v;
    xt(xt'last) := sol.t;
    DoblDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 2"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                xt : in out QuadDobl_Complex_Vectors.Vector;
                sol : out QuadDobl_Complex_Solutions.Link_to_Solution ) is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_IncFix_Continuation;

    sols : Solution_List := Create(xt);
    one : constant quad_double := create(1.0);
    nbequ : constant integer32 := h'last;

    procedure Track is
      new Silent_Continue(Max_Norm,QuadDobl_Homotopy.Eval,
                          QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    QuadDobl_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(sols,target=>QuadDobl_Complex_Numbers.Create(one));
    else      -- overdetermined homotopy
      Track(sols,nbq=>nbequ,target=>QuadDobl_Complex_Numbers.Create(one));
    end if;
    sol := Head_Of(sols);
    xt(xt'first..xt'last-1) := sol.v;
    xt(xt'last) := sol.t;
    QuadDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 3"); raise;
  end Call_Path_Trackers;

-- TRACKING ONE PATH WITH OUTPUT TO FILE :

  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                xt : in out Standard_Complex_Vectors.Vector;
                sol : out Standard_Complex_Solutions.Link_to_Solution ) is

    use Standard_Complex_Solutions;
    use Standard_IncFix_Continuation;

    sols : Solution_List := Create(xt);
    nbequ : constant integer32 := h'last;

    procedure Track is
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                             Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(file,sols,false,target=>Standard_Complex_Numbers.Create(1.0));
    else      -- overdetermined homotopy
      Track(file,sols,false,nbq=>nbequ,
            target=>Standard_Complex_Numbers.Create(1.0));
    end if;
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    Standard_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 4"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                xt : in out DoblDobl_Complex_Vectors.Vector;
                sol : out DoblDobl_Complex_Solutions.Link_to_Solution ) is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_IncFix_Continuation;

    sols : Solution_List := Create(xt);
    one : constant double_double := create(1.0);
    nbequ : constant integer32 := h'last;

    procedure Track is
      new Reporting_Continue(Max_Norm,DoblDobl_Homotopy.Eval,
                             DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

  begin
    DoblDobl_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(file,sols,target=>DoblDobl_Complex_Numbers.Create(one));
    else      -- overdetermined homotopy
      Track(file,sols,target=>DoblDobl_Complex_Numbers.Create(one));
    end if;
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    DoblDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 5"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                xt : in out QuadDobl_Complex_Vectors.Vector;
                sol : out QuadDobl_Complex_Solutions.Link_to_Solution ) is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_IncFix_Continuation;

    sols : Solution_List := Create(xt);
    one : constant quad_double := create(1.0);
    nbequ : constant integer32 := h'last;

    procedure Track is
      new Reporting_Continue(Max_Norm,QuadDobl_Homotopy.Eval,
                             QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    QuadDobl_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(file,sols,target=>QuadDobl_Complex_Numbers.Create(one));
    else      -- overdetermined homotopy
      Track(file,sols,nbq=>nbequ,target=>QuadDobl_Complex_Numbers.Create(one));
    end if;
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    QuadDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 6"); raise;
  end Call_Path_Trackers;

-- UPDATING SOLUTION LISTS :

  procedure Update ( sols : in out Standard_Complex_Solutions.Solution_List;
                     xtsols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Given in the list xtsols the coordinates of the solutions
  --   jointly with the value of continuation parameter,
  --   the list sols is updated, separating the continuation parameter t
  --   from the last coordinate of the solution vectors in xtsols.

    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;
    xtp : Solution_List := xtsols;
    ls,xtls : Link_to_Solution;

  begin
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
  end Update;

  procedure Update ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                     xtsols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Given in the list xtsols the coordinates of the solutions
  --   jointly with the value of continuation parameter,
  --   the list sols is updated, separating the continuation parameter t
  --   from the last coordinate of the solution vectors in xtsols.

    use DoblDobl_Complex_Solutions;
    tmp : Solution_List := sols;
    xtp : Solution_List := xtsols;
    ls,xtls : Link_to_Solution;

  begin
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
  end Update;

  procedure Update ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                     xtsols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Given in the list xtsols the coordinates of the solutions
  --   jointly with the value of continuation parameter,
  --   the list sols is updated, separating the continuation parameter t
  --   from the last coordinate of the solution vectors in xtsols.

    use QuadDobl_Complex_Solutions;
    tmp : Solution_List := sols;
    xtp : Solution_List := xtsols;
    ls,xtls : Link_to_Solution;

  begin
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
  end Update;

-- TRACKING MANY PATHS WITHOUT OUTPUT :

  procedure Call_Path_Trackers
              ( n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out Standard_Complex_Solutions.Solution_List;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_IncFix_Continuation;

    nbequ : constant integer32 := h'last;

    procedure Track is
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(xtsols,false,target=>Standard_Complex_Numbers.Create(1.0));
    else      -- overdetermined homotopy
      Track(xtsols,false,nbq=>nbequ,
            target=>Standard_Complex_Numbers.Create(1.0));
    end if;
    Update(sols,xtsols);
    Standard_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 7"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_IncFix_Continuation;

    one : constant double_double := create(1.0);
    nbequ : constant integer32 := h'last;

    procedure Track is
      new Silent_Continue(Max_Norm,DoblDobl_Homotopy.Eval,
                          DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

  begin
    DoblDobl_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(xtsols,target=>DoblDobl_Complex_Numbers.Create(one));
    else      -- overdetermined homotopy
      Track(xtsols,nbq=>nbequ,target=>DoblDobl_Complex_Numbers.Create(one));
    end if;
    Update(sols,xtsols);
    DoblDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 8"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_IncFix_Continuation;

    one : constant quad_double := create(1.0);
    nbequ : constant integer32 := h'last;

    procedure Track is
      new Silent_Continue(Max_Norm,QuadDobl_Homotopy.Eval,
                          QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    QuadDobl_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(xtsols,target=>QuadDobl_Complex_Numbers.Create(one));
    else      -- overdetermined homotopy
      Track(xtsols,nbq=>nbequ,target=>QuadDobl_Complex_Numbers.Create(one));
    end if;
    Update(sols,xtsols);
    QuadDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 9"); raise;
  end Call_Path_Trackers;

-- TRACKING MANY PATHS WITH OUTPUT TO FILE :

  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out Standard_Complex_Solutions.Solution_List;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
    use Standard_IncFix_Continuation;

    nbequ : constant integer32 := h'last;

    procedure Track is
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                             Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(file,xtsols,false,target=>Standard_Complex_Numbers.Create(1.0));
    else      -- overdetermined homotopy
      Track(file,xtsols,false,nbq=>nbequ,
            target=>Standard_Complex_Numbers.Create(1.0));
    end if;
   -- put_line(file,"In Call_Path_Trackers ...");
    put(file,"Number of solutions in sols   : ");
    put(file,Length_Of(sols),1); new_line(file);
    put(file,"Number of solutions in xtsols : ");
    put(file,Length_Of(xtsols),1); new_line(file);
    Update(sols,xtsols);
    Standard_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 10"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_IncFix_Continuation;

    one : constant double_double := create(1.0);
    nbequ : constant integer32 := h'last;

    procedure Track is
      new Reporting_Continue(Max_Norm,DoblDobl_Homotopy.Eval,
                             DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

  begin
    DoblDobl_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(file,xtsols,target=>DoblDobl_Complex_Numbers.Create(one));
    else      -- overdetermined homotopy
      Track(file,xtsols,nbq=>nbequ,
            target=>DoblDobl_Complex_Numbers.Create(one));
    end if;
   -- put_line(file,"In Call_Path_Trackers ...");
    put(file,"Number of solutions in sols   : ");
    put(file,Length_Of(sols),1); new_line(file);
    put(file,"Number of solutions in xtsols : ");
    put(file,Length_Of(xtsols),1); new_line(file);
    Update(sols,xtsols);
    DoblDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 11"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_IncFix_Continuation;

    one : constant quad_double := create(1.0);
    nbequ : constant integer32 := h'last;

    procedure Track is
      new Reporting_Continue(Max_Norm,QuadDobl_Homotopy.Eval,
                             QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    QuadDobl_Homotopy.Create(h,n+1);
    if nbequ = n then -- square homotopy
      Track(file,xtsols,target=>QuadDobl_Complex_Numbers.Create(one));
    else      -- overdetermined homotopy
      Track(file,xtsols,nbq=>nbequ,
            target=>QuadDobl_Complex_Numbers.Create(one));
    end if;
   -- put_line(file,"In Call_Path_Trackers ...");
    put(file,"Number of solutions in sols   : ");
    put(file,Length_Of(sols),1); new_line(file);
    put(file,"Number of solutions in xtsols : ");
    put(file,Length_Of(xtsols),1); new_line(file);
    Update(sols,xtsols);
    QuadDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 12"); raise;
  end Call_Path_Trackers;

-- MULTITASKED TRACKING OF MANY PATHS WITHOUT OUTPUT :

  procedure Multitasked_Path_Trackers
              ( nv,nt : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out Standard_Complex_Solutions.Solution_List;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    nbequ : constant integer32 := h'last;

  begin
    Standard_Homotopy.Create(h,nv+1);
    if nbequ = nv then -- square homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt);
    else      -- overdetermined homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt,nbq=>nbequ);
    end if;
    Update(sols,xtsols);
    Standard_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 7"); raise;
  end Multitasked_Path_Trackers;

  procedure Multitasked_Path_Trackers
              ( nv,nt : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    nbequ : constant integer32 := h'last;

  begin
    DoblDobl_Homotopy.Create(h,nv+1);
    if nbequ = nv then -- square homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt);
    else      -- overdetermined homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt,nbq=>nbequ);
    end if;
    Update(sols,xtsols);
    DoblDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 8"); raise;
  end Multitasked_Path_Trackers;

  procedure Multitasked_Path_Trackers
              ( nv,nt : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    nbequ : constant integer32 := h'last;

  begin
    QuadDobl_Homotopy.Create(h,nv+1);
    if nbequ = nv then -- square homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt);
    else      -- overdetermined homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt,nbq=>nbequ);
    end if;
    Update(sols,xtsols);
    QuadDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 9"); raise;
  end Multitasked_Path_Trackers;

-- MULTITASKED TRACKING OF MANY PATHS WITH OUTPUT TO FILE :

  procedure Multitasked_Path_Trackers
              ( file : in file_type; nv,nt : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out Standard_Complex_Solutions.Solution_List;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    nbequ : constant integer32 := h'last;

  begin
    Standard_Homotopy.Create(h,nv+1);
    if nbequ = nv then -- square homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt);
    else      -- overdetermined homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt,nbq=>nbequ);
    end if;
   -- put_line(file,"In Call_Path_Trackers ...");
    put(file,"Number of solutions in sols   : ");
    put(file,Length_Of(sols),1); new_line(file);
    put(file,"Number of solutions in xtsols : ");
    put(file,Length_Of(xtsols),1); new_line(file);
    Update(sols,xtsols);
    Standard_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 10"); raise;
  end Multitasked_Path_Trackers;

  procedure Multitasked_Path_Trackers
              ( file : in file_type; nv,nt : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    nbequ : constant integer32 := h'last;

  begin
    DoblDobl_Homotopy.Create(h,nv+1);
    if nbequ = nv then -- square homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt);
    else      -- overdetermined homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt,nbq=>nbequ);
    end if;
   -- put_line(file,"In Call_Path_Trackers ...");
    put(file,"Number of solutions in sols   : ");
    put(file,Length_Of(sols),1); new_line(file);
    put(file,"Number of solutions in xtsols : ");
    put(file,Length_Of(xtsols),1); new_line(file);
    Update(sols,xtsols);
    DoblDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 11"); raise;
  end Multitasked_Path_Trackers;

  procedure Multitasked_Path_Trackers
              ( file : in file_type; nv,nt : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    nbequ : constant integer32 := h'last;

  begin
    QuadDobl_Homotopy.Create(h,nv+1);
    if nbequ = nv then -- square homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt);
    else      -- overdetermined homotopy
      Silent_Multitasking_Path_Tracker(xtsols,nt,nbq=>nbequ);
    end if;
   -- put_line(file,"In Call_Path_Trackers ...");
    put(file,"Number of solutions in sols   : ");
    put(file,Length_Of(sols),1); new_line(file);
    put(file,"Number of solutions in xtsols : ");
    put(file,Length_Of(xtsols),1); new_line(file);
    Update(sols,xtsols);
    QuadDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers 12"); raise;
  end Multitasked_Path_Trackers;

end Wrapped_Path_Trackers;
