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

    procedure Track is
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(h,n+1);
    Track(sols,false,Standard_Complex_Numbers.Create(1.0));
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    Standard_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
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

    procedure Track is
      new Silent_Continue(Max_Norm,DoblDobl_Homotopy.Eval,
                          DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

  begin
    DoblDobl_Homotopy.Create(h,n+1);
    Track(sols,DoblDobl_Complex_Numbers.Create(one));
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    DoblDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
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

    procedure Track is
      new Silent_Continue(Max_Norm,QuadDobl_Homotopy.Eval,
                          QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    QuadDobl_Homotopy.Create(h,n+1);
    Track(sols,QuadDobl_Complex_Numbers.Create(one));
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    QuadDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
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

    procedure Track is
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                             Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(h,n+1);
    Track(file,sols,false,Standard_Complex_Numbers.Create(1.0));
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    Standard_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
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

    procedure Track is
      new Reporting_Continue(Max_Norm,DoblDobl_Homotopy.Eval,
                             DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

  begin
    DoblDobl_Homotopy.Create(h,n+1);
    Track(file,sols,DoblDobl_Complex_Numbers.Create(one));
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    DoblDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
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

    procedure Track is
      new Reporting_Continue(Max_Norm,QuadDobl_Homotopy.Eval,
                             QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    QuadDobl_Homotopy.Create(h,n+1);
    Track(file,sols,QuadDobl_Complex_Numbers.Create(one));
    xt(xt'first..xt'last-1) := Head_Of(sols).v;
    xt(xt'last) := Head_Of(sols).t;
    sol := Head_Of(sols);
    QuadDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
  end Call_Path_Trackers;

-- TRACKING MANY PATHS WITHOUT OUTPUT :

  procedure Call_Path_Trackers
              ( n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out Standard_Complex_Solutions.Solution_List;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
    use Standard_IncFix_Continuation;

    xtp,tmp : Solution_List;
    xtls,ls : Link_to_Solution;

    procedure Track is
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(h,n+1);
    Track(xtsols,false,Standard_Complex_Numbers.Create(1.0));
    tmp := sols;
    xtp := xtsols;
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
    Standard_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_IncFix_Continuation;

    xtp,tmp : Solution_List;
    xtls,ls : Link_to_Solution;
    one : constant double_double := create(1.0);

    procedure Track is
      new Silent_Continue(Max_Norm,DoblDobl_Homotopy.Eval,
                          DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

  begin
    DoblDobl_Homotopy.Create(h,n+1);
    Track(xtsols,DoblDobl_Complex_Numbers.Create(one));
    tmp := sols;
    xtp := xtsols;
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
    DoblDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_IncFix_Continuation;

    xtp,tmp : Solution_List;
    xtls,ls : Link_to_Solution;
    one : constant quad_double := create(1.0);

    procedure Track is
      new Silent_Continue(Max_Norm,QuadDobl_Homotopy.Eval,
                          QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    QuadDobl_Homotopy.Create(h,n+1);
    Track(xtsols,QuadDobl_Complex_Numbers.Create(one));
    tmp := sols;
    xtp := xtsols;
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
    QuadDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
  end Call_Path_Trackers;

-- TRACKING MANY PATHS WITH OUTPUT TO FILE :

  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out Standard_Complex_Solutions.Solution_List;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
    use Standard_IncFix_Continuation;

    xtp,tmp : Solution_List;
    xtls,ls : Link_to_Solution;

    procedure Track is
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                             Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(h,n+1);
    Track(file,xtsols,false,Standard_Complex_Numbers.Create(1.0));
    tmp := sols;
    xtp := xtsols;
   -- put_line(file,"In Call_Path_Trackers ...");
    put(file,"Number of solutions in sols   : ");
    put(file,Length_Of(sols),1); new_line(file);
    put(file,"Number of solutions in xtsols : ");
    put(file,Length_Of(xtsols),1); new_line(file);
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
    Standard_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_IncFix_Continuation;

    xtp,tmp : Solution_List;
    xtls,ls : Link_to_Solution;
    one : constant double_double := create(1.0);

    procedure Track is
      new Reporting_Continue(Max_Norm,DoblDobl_Homotopy.Eval,
                             DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

  begin
    DoblDobl_Homotopy.Create(h,n+1);
    Track(file,xtsols,DoblDobl_Complex_Numbers.Create(one));
    tmp := sols;
    xtp := xtsols;
   -- put_line(file,"In Call_Path_Trackers ...");
    put(file,"Number of solutions in sols   : ");
    put(file,Length_Of(sols),1); new_line(file);
    put(file,"Number of solutions in xtsols : ");
    put(file,Length_Of(xtsols),1); new_line(file);
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
    DoblDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
  end Call_Path_Trackers;

  procedure Call_Path_Trackers
              ( file : in file_type; n : in integer32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                xtsols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_IncFix_Continuation;

    xtp,tmp : Solution_List;
    xtls,ls : Link_to_Solution;
    one : constant quad_double := create(1.0);

    procedure Track is
      new Reporting_Continue(Max_Norm,QuadDobl_Homotopy.Eval,
                             QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    QuadDobl_Homotopy.Create(h,n+1);
    Track(file,xtsols,QuadDobl_Complex_Numbers.Create(one));
    tmp := sols;
    xtp := xtsols;
   -- put_line(file,"In Call_Path_Trackers ...");
    put(file,"Number of solutions in sols   : ");
    put(file,Length_Of(sols),1); new_line(file);
    put(file,"Number of solutions in xtsols : ");
    put(file,Length_Of(xtsols),1); new_line(file);
    while not Is_Null(xtp) loop
      xtls := Head_Of(xtp);
      ls := Head_Of(tmp);
      ls.v := xtls.v(ls.v'range);
      ls.t := xtls.t;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      xtp := Tail_Of(xtp);
    end loop;
    QuadDobl_Homotopy.Clear;
  exception -- adding this exception handler caused no longer exception ...
    when others => put_line("exception in Call Path Trackers"); raise;
  end Call_Path_Trackers;

end Wrapped_Path_Trackers;
