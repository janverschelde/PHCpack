with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Solution_Splitters;        use Standard_Solution_Splitters;
with Standard_Scaling;                   use Standard_Scaling;
with Black_Box_Root_Counters;            use Black_Box_Root_Counters;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;

package body Witness_Generate_and_Classify is

-- AUXILIARIES :

  procedure Down_Continuation
              ( file : in file_type;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                level : in integer32; sols : in out Solution_List;
                pocotime : out duration ) is

  -- DESCRIPTION :
  --   Performs a continuation to remove the slice from the embedded system.
  --   On entry, sols contains the start solutions, on return, the
  --   computed solutions are in the list sols.

    target : constant Standard_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(0.0));
    Black_Box_Polynomial_Continuation(file,true,target,embsys,sols,pocotime);
  end Down_Continuation;

  procedure Timing_Summary ( file : in file_type;
                             roco,hoco,poco,total : in duration ) is

  -- DESCRIPTION :
  --   Displays the timing summary for the black-box solver.

    b0 : constant string :=
     "  ---------------------------------------------------------------------";
    b00 : constant string :=
     "  =====================================================================";
    b1 : constant string :=
     "  |         TIMING INFORMATION SUMMARY for Black-Box Solver           |";
    b2 : constant string :=
     "  |   root counts  |  start system  |  continuation  |   total time   |";

  begin
    put_line(file,b0);
    put_line(file,b1);
    put_line(file,b0);
    put_line(file,b2);
    put_line(file,b00);
    put(file,"  | ");
    print_hms(file,roco); put(file," | ");
    print_hms(file,hoco); put(file," | ");
    print_hms(file,poco); put(file," | ");
    print_hms(file,total); put_line(file," |");
    put_line(file,b0);
  end Timing_Summary;

  procedure Black_Box_Solver
              ( file : in file_type; sys : in Poly_Sys;
                deg : in boolean; sols : out Solution_List;
                rc : out natural32; totaltime : out duration ) is

  -- DESCRIPTION :
  --   The black-box solver consists of the following stages:
  --   root-counting, construction of start system and path following.
  --   The flag "deg" indicates when true that root counting should
  --   only be based on the degrees.

    timer : Timing_Widget;
    q : Poly_Sys(sys'range);
    roco,hoco,poco,total : duration;
    qsols0 : Solution_List;

  begin
    tstart(timer);
    put_line(file,"BLACK BOX SOLVER ON : ");
    put(file,sys'last,1); new_line(file);
    put(file,sys);
    declare
      pp : Poly_Sys(sys'range);
    begin
      Copy(sys,pp);
      Black_Box_Root_Counting(file,0,pp,deg,rc,q,sols,qsols0,roco,hoco);
      if rc /= 0 then
        Scale(pp);
        Black_Box_Polynomial_Continuation(file,true,pp,q,sols,poco);
      end if;
      Clear(pp);
    end;
    tstop(timer);
    total := Elapsed_User_Time(timer);
    new_line(file);
    print_times(file,timer,"Solving the polynomial system");
    new_line(file);
    Timing_Summary(file,roco,hoco,poco,total);
    totaltime := total;
  end Black_Box_Solver;

  procedure Update_Flow_Table
              ( tab : in out Standard_Natural_Matrices.Matrix;
                i : in integer32; sols,sols0,sols1 : in Solution_List ) is

  -- DESCRIPTION :
  --   Updates the i-th row in the table with the length of the solution
  --   lists sols (all end points of the paths), sols0 (paths ended with
  --   zz = 0), and sols1 (paths ended at zz /= 0).

  begin
    tab(i,1) := Length_Of(sols);
    tab(i,2) := Length_Of(sols1);
    tab(i,3) := Length_Of(sols0);
    tab(i,4) := tab(i,1) - tab(i,2) - tab(i,3);
  end Update_Flow_Table;

  procedure Update_Flow_Table
              ( tab : in out Standard_Natural_Matrices.Matrix;
                i : in integer32; junkcnt : in natural32 ) is

  -- DESCRIPTION :
  --   Updates the flow table with the information of the classification
  --   of points as junk.

  begin
    tab(i,2) := tab(i,2) - junkcnt;
    tab(i,3) := tab(i,3) + junkcnt;
  end Update_Flow_Table;
                
-- TARGET PROCEDURES :

  procedure Witness_Generate
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 k : in integer32; degroco : in boolean;
                 zerotol : in double_float;
                 gentims : out Array_of_Duration;
                 flowtab : out Standard_Natural_Matrices.Matrix;
                 dc : out Standard_Irreducible_Decomposition ) is

    ep : constant Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..k)
       := Slice_and_Embed(p,natural32(k));
    sols,sols0,sols1 : Solution_List;
    rc : natural32;
    n : constant integer32 := p'last;

  begin
    Add_Embed_Symbols(natural32(k));
    gentims := (gentims'range => 0.0);
    Black_Box_Solver(file,ep(k).all,degroco,sols,rc,gentims(integer(k)));
    dc := Create(ep);
    Filter_and_Split_Solutions(file,sols,n,k,zerotol,sols0,sols1);
    Update_Flow_Table(flowtab,k,sols,sols0,sols1);
    if not Is_Null(sols0) then
      declare
        cpsols : Solution_List;
      begin
        Copy(sols0,cpsols);
        Add_Generic_Points(dc,k,cpsols);
      end;
    end if;
    if not Is_Null(sols1) then
      Copy(sols1,sols);
      for i in reverse 1..k loop
        Down_Continuation(file,ep(i).all,i,sols,gentims(integer(i)-1));
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions(file,sols,n,i-1,zerotol,sols0,sols1);
        Update_Flow_Table(flowtab,i-1,sols,sols0,sols1);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Add_Generic_Points(dc,0,rsols1);
            end;
          end if;
        else 
          if not Is_Null(sols0) then
            declare
              rsols0 : constant Solution_List := Remove_Component(sols0);
            begin
              Add_Generic_Points(dc,i-1,rsols0);
            end;
          end if;
        end if;
        Clear(sols);
        exit when Is_Null(sols1);
        sols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
  end Witness_Generate;

  procedure Witness_Classify
               ( file : in file_type; full_output : in boolean;
                 dc : in out Standard_Irreducible_Decomposition;
                 method : in natural32; stoptol,membtol : in double_float;
                 clatims : out Array_of_Duration;
                 fp,fp_last : in out List ) is

    k : constant integer32 := Top_Dimension(dc);
    threshold : constant natural32 := natural32(stoptol);
    timer : Timing_Widget;
    cnt : natural32;

  begin
    tstart(timer);
    if method < 4
     then Breakup(file,full_output,dc,k,method,stoptol,membtol,fp,fp_last);
     else Monodromy_Breakup(file,dc,k,threshold,membtol,fp,fp_last);
    end if;
    tstop(timer);
    clatims(integer(k)) := Elapsed_User_Time(timer);
    for i in reverse 1..(k-1) loop
      tstart(timer);
      if method < 4 then
        Filter(file,dc,i,membtol,fp,fp_last,cnt);
        Breakup(file,full_output,dc,i,method,stoptol,membtol,fp,fp_last);
      else
        Homotopy_Filter(file,dc,i,membtol,fp,fp_last,cnt);
        Monodromy_Breakup(file,dc,i,threshold,membtol,fp,fp_last);
      end if;
      tstop(timer);
      clatims(integer(i)) := Elapsed_User_Time(timer);
    end loop;
    tstart(timer);
    if method < 4
     then Filter(file,dc,0,membtol,fp,fp_last,cnt);
     else Homotopy_Filter(file,dc,0,membtol,fp,fp_last,cnt);
    end if;
    tstop(timer);
    clatims(0) := Elapsed_User_Time(timer);
  end Witness_Classify;

  procedure Witness_Classify
              ( file : in file_type; full_output : in boolean;
                dc : in out Multprec_Irreducible_Decomposition;
                method,size : in natural32; stoptol,membtol : in double_float;
                clatims : out Array_of_Duration;
                fp,fp_last : in out List ) is

    k : constant integer32 := Top_Dimension(dc);
    timer : Timing_Widget;
    cnt : natural32;

  begin
    tstart(timer);
    Breakup(file,full_output,dc,k,method,size,stoptol,membtol,fp,fp_last);
    tstop(timer);
    clatims(integer(k)) := Elapsed_User_Time(timer);
    for i in reverse 1..(k-1) loop
      tstart(timer);
      Filter(file,dc,i,membtol,fp,fp_last,cnt);
      Breakup(file,full_output,dc,i,method,size,stoptol,membtol,fp,fp_last);
      tstop(timer);
      clatims(integer(i)) := Elapsed_User_Time(timer);
    end loop;
    tstart(timer);
    Filter(file,dc,0,membtol,fp,fp_last,cnt);
    tstop(timer);
    clatims(0) := Elapsed_User_Time(timer);
  end Witness_Classify;

  procedure Witness_Generate_Classify
              ( file : in file_type; full_output : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                k : in integer32; method : in natural32;
                degroco : in boolean;
                zerotol,stoptol,membtol : in double_float;
                gentims,clatims : out Array_of_Duration;
                flowtab : out Standard_Natural_Matrices.Matrix;
                dc : out Standard_Irreducible_Decomposition;
                fp,fp_last : in out List ) is

    ep : constant Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..k)
       := Slice_and_Embed(p,natural32(k));
    sols,sols0,sols1 : Solution_List;
    timer : Timing_Widget;
    rc,cnt : natural32;
    n : constant integer32 := p'last;

  begin
    Add_Embed_Symbols(natural32(k));
    gentims := (gentims'range => 0.0);
    Black_Box_Solver(file,ep(k).all,degroco,sols,rc,gentims(integer(k)));
    dc := Create(ep);
    Filter_and_Split_Solutions(file,sols,n,k,zerotol,sols0,sols1);
    Update_Flow_Table(flowtab,k,sols,sols0,sols1);
    if not Is_Null(sols0) then
      declare
        cpsols : Solution_List;
      begin
        Copy(sols0,cpsols);
        Add_Generic_Points(dc,k,cpsols);
      end;
      tstart(timer);
      Breakup(file,full_output,dc,k,method,stoptol,membtol,fp,fp_last);
      tstop(timer);
      clatims(integer(k)) := Elapsed_User_Time(timer);
    end if;
    if not Is_Null(sols1) then
      Copy(sols1,sols);
      for i in reverse 1..k loop
        Down_Continuation(file,ep(i).all,i,sols,gentims(integer(i)-1));
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions(file,sols,n,i-1,zerotol,sols0,sols1);
        Update_Flow_Table(flowtab,i-1,sols,sols0,sols1);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Add_Generic_Points(dc,0,rsols1);
            end;
            tstart(timer);
            Filter(file,dc,0,membtol,fp,fp_last,cnt);
            tstop(timer);
            clatims(0) := Elapsed_User_Time(timer);
            Update_Flow_Table(flowtab,0,cnt);
          end if;
        else 
          if not Is_Null(sols0) then
            declare
              rsols0 : constant Solution_List := Remove_Component(sols0);
            begin
              Add_Generic_Points(dc,i-1,rsols0);
            end;
            tstart(timer);
            Filter(file,dc,i-1,membtol,fp,fp_last,cnt);
            Breakup(file,full_output,dc,i-1,method,stoptol,membtol,fp,fp_last);
            tstop(timer);
            clatims(integer(i)-1) := Elapsed_User_Time(timer);
            Update_Flow_Table(flowtab,i-1,cnt);
          end if;
        end if;
        Clear(sols);
        exit when Is_Null(sols1);
        sols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
  end Witness_Generate_Classify;

  procedure Witness_Generate_Classify
              ( file : in file_type; full_output : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mp : in Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                k : in integer32; method,size : in natural32;
                degroco : in boolean;
                zerotol,stoptol,membtol : in double_float;
                gentims,clatims : out Array_of_Duration;
                flowtab : out Standard_Natural_Matrices.Matrix;
                dc : out Multprec_Irreducible_Decomposition;
                fp,fp_last : in out List ) is

    ep : constant Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..k)
       := Slice_and_Embed(p,natural32(k));
    sols,sols0,sols1 : Solution_List;
    timer : Timing_Widget;
    rc,cnt : natural32;
    n : constant integer32 := p'last;

  begin
    Add_Embed_Symbols(natural32(k));
    gentims := (gentims'range => 0.0);
    Black_Box_Solver(file,ep(k).all,degroco,sols,rc,gentims(integer(k)));
    dc := Create(ep);
    Add_Original(dc,mp);
    Filter_and_Split_Solutions(file,sols,n,k,zerotol,sols0,sols1);
    Update_Flow_Table(flowtab,k,sols,sols0,sols1);
    if not Is_Null(sols0) then
      declare
        cpsols : Solution_List;
      begin
        Copy(sols0,cpsols);
        Add_Generic_Points(dc,k,cpsols);
      end;
      tstart(timer);
      Breakup(file,full_output,dc,k,method,size,stoptol,membtol,fp,fp_last);
      tstop(timer);
      clatims(integer(k)) := Elapsed_User_Time(timer);
    end if;
    if not Is_Null(sols1) then
      Copy(sols1,sols);
      for i in reverse 1..k loop
        Down_Continuation(file,ep(i).all,i,sols,gentims(integer(i)-1));
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions(file,sols,n,i-1,zerotol,sols0,sols1);
        Update_Flow_Table(flowtab,i-1,sols,sols0,sols1);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Add_Generic_Points(dc,0,rsols1);
            end;
            tstart(timer);
            Filter(file,dc,0,membtol,fp,fp_last,cnt);
            tstop(timer);
            clatims(0) := Elapsed_User_Time(timer);
            Update_Flow_Table(flowtab,0,cnt);
          end if;
        else
          if not Is_Null(sols0) then
            declare
              rsols0 : constant Solution_List := Remove_Component(sols0);
            begin
              Add_Generic_Points(dc,i-1,rsols0);
            end;
            tstart(timer);
            Filter(file,dc,i-1,membtol,fp,fp_last,cnt);
            Breakup(file,full_output,dc,i-1,method,size,
                    stoptol,membtol,fp,fp_last);
            tstop(timer);
            clatims(integer(i)-1) := Elapsed_User_Time(timer);
            Update_Flow_Table(flowtab,i-1,cnt);
          end if;
        end if;
        Clear(sols);
        exit when Is_Null(sols1);
        sols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
  end Witness_Generate_Classify;

end Witness_Generate_and_Classify;
