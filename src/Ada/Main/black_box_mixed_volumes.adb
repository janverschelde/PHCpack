with Timing_Package;                     use Timing_Package;
with Standard_Complex_Numbers;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Floating_Lifting_Functions;
with Apply_Induced_Permutations;         use Apply_Induced_Permutations;
with Black_Mixed_Volume_Computations;    use Black_Mixed_Volume_Computations;
with Black_Polyhedral_Continuations;     use Black_Polyhedral_Continuations;
with Root_Counters_Output;
with Pipelined_Polyhedral_Drivers;
with Black_Box_Root_Counters;

package body Black_Box_Mixed_Volumes is

  procedure Mixed_Volume
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                mivo,stmv : out natural32;
                stlb : out double_float;
                lifsup : out Link_to_Array_of_Lists;
                mix,perm,iprm : out Link_to_Vector;
                orgmcc,stbmcc : out Mixed_Subdivision;
                rocotime : out duration; verbose : in integer32 := 0 ) is

    timer : Timing_Widget;
    tmv : natural32;
    mixsub : Mixed_Subdivision;
    orgcnt,stbcnt : natural32;

  begin
    if verbose > 0 then
      put_line("-> in black_box_mixed_volumes.Mixed_Volume ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (p,mix,perm,iprm,stlb,lifsup,mixsub,orgmcc,stbmcc,
       mivo,stmv,tmv,orgcnt,stbcnt,verbose-1);
    tstop(timer);
    rocotime := Elapsed_User_Time(timer);
  end Mixed_Volume;

  procedure Mixed_Volume
              ( file : in file_type;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                mivo,stmv : out natural32;
                stlb : out double_float;
                lifsup : out Link_to_Array_of_Lists;
                mix,perm,iprm : out Link_to_Vector;
                orgmcc,stbmcc : out Mixed_Subdivision;
                rocotime : out duration; verbose : in integer32 := 0 ) is

    timer : Timing_Widget;
    tmv : natural32;
    mixsub : Mixed_Subdivision;
    orgcnt,stbcnt : natural32;

  begin
    if verbose > 0 then
      put_line("-> in black_box_mixed_volumes.Mixed_Volume ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (p,mix,perm,iprm,stlb,lifsup,mixsub,orgmcc,stbmcc,
       mivo,stmv,tmv,orgcnt,stbcnt,verbose-1);
    tstop(timer);
    Root_Counters_Output.Write_Mixed_Volumes(file,mivo,stmv);
    new_line(file);
    print_times(file,timer,"Mixed Volume Computation");
    flush(file);
    rocotime := Elapsed_User_Time(timer);
  end Mixed_Volume;

  procedure Construct_Start_System
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                mix : in Link_to_Vector;
                stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                orgmcc,stbmcc : in Mixed_Subdivision;
                q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : in out Standard_Complex_Solutions.Solution_List;
                hocotime : out duration; verbose : in integer32 := 0 ) is

    timer : timing_widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_mixed_volumes.Construct_Start_System 1 ...");
    end if;
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (0,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
  end Construct_Start_System;

  procedure Construct_Start_System
              ( file : in file_type;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                mix : in Link_to_Vector;
                stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                orgmcc,stbmcc : in Mixed_Subdivision;
                q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : in out Standard_Complex_Solutions.Solution_List;
                hocotime : out duration; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    timer : timing_widget;
    wsols : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_mixed_volumes.Construct_Start_System 2 ...");
    end if;
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (0,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
    tstop(timer);
    new_line(file);
    put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
    put_line(file,q);
    new_line(file);
    put_line(file,"START SOLUTIONS :");
    if Is_Null(qsols0) then
      put(file,Length_Of(qsols),natural32(q'last),qsols);
    else
      Push(qsols,wsols); Push(qsols0,wsols);
      put(file,Length_Of(wsols),natural32(q'last),wsols);
    end if;
    new_line(file);
    print_times(file,timer,"Polyhedral Continuation");
    flush(file);
    hocotime := Elapsed_User_Time(timer);
  end Construct_Start_System;

  procedure Black_Box_Polyhedral_Homotopies
              ( silent : in boolean;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                rc : out natural32;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                rocotime,hocotime : out duration;
                verbose : in integer32 := 0 ) is

    mv,smv : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;

  begin
    if verbose > 0 then
      put("-> in black_box_mixed_volumes.");
      put_line("Black_Box_Polyhedral_Homotopies 1 ...");
    end if;
    Mixed_Volume
      (p,mv,smv,stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    if (not silent) or (verbose > 0) then
      Root_Counters_Output.Write_Mixed_Volumes(standard_output,mv,smv);
    end if;
    if smv > mv
     then rc := smv;
     else rc := mv;
    end if;
    Construct_Start_System
      (p,mix,stlb,lifsup,orgmcc,stbmcc,q,qsols,qsols0,hocotime,verbose-1);
  end Black_Box_Polyhedral_Homotopies;

  procedure Black_Box_Polyhedral_Homotopies
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    mv,smv : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;

  begin
    if verbose > 0 then
      put("-> in black_box_mixed_volumes.");
      put_line("Black_Box_Polyhedral_Homotopies 2 ...");
    end if;
    Mixed_Volume
      (p,mv,smv,stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    if smv > mv
     then rc := smv;
     else rc := mv;
    end if;
    declare
      rcs : constant string
          := Root_Counters_Output.Mixed_Volumes_to_String(mv,smv);
    begin
      rocos := new string'(rcs);
      if verbose > 0
       then put_line(rcs);
      end if;
    end;
    Construct_Start_System
      (p,mix,stlb,lifsup,orgmcc,stbmcc,q,qsols,qsols0,hocotime,verbose-1);
  end Black_Box_Polyhedral_Homotopies;

  procedure Black_Box_Polyhedral_Homotopies
               ( nt : in natural32; silent : in boolean;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    mv,smv,tmv : natural32;
    r : integer32;
    mix,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    lq : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    timer : Timing_Widget;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

    use Standard_Complex_Solutions;

  begin
    if verbose > 0 then
      put("-> in black_box_mixed_volumes.");
      put_line("Black_Box_Polyhedral_Homotopies 3 ...");
    end if;
    stlb := Floating_Lifting_Functions.Lifting_Bound(p);
    tstart(timer);
    Pipelined_Polyhedral_Drivers.Pipelined_Polyhedral_Homotopies
      (integer32(nt),true,stlb,p,r,mix,perm,lifsup,mcc,tmv,lq,q,qsols);
    Set_Continuation_Parameter(qsols,zero);
    Apply_Induced_Permutation(p,stlb,r,perm,mix,mcc);
    Black_Box_Root_Counters.Pipelined_Stable_Continuation
        (true,r,mix,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
    Set_Continuation_Parameter(qsols0,zero);
    tstop(timer);
    if not silent
     then Root_Counters_Output.Write_Mixed_Volumes(standard_output,mv,smv);
    end if;
    rc := smv;
    rocotime := 0.0;
    hocotime := Elapsed_User_Time(timer);
  end Black_Box_Polyhedral_Homotopies;

  procedure Black_Box_Polyhedral_Homotopies
               ( nt : in natural32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    mv,smv,tmv : natural32;
    r : integer32;
    mix,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    lq : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    timer : Timing_Widget;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

    use Standard_Complex_Solutions;

  begin
    if verbose > 0 then
      put("-> in black_box_mixed_volumes.");
      put_line("Black_Box_Polyhedral_Homotopies 4 ...");
    end if;
    stlb := Floating_Lifting_Functions.Lifting_Bound(p);
    tstart(timer);
    Pipelined_Polyhedral_Drivers.Pipelined_Polyhedral_Homotopies
      (integer32(nt),true,stlb,p,r,mix,perm,lifsup,mcc,tmv,lq,q,qsols);
    Set_Continuation_Parameter(qsols,zero);
    Apply_Induced_Permutation(p,stlb,r,perm,mix,mcc);
    Black_Box_Root_Counters.Pipelined_Stable_Continuation
        (true,r,mix,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
    Set_Continuation_Parameter(qsols0,zero);
    tstop(timer);
    declare
      rcs : constant string
          := Root_Counters_Output.Mixed_Volumes_to_String(mv,smv);
    begin
      rocos := new string'(rcs);
      if verbose > 0
       then put_line(rcs);
      end if;
    end;
    rc := smv;
    rocotime := 0.0;
    hocotime := Elapsed_User_Time(timer);
  end Black_Box_Polyhedral_Homotopies;

  procedure Black_Box_Polyhedral_Homotopies
               ( file : in file_type;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    mv,smv : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;

  begin
    if verbose > 0 then
      put("-> in black_box_mixed_volumes.");
      put_line("Black_Box_Polyhedral_Homotopies 5 ...");
    end if;
    Mixed_Volume
      (file,p,mv,smv,stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,
       rocotime,verbose-1);
    if smv > mv
     then rc := smv;
     else rc := mv;
    end if;
    Construct_Start_System
      (file,p,mix,stlb,lifsup,orgmcc,stbmcc,q,qsols,qsols0,hocotime,verbose-1);
  end Black_Box_Polyhedral_Homotopies;

  procedure Black_Box_Polyhedral_Homotopies
               ( file : in file_type; nt : in natural32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    mv,smv,tmv : natural32;
    r : integer32;
    mix,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    lq : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    timer : Timing_Widget;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    wsols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in black_box_mixed_volumes.");
      put_line("Black_Box_Polyhedral_Homotopies 6 ...");
    end if;
    stlb := Floating_Lifting_Functions.Lifting_Bound(p);
    tstart(timer);
    Pipelined_Polyhedral_Drivers.Pipelined_Polyhedral_Homotopies
      (integer32(nt),true,stlb,p,r,mix,perm,lifsup,mcc,tmv,lq,q,qsols);
    Set_Continuation_Parameter(qsols,zero);
    Apply_Induced_Permutation(p,stlb,r,perm,mix,mcc);
    Black_Box_Root_Counters.Pipelined_Stable_Continuation
        (true,r,mix,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
    Set_Continuation_Parameter(qsols0,zero);
    tstop(timer);
    Root_Counters_Output.Write_Mixed_Volumes(file,mv,smv);
    new_line(file);
    put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
    put_line(file,q);
    new_line(file);
    put_line(file,"START SOLUTIONS :");
    if Is_Null(qsols0) then
      put(file,Length_Of(qsols),natural32(q'last),qsols);
    else
      Push(qsols,wsols); Push(qsols0,wsols);
      put(file,Length_Of(wsols),natural32(q'last),wsols);
    end if;
    new_line(file);
    print_times(file,timer,"Pipelined Polyhedral Continuation");
    flush(file);
    rc := smv;
    rocotime := 0.0;
    hocotime := Elapsed_User_Time(timer);
  end Black_Box_Polyhedral_Homotopies;

end Black_Box_Mixed_Volumes;
