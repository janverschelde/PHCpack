with Timing_Package;                     use Timing_Package;
with Black_Mixed_Volume_Computations;    use Black_Mixed_Volume_Computations;
with Black_Polyhedral_Continuations;     use Black_Polyhedral_Continuations;
with Root_Counters_Output;

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
    rocotime := Elapsed_User_Time(timer);
  end Mixed_Volume;

  procedure Construct_Start_System
              ( nt : in integer32;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                mix : in Link_to_Vector;
                stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                orgmcc,stbmcc : in Mixed_Subdivision;
                q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : in out Standard_Complex_Solutions.Solution_List;
                hocotime : out duration; verbose : in integer32 := 0 ) is

    timer : timing_widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_mixed_volumes.Construct_Start_System ...");
    end if;
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
  end Construct_Start_System;

  procedure Construct_Start_System
              ( file : in file_type; nt : in integer32;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                mix : in Link_to_Vector;
                stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                orgmcc,stbmcc : in Mixed_Subdivision;
                q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : in out Standard_Complex_Solutions.Solution_List;
                hocotime : out duration; verbose : in integer32 := 0 ) is

    timer : timing_widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_mixed_volumes.Construct_Start_System ...");
    end if;
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
  end Construct_Start_System;

  procedure Black_Box_Polyhedral_Homotopies
              ( nt : in integer32; silent : in boolean;
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
      put("-> in black_box_mixed_volumes");
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
      (nt,p,mix,stlb,lifsup,orgmcc,stbmcc,q,qsols,qsols0,hocotime,verbose-1);
  end Black_Box_Polyhedral_Homotopies;

  procedure Black_Box_Polyhedral_Homotopies
               ( nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
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
      (nt,p,mix,stlb,lifsup,orgmcc,stbmcc,q,qsols,qsols0,hocotime,verbose-1);
  end Black_Box_Polyhedral_Homotopies;

end Black_Box_Mixed_Volumes;
