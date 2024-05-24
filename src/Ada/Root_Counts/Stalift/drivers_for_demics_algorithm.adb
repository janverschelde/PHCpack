with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_Poly_Laur_Convertors;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with Main_Poly_Continuation;
with Floating_Lifting_Functions;
with Floating_Mixed_Subdivisions_io;
with Drivers_for_Static_Lifting;
with Drivers_for_MixedVol_Algorithm;
with DEMiCs_Output_Data;
with DEMiCs_Algorithm;                   use DEMiCs_Algorithm;
with Pipelined_Polyhedral_Homotopies;    use Pipelined_Polyhedral_Homotopies;

with Standard_Floating_Numbers_io;
 use Standard_Floating_Numbers_io;

package body Drivers_for_DEMiCs_Algorithm is

  procedure DEMiCs_Algorithm_Info is

    i : array(1..10) of string(1..65);

  begin
    i(1) :="  The DEMiCs Algorithm calls the code of Tomohiko Mizutani,      ";
    i(2) :="Akiko Takeda, and Masakazu Kojima.  The algorithm is described in";
    i(3) :="Discrete Comput. Geom. 37(3):351-367, 2007.  The software DEMiCs ";
    i(4) :="is published in Software for Algebraic Geometry, Springer, 2008. ";
    i(5) :="DEMiCs stands for Dynamic Enumeration of Mixed Cells and applies ";
    i(6) :="a greedy strategy to run through the tree of face combinations   ";
    i(7) :="which span all mixed cells.  For many different Newton polytopes ";
    i(8) :="DEMiCs is faster than MixedVol, producing cells at a faster pace.";
    i(9) :="Compared to other lift-and-prune strategies, only random lifting ";
   i(10) :="values on the supports are supported.                            ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end DEMiCs_Algorithm_Info;

  procedure BlackBox_DEMiCs_Algorithm
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
              lif : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                vrblvl : in integer32 := 0 ) is  

    dim : constant integer32 := p'last;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_demics_algorithm.");
      put_line("Blackbox_DEMiCs_Algorithm 1 ...");
    end if;
    Extract_Supports(p,mix,sup,false); -- verbose is false
    Call_DEMiCs(mix,sup,false);
    declare
      lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    begin
      Process_Output(dim,mix,sup,lifsup,mcc,false);
      lif := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(lifsup);
    end;
    mv := natural32(DEMiCs_Output_Data.mixed_volume);
  end BlackBox_DEMiCs_Algorithm;

  procedure BlackBox_DEMiCs_Algorithm
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
              lif : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                vrblvl : in integer32 := 0 ) is  

    dim : constant integer32 := p'last;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_demics_algorithm.");
      put_line("Blackbox_DEMiCs_Algorithm 2 ...");
    end if;
    Extract_Supports(p,mix,sup,false);  -- verbose is false
    Call_DEMiCs(mix,sup,false);
    declare
      lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    begin
      Process_Output(dim,mix,sup,lifsup,mcc,false);
      lif := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(lifsup);
    end;
    mv := natural32(DEMiCs_Output_Data.mixed_volume);
  end BlackBox_DEMiCs_Algorithm;

  procedure BlackBox_DEMiCs_Algorithm
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
              lif : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv,smv,tmv : out natural32;
                vrblvl : in integer32 := 0 ) is  

    dim : constant integer32 := p'last;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    stlb : constant double_float
         := Floating_Lifting_Functions.Lifting_Bound(p); -- ,1.0e+10);

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_demics_algorithm.");
      put_line("Blackbox_DEMiCs_Algorithm 3 ...");
    end if;
    Extract_Supports(p,mix,sup,false);    -- verbose is false
    Call_DEMiCs(mix,sup,true,stlb,false); -- stable is true
    declare
      lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    begin
      Process_Output(dim,mix,sup,lifsup,mcc,false);
      lif := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(lifsup);
    end;
    Drivers_for_Static_Lifting.Floating_Volume_Computation
      (dim,stlb,mix.all,mcc,mv,smv,tmv);
  end BlackBox_DEMiCs_Algorithm;

  procedure Write_Random_Coefficient_System
              ( file : in file_type; ranfile : in out file_type;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : in Standard_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in drivers_for_demics_algorithm.");
      put_line("Write_Random_Coefficient_System 1 ...");
    end if;
    new_line(file);
    put_line(file,"RANDOM COEFFICIENT SYSTEM :");
    Standard_System_and_Solutions_io.put_line(file,q,qsols);
    Standard_System_and_Solutions_io.put_line(ranfile,q,qsols);
    close(ranfile);
  end Write_Random_Coefficient_System;

  procedure Write_Random_Coefficient_System
              ( file : in file_type; ranfile : in out file_type;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                qsols,qsols0 : in Standard_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in drivers_for_demics_algorithm.");
      put_line("Write_Random_Coefficient_System 2 ...");
    end if;
    new_line(file);
    put_line(file,"RANDOM COEFFICIENT SYSTEM :");
    Standard_System_and_Solutions_io.put_line(file,q,qsols);
    Standard_System_and_Solutions_io.put_line(ranfile,q,qsols);
    if Standard_Complex_Solutions.Is_Null(qsols0) then
      new_line(file);
      put_line(file,"No solutions with zero coordinate values found.");
    else
      new_line(file);
      put_line(file,"THE SOLUTIONS WITH ZERO COORDINATE VALUES : ");
      put(file,Standard_Complex_Solutions.Length_Of(qsols0),
          natural32(Standard_Complex_Solutions.Head_Of(qsols0).n),qsols0);
      new_line(ranfile);
      put_line(ranfile,"THE SOLUTIONS WITH ZERO COORDINATE VALUES : ");
      put(ranfile,Standard_Complex_Solutions.Length_Of(qsols0),
          natural32(Standard_Complex_Solutions.Head_Of(qsols0).n),qsols0);
    end if;
    close(ranfile);
  end Write_Random_Coefficient_System;

  procedure Run_Polyhedral_Homotopies
              ( file : in file_type;
                mcc2file,ranstart,contrep : in boolean;
                subfile,ranfile : in out file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in integer32;
                mix,perm : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                orgsup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stable : in boolean; stlb : in double_float;
                timer : in Timing_Widget;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                qsols0 : out Standard_Complex_Solutions.Solution_List;
                mv,smv,tmv : out natural32;
                vrblvl : in integer32 := 0 ) is

    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    mcc,orgmcc,stbmcc : Mixed_Subdivision;
    orgcnt,stbcnt : natural32;
    qtimer : Timing_Widget;

    use Drivers_for_MixedVol_Algorithm;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_demics_algorithm");
      put_line("Run_Polyhedral_Homotopies ...");
    end if;
    Process_Output(dim,mix,sup,lifsup,mcc,false);
    new_line(file);
    put_line(file,"The lifted supports :");
    Floating_Mixed_Subdivisions_io.put(file,lifsup);
    if mcc2file then
      Floating_Mixed_Subdivisions_io.put(subfile,natural32(dim),mix.all,mcc);
      close(subfile);
    end if;
    if not stable then
      put_line(file,"A mixed-cell configuration :");
      Floating_Mixed_Subdivisions_io.put(file,natural32(dim),mix.all,mcc,mv);
      put(file,"The mixed volume : "); put(file,mv,1); new_line(file);
    else
      Drivers_for_Static_Lifting.Floating_Volume_Computation
        (file,dim,stlb,mix.all,mcc,mv,smv,tmv);
    end if;
    new_line(file);
    print_times(file,timer,"DEMiCs Algorithm");
    if ranstart then
      if not stable then
        tstart(qtimer);
        Random_Coefficient_System(0,dim,mix.all,lifsup,mcc,q,qsols);
        tstop(qtimer);
        Write_Random_Coefficient_System(file,ranfile,q,qsols);
      else
        Split_Original_Cells(mcc,stlb,orgmcc,stbmcc,orgcnt,stbcnt);
        tstart(qtimer);
        Polyhedral_Continuation
          (file,0,stable,contrep,dim,mix'last,stlb,mix,perm,
           p,orgsup,mcc,orgmcc,stbmcc,q,qsols,qsols0);
        tstop(qtimer);
        Write_Random_Coefficient_System(file,ranfile,q,qsols,qsols0);
      end if;
      new_line(file);
      print_times(file,qtimer,"polyhedral homotopies");
    end if;
  end Run_Polyhedral_Homotopies;

  procedure Run_DEMiCs_Algorithm
              ( file : in file_type; nt : in integer32;
                mcc2file,ranstart : in boolean;
                subfile,ranfile : in out file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in integer32;
                mix,perm : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stable : in boolean; stlb : in double_float;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                qsols0 : out Standard_Complex_Solutions.Solution_List;
                mv,smv,tmv : out natural32;
                vrblvl : in integer32 := 0 ) is

    timer : Timing_Widget;
    lif : Standard_Floating_VecVecs.Link_to_VecVec;
    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    orgsup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(sup'range);
    mcc : Mixed_Subdivision;
    oc : natural32;
    contrep : boolean;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_demics_algorithm.");
      put_line("Run_DEMiCs_Algorithm ...");
    end if;
   -- start system based on supports without artificial origin: orgsup
    if stable
     then Arrays_of_Integer_Vector_Lists.Copy(sup,orgsup);
     else orgsup := sup;
    end if;
    if ranstart then -- allow for tuning of continuation parameters
      new_line;
      Main_Poly_Continuation.Driver_for_Continuation_Parameters(file);
      new_line;
      Main_Poly_Continuation.Driver_for_Process_io(file,oc);
      contrep := (oc /= 0);
    end if;
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt < 2 or not ranstart then -- no multitasking or no path tracking
      if stable then
        tstart(timer);
        Call_DEMiCs(mix,sup,stable,stlb,false);
        tstop(timer);
      else
        tstart(timer);
        Call_DEMiCs(mix,sup,false);
        tstop(timer);
      end if;
      Run_Polyhedral_Homotopies
        (file,mcc2file,ranstart,contrep,subfile,ranfile,p,dim,mix,perm,
         sup,orgsup,stable,stlb,timer,q,qsols,qsols0,mv,smv,tmv);
    else -- if random coefficient system is needed with pipelining
      tstart(timer);
      lif := Random_Lifting(mix,sup);
      Pipeline_Cells_to_Paths(dim,nt,mix,sup,lif,q,qsols,false);
      tstop(timer);
      Process_Output(dim,mix,sup,lifsup,mcc,false);
      new_line(file);
      put_line(file,"The lifted supports :");
      Floating_Mixed_Subdivisions_io.put(file,lifsup);
      put_line(file,"A mixed-cell configuration :");
      Floating_Mixed_Subdivisions_io.put(file,natural32(dim),mix.all,mcc,mv);
      put(file,"The mixed volume : "); put(file,mv,1); new_line(file);
      Write_Random_Coefficient_System(file,ranfile,q,qsols);
      new_line(file);
      print_times(file,timer,"pipelined polyhedral homotopies");
    end if;
    if stable
     then Arrays_of_Integer_Vector_Lists.Deep_Clear(orgsup);
    end if;
  end Run_DEMiCs_Algorithm;

  procedure Driver_for_DEMiCs_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                qsols0 : out Standard_Complex_Solutions.Solution_List;
                mv,smv,tmv : out natural32;
                vrblvl : in integer32 := 0 ) is

    lp,lq : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

    use Standard_Poly_Laur_Convertors;
    use Standard_Laur_Poly_Convertors;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_demics_algorithm.");
      put_line("Driver_for_DEMiCs_Algorithm 1 ...");
    end if;
    lp := Polynomial_to_Laurent_System(p);
    Driver_for_DEMiCs_Algorithm(file,nt,lp,lq,qsols,qsols0,mv,smv,tmv);
    q := Positive_Laurent_Polynomial_System(lq);
    Standard_Complex_Laur_Systems.Clear(lp);
    Standard_Complex_Laur_Systems.Clear(lq);
  end Driver_for_DEMiCs_Algorithm;

  procedure Driver_for_DEMiCs_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                qsols0 : out Standard_Complex_Solutions.Solution_List;
                mv,smv,tmv : out natural32;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := p'last;
    mix,perm : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    timer : Timing_Widget;
    mcc2file : boolean := false;
    ranstart : boolean := false;
    subfile,ranfile : file_type;
    ans : character;
    genuine : constant boolean
            := Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(p);
    stable : boolean;
    stlb : double_float;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_demics_algorithm.");
      put_line("Driver_for_DEMiCs_Algorithm 2 ...");
    end if;
    new_line;
    Drivers_for_Static_Lifting.Prompt_for_File(mcc2file,subfile);
    if genuine then
      stable := false; -- no stable mixed volume for genuine Laurent system
    else
      new_line;
      put("Do you want to compute the stable mixed volume ? (y/n) ");
      Ask_Yes_or_No(ans);
      stable := (ans = 'y');
      if stable then
        stlb := Floating_Lifting_Functions.Lifting_Bound(p); -- ,1.0e+10);
        if vrblvl > 0
         then put("The lifting bound :"); put(stlb); new_line;
        end if;
      else
        stlb := 0.0;
      end if;
    end if;
    new_line;
    put("Monitor the progress of the computations on screen ? (y/n) ");
    Ask_Yes_or_No(ans);
    DEMiCs_Output_Data.monitor := (ans = 'y');
    new_line;
    put("Run polyhedral homotopies ");
    put("to solve a random coefficient system ? (y/n) ");
    Ask_Yes_or_No(ans);
    ranstart := (ans = 'y');
    if ranstart then
      put("Reading the name of the file ");
      put_line("for the random coefficient system ...");
      Read_Name_and_Create_File(ranfile);
    end if;
    Extract_Supports(p,mix,perm,sup,false); -- verbose is false
    Run_DEMiCs_Algorithm
      (file,nt,mcc2file,ranstart,subfile,ranfile,
       p,dim,mix,perm,sup,stable,stlb,q,qsols,qsols0,mv,smv,tmv);
  end Driver_for_DEMiCs_Algorithm;

  procedure Driver_for_DEMiCs_Algorithm
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                vrblvl : in integer32 := 0 ) is

    q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols,qsols0 : Standard_Complex_Solutions.Solution_List;
    mv,smv,tmv : natural32;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_demics_algorithm.");
      put_line("Driver_for_DEMiCs_Algorithm 3 ...");
    end if;
    Driver_for_DEMiCs_Algorithm
      (file,nt,p,q,qsols,qsols0,mv,smv,tmv,vrblvl-1);
  end Driver_for_DEMiCs_Algorithm;

end Drivers_for_DEMiCs_Algorithm;
