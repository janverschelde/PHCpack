with Ada.Calendar;
with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Arrays_of_Floating_Vector_Lists;  
with Floating_Mixed_Subdivisions;
with Black_Mixed_Volume_Computations;    use Black_Mixed_Volume_Computations;
with Black_Polyhedral_Continuations;     use Black_Polyhedral_Continuations;

procedure babldmvc ( nt : in natural32; infilename,outfilename : in string ) is

  start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
  ended_moment : Ada.Calendar.Time;

  procedure Read_System
              ( file : in out file_type; filename : in string;
                lp : out Link_to_Poly_Sys ) is

  -- DESCRIPTION :
  --   Attempts to read a polynomial system from file.

  begin
    if filename /= "" then
      Open_Input_File(file,filename);
      get(file,lp);
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   lp := null; return;
  end Read_System;

  procedure Call_MixedVol
      ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
        mivo,stmv : out natural32;
        stlb : out double_float;
        lifsup : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
        mix,perm,iprm : out Link_to_Vector;
        orgmcc,stbmcc : out Floating_Mixed_Subdivisions.Mixed_Subdivision ) is

    mixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    mv,smv,tmv,orgcnt,stbcnt : natural32;

  begin
    Black_Box_Mixed_Volume_Computation
      (p,mix,perm,iprm,stlb,lifsup,mixsub,orgmcc,stbmcc,
       mv,smv,tmv,orgcnt,stbcnt);
    mivo := mv; stmv := smv;
  end Call_MixedVol;

  procedure Polyhedral_Solver ( file : in file_type; p : in out Poly_Sys ) is

  -- DESCRIPTION :
  --   Calls the MixedVol algorithm and solves a random coefficient system.

    timer : Timing_Widget;
    q : Poly_Sys(p'range);
    qsols,qsols0 : Solution_List;
    mix,perm,iprm : Link_to_Vector;
    lifsup : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
    stlb : double_float;
    orgmcc,stbmcc : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    mv,stmv : natural32;

  begin
    tstart(timer);
    Call_MixedVol(p,mv,stmv,stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc);
    tstop(timer);
    new_line(file);
    put(file,"mixed volume : "); put(file,mv,1); new_line(file);
    put(file,"stable mixed volume : ");
    put(file,stmv,1); new_line(file);
    new_line(file);
    print_times(file,timer,"Mixed-Volume Computation");
    if mv > 0 then
      tstart(timer);
      Black_Box_Polyhedral_Continuation
        (integer32(nt),p,mix,stlb,lifsup.all,orgmcc,stbmcc,q,qsols,qsols0);
      tstop(timer);
      new_line(file);
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      new_line(file);
      put_line(file,q);
      new_line(file);
      put_line(file,"SOLUTIONS with nonzero coordinates:");
      new_line(file);
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      new_line(file);
      if not Is_Null(qsols0) then
        put_line(file,"SOLUTIONS with zero coordinates :");
        new_line(file);
        put(file,Length_Of(qsols0),natural32(Head_Of(qsols0).n),qsols0);
      end if;
      print_times(file,timer,"Polyhedral Continuation");
    end if;
  end Polyhedral_Solver;

  procedure Main is

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;

  begin
    Read_System(infile,infilename,lp);
    if lp = null then
      new_line;
      get(lp);
    end if;
    Close(infile);
    Create_Output_File(outfile,outfilename);
    put(outfile,natural32(lp'last),lp.all);
    Polyhedral_Solver(outfile,lp.all);
    ended_moment := Ada.Calendar.Clock;
    new_line(outfile);
    put(outfile,"PHC ran from ");
    Write_Time_Stamp(outfile,start_moment);
    put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
    put_line(outfile,".");
    Write_Elapsed_Time(outfile,start_moment,ended_moment);
    Close(outfile);
  end Main;

begin
  Main;
end babldmvc;
