with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Numbers_io;                         use Numbers_io;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Homotopy;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Continuation_Parameters;            use Continuation_Parameters;
with Standard_Continuation_Data;         use Standard_Continuation_Data;
with Standard_Path_Trackers;             use Standard_Path_Trackers;
with Drivers_for_Homotopy_Creation;      use Drivers_for_Homotopy_Creation;
with Main_Poly_Continuation;             use Main_Poly_Continuation;

procedure Driver_for_Winding_Numbers
             ( file : in file_type; p : in Poly_Sys;
               sols : in out Solution_List ) is

  procedure Write_Statistics_and_Condition
                ( file : in file_type; i : in integer32;
                  nstep,nfail,niter,nsyst : in natural32;
                  rcond : in double_float ) is

  -- DESCRIPTION :
  --   Writes the computing statistics of the ith path on file, followed
  --   by the estimate for the inverse of the condition number.

  begin
    put(file,"========================================");
    put_line(file,"===================================");
    put(file,"== "); put(file,i,1); put(file," = ");
    put(file," #step : "); put(file,nstep,3);
    put(file," #fail : "); put(file,nfail,2);
    put(file," #iter : "); put(file,niter,3);
    if nsyst /= niter
     then put(file," #syst : "); put(file,nsyst,3);
    end if;
    put(file," = ");
    put(file," rco : "); put(file,rcond,2,3,3);
    put_line(file," ==");
  end Write_Statistics_and_Condition;

  procedure Winding_Numbers
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean; target : in Complex_Number;
                 mwc : in natural32 ) is

  -- DESCRIPTION :
  --   Computes the winding number for the given list of solutions.

  -- REQUIRED :
  --   The homotopy is already defined and stored in the package homotopy.

    sa : Solu_Info_Array(1..integer32(Length_Of(sols))) := Deep_Create(sols);
    pen : constant Pred_Pars := Continuation_Parameters.Create_End_Game;
    cen : constant Corr_Pars := Continuation_Parameters.Create_End_Game;
    tol : constant double_float := 1.0E-10;
    epslop : constant double_float := 1.0E-6;
    wc : natural32;
    sum,allsum : Standard_Complex_Vectors.Vector(sa(sa'first).sol.v'range);

    procedure CCont2 is
      new Circular_Single_Conditioned_Reporting_Continue
            (Max_Norm,Standard_Homotopy.Eval,
             Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    for i in sa'range loop
      if AbsVal(Create(1.0) - sa(i).sol.t) > epslop then
        CCont2(file,sa(i),target,tol,epslop,wc,mwc,sum,allsum,false,pen,cen);
        sa(i).sol.m := integer32(wc);
        sa(i).sol.v := allsum;
        sa(i).sol.t := target;
        Write_Statistics_and_Condition
          (file,i,sa(i).nstep,sa(i).nfail,sa(i).niter,sa(i).nsyst,sa(i).rcond);
        put(file,sa(i).sol.all);
      end if;
    end loop;
    Clear(sols);
    sols := Shallow_Create(sa);
  end Winding_Numbers;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the parameters,
  --   constructs the homotopy and then tracks the paths.

    ls : String_Splitters.Link_to_Array_of_Strings := null;
    infile : file_type;
    timer : timing_widget;
    pp,q : Poly_Sys(p'range);
    proj : boolean := false;
    target : Complex_Number;
    ans : character;
    oc,deci,max_wc : natural32 := 0;

  begin
   -- READING MAXIMUM WINDING NUMBER :
    new_line;
    put("Give the maximum winding number : "); Read_Natural(max_wc);
   -- READING THE START SYSTEM :
    new_line;
    put_line("Reading the name of the file for start system.");
    Read_Name_and_Open_File(infile); get(infile,q);
    new_line;
    put_line(file,"THE START SYSTEM :");
    new_line(file); put(file,q); new_line(file);
   -- CONSTRUCTING THE HOMOTOPY AND TUNING OF PARAMETERS :
    Copy(p,pp);
    Driver_for_Homotopy_Construction(file,ls,pp,q,sols,target,deci);
    proj := (Number_of_Unknowns(q(q'first)) > natural32(q'last));
    Driver_for_Continuation_Parameters(file);
    new_line;
    put("Do you want intermediate output during continuation ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Driver_for_Process_io(file,oc);
    end if;
   -- COMPUTATION OF WINDING NUMBERS :
    new_line;
    put_line("No more input desired.  Computing winding numbers ...");
    put_line("The program can now run in the background.");
    new_line;
    tstart(timer);
    Winding_Numbers(file,sols,proj,target,max_wc);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"computation of winding numbers");
    new_line(file);
  end Main;

begin
  Main;
end Driver_for_Winding_Numbers;
