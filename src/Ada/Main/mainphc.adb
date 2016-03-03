with Ada.Calendar;                       use Ada.Calendar;
with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Time_Stamps;                        use Time_Stamps;
with String_Splitters;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Write_Seed_Number;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Symbol_Table;
with Standard_Complex_Laur_Strings;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Laur_Poly_Convertors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Floating_Numbers;
with Multprec_Complex_Solutions;
with Drivers_for_Scaling;                use Drivers_for_Scaling;
with Drivers_for_Reduction;              use Drivers_for_Reduction;
with Drivers_for_Poly_Continuation;      use Drivers_for_Poly_Continuation;
with Drivers_for_Root_Counts;            use Drivers_for_Root_Counts;
with Drivers_for_Homotopy_Creation;      use Drivers_for_Homotopy_Creation;
with Driver_for_Root_Refining;
with String_System_Readers;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with Greeting_Banners;
--with Bye_Bye_Message;

procedure mainphc ( nt : in natural32; infilename,outfilename : in string ) is

  start_moment : constant Time := Clock;
  ended_moment : Time;

  procedure Display_Options is

  -- DESCRIPTION :
  --   Displays an overview of all options on screen.

    o : array(0..21) of string(1..65);

  begin
    put("Running full mode,");
    if nt = 0
     then put(" no tasking.");
     else put(" with "); put(nt,1); put(" tasks.");
    end if;
    put_line("  Note also the following options:");
    o(0) := "  phc -0 : random numbers with fixed seed for repeatable runs    ";
    o(1) := "  phc -a : Solving polynomial systems equation-by-equation       ";
    o(2) := "  phc -b : Batch or black-box processing                         ";
    o(3) := "  phc -c : Irreducible decomposition for solution components     ";
    o(4) := "  phc -d : Linear and nonlinear Reduction w.r.t. the total degree";
    o(5) := "  phc -e : SAGBI/Pieri homotopies to intersect linear subspaces  ";
    o(6) := "  phc -f : Factor pure dimensional solution set into irreducibles";
    o(7) := "  phc -g : Checking whether an input system has the right syntax ";
    o(8) := "  phc -k : realization of dynamic output feedback placing poles  ";
    o(9) := "  phc -l : Witness Set for Hypersurface cutting with Random Line ";
    o(10):= "  phc -m : Mixed-Volume Computation via lift+prune and MixedVol  ";
    o(11):= "  phc -o : write order of symbols after parsing polynomial system";
    o(12):= "  phc -p : Polynomial Continuation by a homotopy in one parameter";
    o(13):= "  phc -q : Tracking Solution Paths with incremental read/write   ";
    o(14):= "  phc -r : Root counting and Construction of start systems       ";
    o(15):= "  phc -s : Equation and variable Scaling on system and solutions ";
    o(16):= "  phc -t : Tasking for tracking paths using multiple threads     ";
    o(17):= "  phc -v : Verification, refinement and purification of solutions";
    o(18):= "  phc -w : Witness Set Intersection using Diagonal Homotopies    ";
    o(19):= "  phc -x : convert solutions from PHCpack into Python dictionary ";
    o(20):= "  phc -y : sample points from an algebraic set, given witness set";
    o(21):= "  phc -z : strip phc output solution lists into Maple format     ";
    for i in o'range loop
      put_line(o(i));
    end loop;
    put_line("Options may be combined, e.g.: phc -b -0 or phc -0 -b.");
    put_line("To run the blackbox solver with 8 threads, do phc -b -t8.");
    put_line("Use -b2 or -b4 for double double or quad double precision.");
  end Display_Options;

  procedure Main_Polynomial_Solver 
              ( file : in file_type; p : in Poly_Sys;
                ls : in String_Splitters.Link_to_Array_of_Strings ) is

  -- DESCRIPTION :
  --   This is the main interactive solver for polynomial systems,
  --   running through all drivers.

    q,scalp,projp : Poly_Sys(p'range);
    target : Complex_Number;
    basis,roco,deci,size : natural32 := 0;
    scalvec : Link_to_Vector;
    sols : Solution_List;
    ddsols : DoblDobl_Complex_Solutions.Solution_List;
    qdsols : QuadDobl_Complex_Solutions.Solution_List;
    mpsols : Multprec_Complex_Solutions.Solution_List;
    proj : boolean;

  begin
    Copy(p,scalp);
    Driver_for_Scaling(file,scalp,basis,scalvec);
    Driver_for_Reduction(file,scalp,roco,true);
    Copy(scalp,projp);
    Driver_for_Root_Counts(file,projp,q,true,sols,roco);
    if Length_Of(sols) > 0 then
      Driver_for_Homotopy_Construction(file,ls,projp,q,sols,target,deci);
      proj := (Number_of_Unknowns(p(p'first)) > natural32(p'last));
      if Head_Of(sols).t /= Create(0.0)
       then Set_Continuation_Parameter(sols,Create(0.0));
      end if;
      if deci <= 16 then
        Driver_for_Standard_Continuation(file,sols,proj,target=>target);
        Driver_for_Root_Refining(file,scalp,p,basis,scalvec,sols);
      elsif deci <= 32 then
        ddsols := DoblDobl_Complex_Solutions.Create(sols);
        Driver_for_DoblDobl_Continuation(file,ddsols,target=>target);
      elsif deci <= 64 then
        qdsols := QuadDobl_Complex_Solutions.Create(sols);
        Driver_for_QuadDobl_Continuation(file,qdsols,target=>target);
      else
        mpsols := Multprec_Complex_Solutions.Create(sols);
        size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
        Multprec_Complex_Solutions.Set_Size(mpsols,size);
        Driver_for_Multprec_Continuation(file,mpsols,proj,deci,target);
      end if;
    end if;
  end Main_Polynomial_Solver;

  procedure Main_Laurent_Solver 
              ( file : in file_type; p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   This is the main interactive solver for Laurent systems,
  --   primarily guiding through polyhedral homotopies.

    cp,q : Laur_Sys(p'range);
    sols : Solution_List;
    roco : natural32;
    poco : duration := 0.0;

  begin
    Copy(p,cp);
    Driver_for_Root_Counts(file,cp,q,sols,roco);
    if Length_Of(sols) > 0 then
      if Head_Of(sols).t /= Create(0.0)
       then Set_Continuation_Parameter(sols,Create(0.0));
      end if;
      if nt = 0
       then Black_Box_Polynomial_Continuation(file,cp,q,sols,poco);
       else Black_Box_Polynomial_Continuation
              (file,integer32(nt),cp,q,sols,poco);
      end if;
    end if;
  end Main_Laurent_Solver;

  procedure Main_Dispatch 
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                ls : in String_Splitters.Link_to_Array_of_Strings ) is

   -- use Standard_Complex_Laur_Systems;

    timer : Timing_Widget;
    outpt : file_type;

  begin
    Create_Output_File(outpt,outfilename);
    put(outpt,q'last,1);
    new_line(outpt);
    for k in ls'range loop
      put_line(outpt,ls(k).all);
    end loop;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(q) then
      tstart(timer);
      Main_Laurent_Solver(outpt,q);
      tstop(timer);
    else
      declare
        use Standard_Laur_Poly_Convertors;
        p : constant Poly_Sys(q'range)
          := Positive_Laurent_Polynomial_System(q);
      begin
        tstart(timer);
        Main_Polynomial_Solver(outpt,p,ls);
        tstop(timer);
      end;
    end if;
    new_line(outpt);
    print_times(outpt,timer,"solving the polynomial system");
    new_line(outpt);
   -- put(outpt,Bye_Bye_Message);
    ended_moment := Clock;
    put(outpt,"PHC ran from "); Write_Time_Stamp(outpt,start_moment);
    put(outpt," till ");        Write_Time_Stamp(outpt,ended_moment);
    put_line(outpt,".");
    Write_Elapsed_Time(outpt,start_moment,ended_moment);
    Write_Seed_Number(outpt);
    put_line(outpt,Greeting_Banners.Version);
    Close(outpt);
  end Main_Dispatch;

  procedure String_Main is

  -- DESCRIPTION :
  --   Reads in the system as a pointer to an array of strings
  --   and converts to a Laurent polynomial system with standard
  --   complex coefficients for a first dispatch.

    inft : file_type;
    ls : String_Splitters.Link_to_Array_of_Strings;
    n,m : natural32;

    use String_Splitters;

  begin
    new_line;
    Display_Options;
    String_System_Readers.Read_System(inft,infilename,n,m,ls);
    if ls = null then
      new_line;
      put_line("Reading the target polynomial system...");
      Read_Name_and_Open_File(inft);
      String_Splitters.get(inft,natural(n),natural(m),ls);
      close(inft);
    end if;
    Symbol_Table.Init(m);
    declare
      q : constant Standard_Complex_Laur_Systems.Laur_Sys(1..integer32(n))
         := Standard_Complex_Laur_Strings.Parse(m,ls.all);
    begin
      Main_Dispatch(q,ls);
    end;
  end String_Main;

begin
  String_Main;
end mainphc;
