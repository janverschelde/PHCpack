with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Time_Stamps;                        use Time_Stamps;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Write_Seed_Number;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Symbol_Table;
with Standard_Complex_Laur_Strings;
with Standard_Laur_Poly_Convertors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Floating_Numbers;
with Multprec_Complex_Solutions;
with Scaling_Methods;
with Main_Reduction;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Main_Root_Counters;
with Drivers_for_Homotopy_Creation;      use Drivers_for_Homotopy_Creation;
with Driver_for_Root_Refining;
with String_System_Readers;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with Greeting_Banners;

package body Polynomial_Homotopy_Continuation is

  procedure Display_all_Options ( nt : in natural32 ) is

    o : array(0..26) of string(1..65);

  begin
    put("Running full mode,");
    if nt = 0
     then put(" no tasking.");
     else put(" with "); put(nt,1); put(" tasks.");
    end if;
    put_line("  Note also the following options:");
    o(0) := "  phc -0 : random numbers with fixed seed for repeatable runs    ";
    o(1) := "  phc -a : solving polynomial systems equation-by-equation       ";
    o(2) := "  phc -b : batch or black-box processing, the blackbox solver    ";
    o(3) := "  phc -B : blackbox numerical irreducible decomposition solver   ";
    o(4) := "  phc -c : irreducible decomposition for solution components     ";
    o(5) := "  phc -d : linear and nonlinear reduction w.r.t. the total degree";
    o(6) := "  phc -e : SAGBI/Pieri/Littlewood-Richardson homotopies          ";
    o(7) := "  phc -f : factor pure dimensional solution set into irreducibles";
    o(8) := "  phc -g : checking whether an input system has the right syntax ";
    o(9) := "  phc -h : displays help, e.g.: phc -h -b or phc -b -h, or --help";
    o(10):= "  phc -j : path tracking with algorithmic differentiation        ";
    o(11):= "  phc -k : realization of dynamic output feedback placing poles  ";
    o(12):= "  phc -l : witness set for hypersurface cutting with random line ";
    o(13):= "  phc -m : mixed volumes via lift+prune, MixedVol, and DEMiCs    ";
    o(14):= "  phc -o : write order of symbols after parsing polynomial system";
    o(15):= "  phc -p : polynomial continuation by a homotopy in one parameter";
    o(16):= "  phc -q : tracking solution paths with incremental read/write   ";
    o(17):= "  phc -r : root counting and construction of start systems       ";
    o(18):= "  phc -s : equation and variable scaling on system and solutions ";
    o(19):= "  phc -t : tasking for tracking paths using multiple threads     ";
    o(20):= "  phc -u : power series, Pade approximants for path tracking     ";
    o(21):= "  phc -v : verification, refinement and purification of solutions";
    o(22):= "  phc -V : in verbose mode, at some given verbose level          ";
    o(23):= "  phc -w : witness set intersection using diagonal homotopies    ";
    o(24):= "  phc -x : convert solutions from PHCpack into Python dictionary ";
    o(25):= "  phc -y : sample points from an algebraic set, given witness set";
    o(26):= "  phc -z : strip phc output solution lists into Maple format     ";
    for i in o'range loop
      put_line(o(i));
    end loop;
    put_line("Options may be combined, e.g.: phc -b -0 or phc -0 -b.");
    put_line("To run the blackbox solver with 8 threads, do phc -b -t8.");
    put_line("Use -b2 or -b4 for double double or quad double precision.");
  end Display_all_Options;

  procedure Main_Polynomial_Solver 
              ( file : in file_type; nt : in natural32; p : in Poly_Sys;
                ls : in String_Splitters.Link_to_Array_of_Strings;
                vrb : in integer32 := 0 ) is

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
    if vrb > 0 then
      put_line("-> in polynomial_homotopy_continuation.");
      put_line("Main_Polynomial_Solver ...");
    end if;
    Copy(p,scalp);
    Scaling_Methods.Main(file,scalp,basis,scalvec,vrb-1);
    Main_Reduction.Reduce(file,scalp,roco,true,vrb-1);
    Copy(scalp,projp);
    Main_Root_Counters.Polynomial_Main(file,nt,projp,q,true,sols,roco,vrb-1);
    if Length_Of(sols) > 0 then
      Driver_for_Homotopy_Construction(file,ls,projp,q,sols,target,deci,vrb-1);
      proj := (Number_of_Unknowns(p(p'first)) > natural32(p'last));
      if Head_Of(sols).t /= Create(0.0)
       then Set_Continuation_Parameter(sols,Create(0.0));
      end if;
      if deci <= 16 then
        Driver_for_Standard_Continuation
          (file,sols,proj,target=>target,verbose=>vrb-1);
        Driver_for_Root_Refining(file,scalp,p,basis,scalvec,sols,vrb-1);
      elsif deci <= 32 then
        ddsols := DoblDobl_Complex_Solutions.Create(sols);
        Driver_for_DoblDobl_Continuation
          (file,ddsols,target=>target,verbose=>vrb-1);
      elsif deci <= 64 then
        qdsols := QuadDobl_Complex_Solutions.Create(sols);
        Driver_for_QuadDobl_Continuation
          (file,qdsols,target=>target,verbose=>vrb-1);
      else
        mpsols := Multprec_Complex_Solutions.Create(sols);
        size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
        Multprec_Complex_Solutions.Set_Size(mpsols,size);
        Driver_for_Multprec_Continuation
          (file,mpsols,proj,deci,target,verbose=>vrb-1);
      end if;
    end if;
  end Main_Polynomial_Solver;

  procedure Main_Laurent_Solver 
              ( file : in file_type; nt : in natural32;  p : in Laur_Sys;
                vrb : in integer32 := 0 ) is

    cp,q : Laur_Sys(p'range);
    sols : Solution_List;
    roco : natural32;
    poco : duration := 0.0;

  begin
    if vrb > 0 then
      put("-> in polynomial_homotopy_continuation.");
      put_line("Main_Laurent_Solver ...");
    end if;
    Copy(p,cp);
    Main_Root_Counters.Laurent_Main(file,nt,cp,q,sols,roco,vrb-1);
    if Length_Of(sols) > 0 then
      if Head_Of(sols).t /= Create(0.0)
       then Set_Continuation_Parameter(sols,Create(0.0));
      end if;
      if nt = 0
       then Black_Box_Polynomial_Continuation(file,cp,q,sols,poco,vrb-1);
       else Black_Box_Polynomial_Continuation
              (file,integer32(nt),cp,q,sols,poco,vrb-1);
      end if;
    end if;
  end Main_Laurent_Solver;

  procedure Start_Main
              ( start_moment : in Ada.Calendar.Time; nt : in natural32;
                outfilename : in string; q : in Laur_Sys;
                ls : in String_Splitters.Link_to_Array_of_Strings;
                vrb : in integer32 := 0 ) is

    use Ada.Calendar;

    ended_moment : Time;
    timer : Timing_Widget;
    outpt : file_type;

  begin
    if vrb > 0
     then put_line("-> in polynomial_homotopy_continuation.Start_Main ...");
    end if;
    Create_Output_File(outpt,outfilename);
    put(outpt,q'last,1);
    new_line(outpt);
    for k in ls'range loop
      put_line(outpt,ls(k).all);
    end loop;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(q) then
      tstart(timer);
      Main_Laurent_Solver(outpt,nt,q,vrb-1);
      tstop(timer);
    else
      declare
        use Standard_Laur_Poly_Convertors;
        p : constant Poly_Sys(q'range)
          := Positive_Laurent_Polynomial_System(q);
      begin
        tstart(timer);
        Main_Polynomial_Solver(outpt,nt,p,ls,vrb-1);
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
  end Start_Main;

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

    inft : file_type;
    ls : String_Splitters.Link_to_Array_of_Strings;
    n,m : natural32;

    use String_Splitters;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in polynomial_homotopy_continuation.Main ...");
    end if;
    new_line;
    Display_all_Options(nt);
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
      Start_Main(start_moment,nt,outfilename,q,ls,verbose-1);
    end;
  end Main;

end Polynomial_Homotopy_Continuation;
