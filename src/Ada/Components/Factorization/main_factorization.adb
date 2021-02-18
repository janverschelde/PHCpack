with Ada.Calendar;
with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Numbers_io;                         use Numbers_io;
with Multprec_Floating_Numbers;
with Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Main_Solution_Filters;
with Sampling_Machine;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Standard_Stacked_Sample_Grids;
with Multprec_Stacked_Sample_Grids;
with Make_Sample_Grids;                  use Make_Sample_Grids;
with Standard_Trace_Interpolators;
with Multprec_Trace_Interpolators;
with Homotopy_Membership_Tests;          use Homotopy_Membership_Tests;
with Multitasking_Membership_Tests;      use Multitasking_Membership_Tests;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Drivers_to_Factor_Components;       use Drivers_to_Factor_Components;
with Drivers_to_Factor_Polynomials;      use Drivers_to_Factor_Polynomials;
with Drivers_to_Breakup_Solutions;       use Drivers_to_Breakup_Solutions;
with Driver_for_Common_Factor;
with Greeting_Banners;
with Write_Seed_Number;

package body Main_Factorization is

  procedure Tune_Member_Tolerances ( restol,homtol : in out double_float ) is

    ans : character;

  begin
    loop
      new_line;
      put_line("TOLERANCES to decide whether a test point");
      put("  1. satisfies one polynomial : "); put(restol,3); new_line;
      put("  2. belongs to a witness set : "); put(homtol,3); new_line;
      put("Type 1 or 2 to change a tolerance, 0 to exit : ");
      Ask_Alternative(ans,"012");
      exit when (ans = '0');
      put("Give new tolerance for point to ");
      if ans = '1' then
        put("satisfy one polynomial : "); Read_Double_Float(restol);
      else
        put("belong to a witness set : "); Read_Double_Float(homtol);
      end if;
    end loop;
  end Tune_Member_Tolerances;

  procedure Standard_Multitasked_Membership_Test
              ( file : in file_type;
                nbtasks,dim : in natural32; homtol : in double_float;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                gpts,sols : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    start_moment,ended_moment : Ada.Calendar.Time;
    isols : Solution_List;
    timer : Timing_Widget;

  begin
    start_moment := Ada.Calendar.Clock;
    tstart(timer);
    isols := Standard_Membership_Filter(nbtasks,dim,homtol,ep,gpts,sols);
    tstop(timer);
    ended_moment := Ada.Calendar.Clock;
    put(file,"Tested ");
    put(file,Length_Of(sols),1);
    put(file," points, with ");
    put(file,nbtasks,1); put_line(file," tasks.");
    put(file,"Number of points not on component : ");
    put(file,Length_Of(isols),1); put_line(file,".");
    new_line(file);
    print_times(file,timer,"multitasked membership test");
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    new_line(file);
    put_line(file,"SOLUTIONS THAT FAILED MEMBERSHIP TEST :");
    put(file,Length_Of(isols),natural32(Head_Of(isols).n),isols);
  end Standard_Multitasked_Membership_Test;

  procedure DoblDobl_Multitasked_Membership_Test
              ( file : in file_type;
                nbtasks,dim : in natural32; homtol : in double_float;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                gpts,sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    start_moment,ended_moment : Ada.Calendar.Time;
    isols : Solution_List;
    timer : Timing_Widget;

  begin
    start_moment := Ada.Calendar.Clock;
    tstart(timer);
    isols := DoblDobl_Membership_Filter(nbtasks,dim,homtol,ep,gpts,sols);
    tstop(timer);
    ended_moment := Ada.Calendar.Clock;
    put(file,"Tested ");
    put(file,Length_Of(sols),1);
    put(file," points, with ");
    put(file,nbtasks,1); put_line(file," tasks.");
    put(file,"Number of points not on component : ");
    put(file,Length_Of(isols),1); put_line(file,".");
    new_line(file);
    print_times(file,timer,"multitasked membership test");
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    new_line(file);
    put_line(file,"SOLUTIONS THAT FAILED MEMBERSHIP TEST :");
    put(file,Length_Of(isols),natural32(Head_Of(isols).n),isols);
  end DoblDobl_Multitasked_Membership_Test;

  procedure QuadDobl_Multitasked_Membership_Test
              ( file : in file_type;
                nbtasks,dim : in natural32; homtol : in double_float;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                gpts,sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    start_moment,ended_moment : Ada.Calendar.Time;
    isols : Solution_List;
    timer : Timing_Widget;

  begin
    start_moment := Ada.Calendar.Clock;
    tstart(timer);
    isols := QuadDobl_Membership_Filter(nbtasks,dim,homtol,ep,gpts,sols);
    tstop(timer);
    ended_moment := Ada.Calendar.Clock;
    put(file,"Tested ");
    put(file,Length_Of(sols),1);
    put(file," points, with ");
    put(file,nbtasks,1); put_line(file," tasks.");
    put(file,"Number of points not on component : ");
    put(file,Length_Of(isols),1); put_line(file,".");
    new_line(file);
    print_times(file,timer,"multitasked membership test");
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    new_line(file);
    put_line(file,"SOLUTIONS THAT FAILED MEMBERSHIP TEST :");
    put(file,Length_Of(isols),natural32(Head_Of(isols).n),isols);
  end QuadDobl_Multitasked_Membership_Test;

  procedure Standard_Homotopy_Membership_Test
              ( nt : in natural32; vrblvl : in integer32 := 0 ) is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : double_float := 1.0E-10;
    homtol : double_float := 1.0E-6;
    timer : Timing_Widget;

  begin
    if vrblvl > 0 then
      put("-> in main_factorization.");
      put_line("Standard_Homotopy_Membership_Test ...");
    end if;
    Standard_Read_Embedding(lp,genpts,dim);
    new_line;
    put_line("The input format of the test points is the solutions format.");
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Tune_Member_Tolerances(restol,homtol);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
      tstop(timer);
      new_line(file);
      print_times(file,timer,"homotopy membership test");
    else
      Standard_Multitasked_Membership_Test
        (file,nt,dim,homtol,lp.all,genpts,sols);
    end if;
    new_line(file);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end Standard_Homotopy_Membership_Test;

  procedure DoblDobl_Homotopy_Membership_Test
              ( nt : in natural32; vrblvl : in integer32 := 0 ) is

    file : file_type;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : DoblDobl_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : double_float := 1.0E-10;
    homtol : double_float := 1.0E-6;
    timer : Timing_Widget;

  begin
    if vrblvl > 0 then
      put("-> in main_factorization.");
      put_line("DoblDobl_Homotopy_Membership_Test ...");
    end if;
    DoblDobl_Read_Embedding(lp,genpts,dim);
    new_line;
    put_line("The input format of the test points is the solutions format.");
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Tune_Member_Tolerances(restol,homtol);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
      tstop(timer);
      new_line(file);
      print_times(file,timer,"homotopy membership test");
    else
      DoblDobl_Multitasked_Membership_Test
        (file,nt,dim,homtol,lp.all,genpts,sols);
    end if;
    new_line(file);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end DoblDobl_Homotopy_Membership_Test;

  procedure QuadDobl_Homotopy_Membership_Test
              ( nt : in natural32; vrblvl : in integer32 := 0 ) is

    file : file_type;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts,sols : QuadDobl_Complex_Solutions.Solution_List;
    dim : natural32;
    restol : double_float := 1.0E-10;
    homtol : double_float := 1.0E-6;
    timer : Timing_Widget;

  begin
    if vrblvl > 0 then
      put("-> in main_factorization.");
      put_line("QuadDobl_Homotopy_Membership_Test ...");
    end if;
    QuadDobl_Read_Embedding(lp,genpts,dim);
    new_line;
    put_line("The input format of the test points is the solutions format.");
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Tune_Member_Tolerances(restol,homtol);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Homotopy_Membership_Test(file,lp.all,dim,genpts,sols,restol,homtol);
      tstop(timer);
      new_line(file);
      print_times(file,timer,"homotopy membership test");
    else
      QuadDobl_Multitasked_Membership_Test
        (file,nt,dim,homtol,lp.all,genpts,sols);
    end if;
    new_line(file);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end QuadDobl_Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test
              ( nt : in natural32; vrblvl : in integer32 := 0 ) is

    ans : character;

  begin
    if vrblvl > 0
     then put_line("-> in main_factorization.Homotopy_Membership_Test ...");
    end if;
    new_line;
    put_line("Membership test with homotopy :");
    put_line("  Input : embedded polynomial system with generic points, and");
    put_line("          list of test points.");
    put_line("  Output : decision whether test point lies on component.");
    new_line;
    put_line("MENU to choose the precision : ");
    put_line("  0. standard double precision homotopy continuation; or");
    put_line("  1. double double precision homotopy continuation; or");
    put_line("  2. quad double precision homotopy continuation.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Homotopy_Membership_Test(nt,vrblvl-1);
      when '1' => DoblDobl_Homotopy_Membership_Test(nt,vrblvl-1);
      when '2' => QuadDobl_Homotopy_Membership_Test(nt,vrblvl-1);
      when others => null;
    end case;
  end Homotopy_Membership_Test;

  procedure Trace_Form_Interpolation ( vrblvl : in integer32 := 0 ) is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    ans : character;

  begin
    if vrblvl > 0
     then put_line("-> in main_factorization.Trace_Form_Interpolation ...");
    end if;
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("MENU to certify monodromy breakup with by interpolation :");
    put_line
      ("  1. on given decomposition : use bootstrapping Newton to certify;");
    put_line
      ("  2.                        : use full trace form to certify;");
    put_line
      ("  3.                        : use Newton identities on trace form;");
    put_line
      ("  4.                        : use linear trace only to certify.");
    put("Type 1, 2, 3, or 4 to make your choice : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Call_Newton_Interpolate(file,lp.all,sols,dim);
      when '2' => Call_Trace_Interpolate(file,lp.all,sols,dim);
      when '3' => Call_Power_Trace_Interpolate(file,lp.all,sols,dim);
      when '4' => Call_Linear_Trace_Interpolate(file,lp.all,sols,dim);
      when others => null;
    end case;
  end Trace_Form_Interpolation;

  procedure Incremental_Interpolation ( vrblvl : in integer32 := 0 ) is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32;
    ans : character;

  begin
    if vrblvl > 0
     then put_line("-> in main_factorization.Incremental_Interpolation ...");
    end if;
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("MENU to decompose with incremental use of interpolation :");
    put_line("  1. massive interpolate with standard arithmetic;");
    put_line("  2. incremental interpolate with standard arithmetic;");
    put_line("  3.  + determine span with standard arithmetic;");
    put_line("  4.  + use central projections;");
    put_line("  5. massive interpolate with multi-precision arithmetic;");
    put_line("  6. incremental interpolate with multi-precision arithmetic;");
    put_line("  7.  + determine span with multi-precision arithmetic;");
    put_line("  8.  + use central projections;");
    Ask_Alternative(ans,"12345678");
    case ans is
      when '1' => Call_Standard_Interpolate(file,lp.all,sols,dim,0);
      when '2' => Call_Standard_Interpolate(file,lp.all,sols,dim,1);
      when '3' => Call_Standard_Interpolate(file,lp.all,sols,dim,2);
      when '4' => Call_Standard_Interpolate(file,lp.all,sols,dim,3);
      when '5' => Call_Multprec_Interpolate(file,lp.all,sols,dim,0);
      when '6' => Call_Multprec_Interpolate(file,lp.all,sols,dim,1);
      when '7' => Call_Multprec_Interpolate(file,lp.all,sols,dim,2);
      when '8' => Call_Multprec_Interpolate(file,lp.all,sols,dim,3);
      when others => null;
    end case;
  end Incremental_Interpolation;

  procedure Standard_Eliminate
                ( file : in file_type;
                  p : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in Standard_Complex_Solutions.Solution_List;
                  dim : in integer32 ) is

    use Standard_Stacked_Sample_Grids,Standard_Trace_Interpolators;

    sli : constant Standard_Complex_VecVecs.VecVec := Slices(p,natural32(dim));
    sps : constant Standard_Sample_List := Create(sols,sli);
    deg : constant integer32 := integer32(Length_Of(sps));
    ip : Standard_Complex_Polynomials.Poly;
 
  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    if dim = 1 then
      declare
        grid : Array_of_Standard_Sample_Lists(0..deg);
        eps,dst : double_float;
        t : Trace_Interpolator1;  
      begin
        Standard_Rectangular_Grid_Creator
          (file,sps,natural32(deg),grid,eps,dst);
        put(file,"Maximal error of the samples on the grid : ");
        put(file,eps,3); new_line(file);
        put(file,"Minimal distance between samples in one list in grid :");
        put(file,dst,3); new_line(file);
        t := Create(file,grid);
        put(file,"Maximal residual of evaluation at the grid : ");
        put(file,Maximal_Error(t,grid),3); new_line(file);    
        ip := Expand(t);
      end;
    else
      declare
        grid : Stacked_Sample_Grid(dim,deg);
        t : Trace_Interpolator;
      begin
        Standard_Stacked_Grid_Creator(file,sps,true,grid); 
        t := Create(file,grid,deg);
        put(file,"Maximal residual of evaluation at the grid : ");
        put(file,Maximal_Error(t,grid),3); new_line(file);    
        ip := Expand(t);
      end;
    end if;
    put_line(file,"The trace interpolator expanded as polynomial : ");
    put_line(file,ip);
  end Standard_Eliminate;

  procedure Multprec_Eliminate
                ( file : in file_type;
                  ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                  mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  sols : in Standard_Complex_Solutions.Solution_List;
                  dim : in integer32; size : in natural32 ) is       

    use Multprec_Stacked_Sample_Grids,Multprec_Trace_Interpolators;

    sli : constant Standard_Complex_VecVecs.VecVec := Slices(ep,natural32(dim));
    sps : Standard_Sample_List := Create(sols,sli);
    deg : constant integer32 := integer32(Length_Of(sps));
    ip : Multprec_Complex_Polynomials.Poly;

  begin
    Sampling_Machine.Initialize(ep,mp,dim,size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    new_line;
    put_line("See the output file for results...");
    new_line;                                                
    if dim = 1 then
      declare
        grid : Array_of_Multprec_Sample_Lists(0..deg);
        eps,dst : double_float;
        t : Trace_Interpolator1;
      begin
        Multprec_Rectangular_Grid_Creator
          (file,sps,natural32(deg),size,grid,eps,dst);
        put(file,"Maximal error of the samples on the grid : ");
        put(file,eps,3); new_line(file);
        put(file,"Minimal distance between samples in one list in grid :");
        put(file,dst,3); new_line(file);
        t := Create(file,grid);
        ip := Expand(t);
      end;
    else
      declare
        grid : Stacked_Sample_Grid(dim,deg);
        t : Trace_Interpolator;
      begin
        Multprec_Stacked_Grid_Creator(file,sps,true,size,grid);
        t := Create(file,grid,deg);
        ip := Expand(t);
      end;
    end if;
    put_line(file,"The trace interpolator expanded as polynomial : ");
    put_line(file,ip);
  end Multprec_Eliminate;

  procedure Eliminate_Variables ( vrblvl : in integer32 := 0 ) is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim,deci,size : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in main_factorization.Eliminate_Variables ...");
    end if;
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    Determine_Order(lp.all,sols);
    new_line;
    put("Give the number of decimal places (<= 16 is standard) : ");
    get(deci);
    new_line;
    if deci > 16 then
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Get_Multprec_System(lp.all,mp,size,dim);
      Multprec_Eliminate(file,lp.all,mp.all,sols,integer32(dim),size);
    else
      Standard_Eliminate(file,lp.all,sols,integer32(dim));
    end if;
  end Eliminate_Variables;

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   vrblvl : in integer32 := 0 ) is

    ans : character;

  begin
    if vrblvl > 0
     then put_line("-> in main_factorization.Main ...");
    end if;
    new_line;
    put_line("MENU to filter junk, factor components, and eliminate :");
    put_line("  0. filter solution lists subject to criteria;");
    put_line("  1. filter junk with homotopy membership test;");
    put_line("  2. breakup a filtered witness point set into irreducibles;");
    put_line("  3. given partition of breakup, compute trace form of filter;");
    put_line("  4. perform tasks 1, 2, and 3 by incremental interpolation;");
    put_line("  5. eliminate variables by interpolation via special slices;");
    put_line("  6. factor a complex polynomial in several variables;");
    put_line("  7. detect a common factor of two Laurent polynomials;");
    put_line("  8. filtered witness set breakup in double double precision;");
    put_line("  9. filtered witness set breakup in quad double precision.");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, or 9 to select a task : ");
    Ask_Alternative(ans,"0123456789");
    case ans is
      when '0' => Main_Solution_Filters.Main(infilename,outfilename);
      when '1' => Homotopy_Membership_Test(nt,vrblvl-1);
      when '2' => Standard_Breakup(infilename,outfilename);
      when '3' => Trace_Form_Interpolation(vrblvl-1);
      when '4' => Incremental_Interpolation(vrblvl-1);
      when '5' => Eliminate_Variables(vrblvl-1);
      when '6' =>
        if infilename = ""
         then Driver_to_Factor_Polynomial;
         else Driver_to_Factor_Polynomial(infilename);
        end if;
      when '7' => Driver_for_Common_Factor;
      when '8' => DoblDobl_Breakup(infilename,outfilename);
      when '9' => QuadDobl_Breakup(infilename,outfilename);
      when others => null;
    end case;
  end Main;

end Main_Factorization;
