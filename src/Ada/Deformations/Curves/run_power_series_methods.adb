with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_System_and_Solutions_io;
with TripDobl_Complex_Polynomials;
-- with TripDobl_Complex_Poly_Systems_io;   use TripDobl_Complex_Poly_Systems_io;
with TripDobl_System_and_Solutions_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_System_and_Solutions_io;
with PentDobl_Complex_Polynomials;
-- with PentDobl_Complex_Poly_Systems_io;   use PentDobl_Complex_Poly_Systems_io;
with PentDobl_System_and_Solutions_io;
with OctoDobl_Complex_Polynomials;
-- with OctoDobl_Complex_Poly_Systems_io;   use OctoDobl_Complex_Poly_Systems_io;
with OctoDobl_System_and_Solutions_io;
with DecaDobl_Complex_Polynomials;
-- with DecaDobl_Complex_Poly_Systems_io;   use DecaDobl_Complex_Poly_Systems_io;
with DecaDobl_System_and_Solutions_io;
with Standard_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;
with Series_and_Solutions;
with Power_Series_Methods;               use Power_Series_Methods;

package body Run_Power_Series_Methods is

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 ) is

    nbrit,maxdeg : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    if vrb > 0
     then put_line("-> in run_power_series_methods.Run_Newton 1 ...");
    end if;
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Give the maximal degree of the series : "); get(maxdeg);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    if echelon then
      new_line;
      put_line("Echelon Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_Echelon_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    else
      new_line;
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"power series Newton in double precision");
  end Run_Newton;

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 ) is

    maxdeg,nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    if vrb > 0
     then put_line("-> in run_power_series_methods.Run_Newton 2 ...");
    end if;
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Give the maximal degree of the series : "); get(maxdeg);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    if echelon then
      new_line;
      put_line("Echelon Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_Echelon_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    else
      new_line;
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"power series Newton in double double precision");
  end Run_Newton;

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in TripDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 ) is

    maxdeg,nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    if vrb > 0
     then put_line("-> in run_power_series_methods.Run_Newton 3 ...");
    end if;
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Give the maximal degree of the series : "); get(maxdeg);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    if echelon then
      new_line;
      put_line("Echelon Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
     -- Run_Echelon_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    else
      new_line;
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"power series Newton in triple double precision");
  end Run_Newton;

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 ) is

    maxdeg,nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    if vrb > 0
     then put_line("-> in run_power_series_methods.Run_Newton 4 ...");
    end if;
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Give the maximal degree of the series : "); get(maxdeg);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    if echelon then
      new_line;
      put_line("Echelon Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_Echelon_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    else
      new_line;
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"power series Newton in quad double precision");
  end Run_Newton;

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in PentDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 ) is

    maxdeg,nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    if vrb > 0
     then put_line("-> in run_power_series_methods.Run_Newton 5 ...");
    end if;
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Give the maximal degree of the series : "); get(maxdeg);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    if echelon then
      new_line;
      put_line("Echelon Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
     -- Run_Echelon_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    else
      new_line;
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"power series Newton in penta double precision");
  end Run_Newton;

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in OctoDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 ) is

    maxdeg,nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    if vrb > 0
     then put_line("-> in run_power_series_methods.Run_Newton 6 ...");
    end if;
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Give the maximal degree of the series : "); get(maxdeg);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    if echelon then
      new_line;
      put_line("Echelon Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
     -- Run_Echelon_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    else
      new_line;
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"power series Newton in octo double precision");
  end Run_Newton;

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in DecaDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 ) is

    maxdeg,nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    if vrb > 0
     then put_line("-> in run_power_series_methods.Run_Newton 7 ...");
    end if;
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Give the maximal degree of the series : "); get(maxdeg);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    if echelon then
      new_line;
      put_line("Echelon Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
     -- Run_Echelon_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    else
      new_line;
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"power series Newton in deca double precision");
  end Run_Newton;

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant Standard_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : Standard_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.Run_Newton_at_Constant 1 ...");
    end if;
    Run_Newton(file,false,srp,srv,vrb-1);
    Standard_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant DoblDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.Run_Newton_at_Constant 2 ...");
    end if;
    Run_Newton(file,false,srp,srv,vrb-1);
    DoblDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in TripDobl_Complex_Poly_Systems.Poly_Sys;
               s : in TripDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 ) is

    use TripDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant TripDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : TripDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.Run_Newton_at_Constant 3 ...");
    end if;
    Run_Newton(file,false,srp,srv,vrb-1);
    TripDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant QuadDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.Run_Newton_at_Constant 4 ...");
    end if;
    Run_Newton(file,false,srp,srv,vrb-1);
    QuadDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in PentDobl_Complex_Poly_Systems.Poly_Sys;
               s : in PentDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 ) is

    use PentDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant PentDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : PentDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.Run_Newton_at_Constant 5 ...");
    end if;
    Run_Newton(file,false,srp,srv,vrb-1);
    PentDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in OctoDobl_Complex_Poly_Systems.Poly_Sys;
               s : in OctoDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 ) is

    use OctoDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant OctoDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : OctoDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.Run_Newton_at_Constant 6 ...");
    end if;
    Run_Newton(file,false,srp,srv,vrb-1);
    OctoDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in DecaDobl_Complex_Poly_Systems.Poly_Sys;
               s : in DecaDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 ) is

    use DecaDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant DecaDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : DecaDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.Run_Newton_at_Constant 7 ...");
    end if;
    Run_Newton(file,false,srp,srv,vrb-1);
    DecaDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Standard_Main_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.Standard_Main_at_Constant ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      Standard_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      Standard_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton_at_Constant(outfile,idx,lp.all,sols,vrb-1);
  end Standard_Main_at_Constant;

  procedure DoblDobl_Main_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.DoblDobl_Main_at_Constant ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      DoblDobl_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      DoblDobl_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton_at_Constant(outfile,idx,lp.all,sols,vrb-1);
  end DoblDobl_Main_at_Constant;

  procedure TripDobl_Main_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use TripDobl_Complex_Polynomials;
    use TripDobl_Complex_Poly_Systems;
    use TripDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.TripDobl_Main_at_Constant ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      TripDobl_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      TripDobl_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton_at_Constant(outfile,idx,lp.all,sols,vrb-1);
  end TripDobl_Main_at_Constant;

  procedure QuadDobl_Main_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.QuadDobl_Main_at_Constant ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      QuadDobl_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      QuadDobl_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton_at_Constant(outfile,idx,lp.all,sols,vrb-1);
  end QuadDobl_Main_at_Constant;

  procedure PentDobl_Main_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use PentDobl_Complex_Polynomials;
    use PentDobl_Complex_Poly_Systems;
    use PentDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.PentDobl_Main_at_Constant ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      PentDobl_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      PentDobl_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton_at_Constant(outfile,idx,lp.all,sols,vrb-1);
  end PentDobl_Main_at_Constant;

  procedure OctoDobl_Main_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use OctoDobl_Complex_Polynomials;
    use OctoDobl_Complex_Poly_Systems;
    use OctoDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.OctoDobl_Main_at_Constant ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      OctoDobl_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      OctoDobl_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton_at_Constant(outfile,idx,lp.all,sols,vrb-1);
  end OctoDobl_Main_at_Constant;

  procedure DecaDobl_Main_at_Constant
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use DecaDobl_Complex_Polynomials;
    use DecaDobl_Complex_Poly_Systems;
    use DecaDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.DecaDobl_Main_at_Constant ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      DecaDobl_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      DecaDobl_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton_at_Constant(outfile,idx,lp.all,sols,vrb-1);
  end DecaDobl_Main_at_Constant;

  procedure Standard_Main_at_Series
               ( infilename,outfilename : in string;
                 vrb : in integer32 := 0 ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    srv : Standard_Complex_Series_Vectors.Link_to_Vector;

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.Standard_Main_at_Series ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading a polynomial system ...");
      get(lp);
    else
      Open_Input_File(infile,infilename);
      get(infile,lp);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    new_line;
    put_line("Reading a series to start Newton's method at ...");
    Complex_Series_and_Polynomials_io.get(srv);
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    declare
      s : Standard_Complex_Series_VecVecs.VecVec(1..1);
      srp : Standard_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(outfile,true,srp,s,vrb-1);
      Standard_CSeries_Poly_Systems.Clear(srp);
    end;
  end Standard_Main_at_Series;

  procedure DoblDobl_Main_at_Series
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    srv : DoblDobl_Complex_Series_Vectors.Link_to_Vector;

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.DoblDobl_Main_at_Series ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading a polynomial system ...");
      get(lp);
    else
      Open_Input_File(infile,infilename);
      get(infile,lp);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    new_line;
    put_line("Reading a series to start Newton's method at ...");
    Complex_Series_and_Polynomials_io.get(srv);
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    declare
      s : DoblDobl_Complex_Series_VecVecs.VecVec(1..1);
      srp : DoblDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(outfile,true,srp,s,vrb-1);
      DoblDobl_CSeries_Poly_Systems.Clear(srp);
    end;
  end DoblDobl_Main_at_Series;

  procedure QuadDobl_Main_at_Series
              ( infilename,outfilename : in string;
                vrb : in integer32 := 0 ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    srv : QuadDobl_Complex_Series_Vectors.Link_to_Vector;

  begin
    if vrb > 0 then
      put_line("-> in run_power_series_methods.QuadDobl_Main_at_Series ...");
    end if;
    if infilename = "" then
      new_line;
      put_line("Reading a polynomial system ...");
      get(lp);
    else
      Open_Input_File(infile,infilename);
      get(infile,lp);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    new_line;
    put_line("Reading a series to start Newton's method at ...");
    Complex_Series_and_Polynomials_io.get(srv);
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    declare
      s : QuadDobl_Complex_Series_VecVecs.VecVec(1..1);
      srp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(outfile,true,srp,s,vrb-1);
      QuadDobl_CSeries_Poly_Systems.Clear(srp);
    end;
  end QuadDobl_Main_at_Series;

end Run_Power_Series_Methods;
