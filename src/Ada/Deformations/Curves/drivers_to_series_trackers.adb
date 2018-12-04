with Communications_with_User;           use Communications_with_User;
with Characters_and_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers_cv;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Root_Refining_Parameters;
with Standard_Root_Refiners;
with DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Series_and_Homotopies;
with Series_and_Trackers;

package body Drivers_to_Series_Trackers is

  procedure Standard_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number ) is

    start : constant Standard_Complex_Poly_Systems.Poly_Sys
          := Standard_Homotopy.Start_System;
    target : constant Standard_Complex_Poly_Systems.Poly_Sys
           := Standard_Homotopy.Target_System;
    q : Standard_Complex_Poly_Systems.Poly_Sys(start'range);
    p : Standard_Complex_Poly_Systems.Poly_Sys(target'range);

  begin
    Standard_Complex_Poly_Systems.Copy(start,q);
    Standard_Complex_Poly_Systems.Copy(target,p);
    Standard_Homotopy.Clear;
    Standard_Homotopy.Create(p,q,1,gamma);
  end Standard_Reset_Gamma;

  procedure DoblDobl_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number ) is

    start : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
          := DoblDobl_Homotopy.Start_System;
    target : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
           := DoblDobl_Homotopy.Target_System;
    q : DoblDobl_Complex_Poly_Systems.Poly_Sys(start'range);
    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(target'range);
    ddgamma : constant DoblDobl_Complex_Numbers.Complex_Number
            := DoblDobl_Complex_Numbers_cv.Standard_to_DoblDobl_Complex(gamma);

  begin
    DoblDobl_Complex_Poly_Systems.Copy(start,q);
    DoblDobl_Complex_Poly_Systems.Copy(target,p);
    DoblDobl_Homotopy.Clear;
    DoblDobl_Homotopy.Create(p,q,1,ddgamma);
  end DoblDobl_Reset_Gamma;

  procedure QuadDobl_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number ) is

    start : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
          := QuadDobl_Homotopy.Start_System;
    target : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
           := QuadDobl_Homotopy.Target_System;
    q : QuadDobl_Complex_Poly_Systems.Poly_Sys(start'range);
    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(target'range);
    qdgamma : constant QuadDobl_Complex_Numbers.Complex_Number
            := QuadDobl_Complex_Numbers_cv.Standard_to_QuadDobl_Complex(gamma);

  begin
    QuadDobl_Complex_Poly_Systems.Copy(start,q);
    QuadDobl_Complex_Poly_Systems.Copy(target,p);
    QuadDobl_Homotopy.Clear;
    QuadDobl_Homotopy.Create(p,q,1,qdgamma);
  end QuadDobl_Reset_Gamma;

  procedure Set_Output
              ( file : in out file_type; verbose,tofile : out boolean ) is

    ans : character;

  begin
    new_line;
    put("Verbose?  Want to see extra output ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    put("Output to file ? (y/n) "); Ask_Yes_or_No(ans);
    tofile := (ans = 'y');
    if tofile then
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file);
    end if;
  end Set_Output;

  procedure Standard_Track
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);

  begin
    Series_and_Trackers.Track_Many_Paths(s,sols,pars);
    Standard_CSeries_Poly_Systems.Clear(s);
    Standard_Complex_Poly_Systems.Clear(h);
  end Standard_Track;

  procedure Standard_Track
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    Standard_Track(nq,sols,p);
  end Standard_Track;

  procedure Standard_Track
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false ) is

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);

  begin
    Series_and_Trackers.Track_Many_Paths(file,s,sols,pars,verbose);
    Standard_CSeries_Poly_Systems.Clear(s);
    Standard_Complex_Poly_Systems.Clear(h);
  end Standard_Track;

  procedure Standard_Track
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                verbose : in boolean := false ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    Standard_Track(file,nq,sols,p,verbose);
  end Standard_Track;

  procedure DoblDobl_Track
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);

  begin
    Series_and_Trackers.Track_Many_Paths(s,sols,pars);
    DoblDobl_CSeries_Poly_Systems.Clear(s);
    DoblDobl_Complex_Poly_Systems.Clear(h);
  end DoblDobl_Track;

  procedure DoblDobl_Track
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    DoblDobl_Track(nq,sols,p);
  end DoblDobl_Track;

  procedure DoblDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false ) is

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);

  begin
    Series_and_Trackers.Track_Many_Paths(file,s,sols,pars,verbose);
    DoblDobl_CSeries_Poly_Systems.Clear(s);
    DoblDobl_Complex_Poly_Systems.Clear(h);
  end DoblDobl_Track;

  procedure DoblDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in boolean := false ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    DoblDobl_Track(file,nq,sols,p,verbose);
  end DoblDobl_Track;

  procedure QuadDobl_Track
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);

  begin
    Series_and_Trackers.Track_Many_Paths(s,sols,pars);
    QuadDobl_CSeries_Poly_Systems.Clear(s);
    QuadDobl_Complex_Poly_Systems.Clear(h);
  end QuadDobl_Track;

  procedure QuadDobl_Track
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    QuadDobl_Track(nq,sols,p);
  end QuadDobl_Track;

  procedure QuadDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false ) is

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);

  begin
    Series_and_Trackers.Track_Many_Paths(file,s,sols,pars,verbose);
    QuadDobl_CSeries_Poly_Systems.Clear(s);
    QuadDobl_Complex_Poly_Systems.Clear(h);
  end QuadDobl_Track;

  procedure QuadDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in boolean := false ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    QuadDobl_Track(file,nq,sols,p,verbose);
  end QuadDobl_Track;

  procedure Write_Timer
              ( file : in file_type;
                numdeg,dendeg,precision : in natural32;
                timer : in Timing_Widget ) is

    s : constant string
      := "[" & Characters_and_Numbers.Convert(integer32(numdeg))
       & "," & Characters_and_Numbers.Convert(integer32(dendeg))
       & "]-Tracking";
 
  begin
    new_line(file);
    case precision is
      when 0 => print_times(file,timer,s & " in double precision.");
      when 1 => print_times(file,timer,s & " in double double precision.");
      when 2 => print_times(file,timer,s & " in quad double precision.");
      when others => null;
    end case;
  end Write_Timer;

  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Target_System;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout : boolean;

    use Root_Refining_Parameters,Standard_Root_Refiners;

  begin
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false;
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
  end Refine_Roots;

  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Target_System;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout : boolean;

    use Root_Refining_Parameters,DoblDobl_Root_Refiners;

  begin
    DoblDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false;
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
  end Refine_Roots;

  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Target_System;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout : boolean;

    use Root_Refining_Parameters,QuadDobl_Root_Refiners;

  begin
    QuadDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false;
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
  end Refine_Roots;

end Drivers_to_Series_trackers;
