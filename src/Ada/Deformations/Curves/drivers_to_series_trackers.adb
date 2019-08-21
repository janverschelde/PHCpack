with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
with Characters_and_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Hessians;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Hessians;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Hessians;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Projective_Transformations;
with Multi_Projective_Transformations;
with Affine_Transformations;
with Standard_Mixed_Residuals;
with DoblDobl_Mixed_Residuals;
with QuadDobl_Mixed_Residuals;
with Root_Refining_Parameters;
with Standard_Root_Refiners;
with DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Singular_Values_of_Hessians;
with Series_and_Homotopies;
with Series_and_Trackers;
with Write_Seed_Number;
with Greeting_Banners;

package body Drivers_to_Series_Trackers is

  procedure Standard_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number;
                tpow : in natural32 := 2 ) is

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
    Standard_Homotopy.Create(p,q,tpow,gamma);
  end Standard_Reset_Gamma;

  procedure DoblDobl_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number;
                tpow : in natural32 := 2 ) is

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
    DoblDobl_Homotopy.Create(p,q,tpow,ddgamma);
  end DoblDobl_Reset_Gamma;

  procedure QuadDobl_Reset_Gamma
              ( gamma : in Standard_Complex_Numbers.Complex_Number;
                tpow : in natural32 := 2 ) is

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
    QuadDobl_Homotopy.Create(p,q,tpow,qdgamma);
  end QuadDobl_Reset_Gamma;

  procedure Set_Output
              ( file : in out file_type;
                monitor,verbose,tofile : out boolean ) is

    ans : character;

  begin
    new_line;
    put("Monitor progress of the path tracker ? (y/n) ");
    Ask_Yes_or_No(ans);
    monitor := (ans = 'y');
    put("Verbose?  Want to see extra output ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    put("Output to file ? (y/n) "); Ask_Yes_or_No(ans);
    tofile := (ans = 'y');
    if tofile then
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file);
      if not monitor then
        new_line;
        put_line("See the output file for results ...");
        new_line;
      end if;
    end if;
  end Set_Output;

  procedure Standard_Track
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    jm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : Standard_Complex_Hessians.Link_to_Array_of_Hessians;

    use Singular_Values_of_Hessians;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.Standard_Track 1 ...");
    end if;
    Standard_Jacobian_Hessians_of_Homotopy(jm,hs);
    Series_and_Trackers.Track_Many_Paths(jm,hs,s,sols,pars,mhom,idz,vrblvl-1);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Hessians.Clear(hs);
    Standard_CSeries_Poly_Systems.Clear(s);
    Standard_Complex_Poly_Systems.Clear(h);
  end Standard_Track;

  procedure Standard_Track
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.Standard_Track 2 ...");
    end if;
    Standard_Track(nq,sols,p,mhom,idz,vrblvl-1);
  end Standard_Track;

  procedure Standard_Track
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    jm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : Standard_Complex_Hessians.Link_to_Array_of_Hessians;

    use Singular_Values_of_Hessians;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.Standard_Track 3 ...");
    end if;
    Standard_Jacobian_Hessians_of_Homotopy(jm,hs);
    if verbose then
      Series_and_Trackers.Track_Many_Paths
        (file,jm,hs,s,sols,pars,mhom,idz,true,true,vrblvl-1);
    else
      Series_and_Trackers.Track_Many_Paths
        (file,jm,hs,s,sols,pars,mhom,idz,vrblvl=>vrblvl-1);
    end if;
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Hessians.Clear(hs);
    Standard_CSeries_Poly_Systems.Clear(s);
    Standard_Complex_Poly_Systems.Clear(h);
  end Standard_Track;

  procedure Standard_Track
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.Standard_Track 4 ...");
    end if;
    Standard_Track(file,nq,sols,p,mhom,idz,verbose,vrblvl-1);
  end Standard_Track;

  procedure DoblDobl_Track
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;

    use Singular_Values_of_Hessians;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.DoblDobl_Track 1 ...");
    end if;
    DoblDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
    Series_and_Trackers.Track_Many_Paths(jm,hs,s,sols,pars,mhom,idz,vrblvl-1);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Hessians.Clear(hs);
    DoblDobl_CSeries_Poly_Systems.Clear(s);
    DoblDobl_Complex_Poly_Systems.Clear(h);
  end DoblDobl_Track;

  procedure DoblDobl_Track
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.DoblDobl_Track 2 ...");
    end if;
    DoblDobl_Track(nq,sols,p,mhom,idz,vrblvl-1);
  end DoblDobl_Track;

  procedure DoblDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;

    use Singular_Values_of_Hessians;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.DoblDobl_Track 3 ...");
    end if;
    DoblDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
    if verbose then
      Series_and_Trackers.Track_Many_Paths
        (file,jm,hs,s,sols,pars,mhom,idz,true,true,vrblvl-1);
    else
      Series_and_Trackers.Track_Many_Paths
        (file,jm,hs,s,sols,pars,mhom,idz,vrblvl=>vrblvl-1);
    end if;
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Hessians.Clear(hs);
    DoblDobl_CSeries_Poly_Systems.Clear(s);
    DoblDobl_Complex_Poly_Systems.Clear(h);
  end DoblDobl_Track;

  procedure DoblDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.DoblDobl_Track 4 ...");
    end if;
    DoblDobl_Track(file,nq,sols,p,mhom,idz,verbose,vrblvl-1);
  end DoblDobl_Track;

  procedure QuadDobl_Track
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    jm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;

    use Singular_Values_of_Hessians;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.QuadDobl_Track 1 ...");
    end if;
    QuadDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
    Series_and_Trackers.Track_Many_Paths(jm,hs,s,sols,pars,mhom,idz,vrblvl-1);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Hessians.Clear(hs);
    QuadDobl_CSeries_Poly_Systems.Clear(s);
    QuadDobl_Complex_Poly_Systems.Clear(h);
  end QuadDobl_Track;

  procedure QuadDobl_Track
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.QuadDobl_Track 2 ...");
    end if;
    QuadDobl_Track(nq,sols,p,mhom,idz,vrblvl-1);
  end QuadDobl_Track;

  procedure QuadDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    jm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;

    use Singular_Values_of_Hessians;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.QuadDobl_Track 3 ...");
    end if;
    QuadDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
    if verbose then
      Series_and_Trackers.Track_Many_Paths
        (file,jm,hs,s,sols,pars,mhom,idz,true,true,vrblvl-1);
    else
      Series_and_Trackers.Track_Many_Paths
        (file,jm,hs,s,sols,pars,mhom,idz,vrblvl=>vrblvl-1);
    end if;
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Hessians.Clear(hs);
    QuadDobl_CSeries_Poly_Systems.Clear(s);
    QuadDobl_Complex_Poly_Systems.Clear(h);
  end QuadDobl_Track;

  procedure QuadDobl_Track
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    p : constant Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;

  begin
    if vrblvl > 0 then
      put_line("-> in drivers_to_series_trackers.QuadDobl_Track 4 ...");
    end if;
    QuadDobl_Track(file,nq,sols,p,mhom,idz,verbose,vrblvl-1);
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

  procedure Write_Conclusion 
              ( file : in file_type; start_moment : in Ada.Calendar.Time ) is

    ended_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

  begin
    new_line(file);
    put(file,"PHC ran from ");
    Time_Stamps.Write_Time_Stamp(file,start_moment);
    put(file," till ");
    Time_Stamps.Write_Time_Stamp(file,ended_moment);
    put_line(file,".");
    Time_Stamps.Write_Elapsed_Time(file,start_moment,ended_moment);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end Write_Conclusion;

  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Target_System;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout : boolean;

    use Root_Refining_Parameters,Standard_Root_Refiners;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_to_series_trackers.Refine_Roots 1 ...");
    end if;
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false;
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
  end Refine_Roots;

  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Target_System;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout : boolean;

    use Root_Refining_Parameters,DoblDobl_Root_Refiners;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_to_series_trackers.Refine_Roots 2 ...");
    end if;
    DoblDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false;
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
  end Refine_Roots;

  procedure Refine_Roots
              ( file : in file_type; nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Target_System;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout : boolean;

    use Root_Refining_Parameters,QuadDobl_Root_Refiners;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_to_series_trackers.Refine_Roots 3 ...");
    end if;
    QuadDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false;
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
  end Refine_Roots;

  procedure Refine_Roots
              ( file : in file_type;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use Root_Refining_Parameters,Standard_Root_Refiners;

    p : constant Standard_Complex_Poly_Systems.Poly_Sys
      := Standard_Homotopy.Target_System;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout : boolean;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_to_series_trackers.Refine_Roots 4 ...");
    end if;
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false;
    if mhom = 0 then
      Reporting_Root_Refiner
        (file,p,abh,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
    else -- dehomogenize p and the solutions
      if mhom = 1 then
        Projective_Transformations.Affine_Transformation(sols);
        declare
          q : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last-1);
          abq : Standard_Complex_Poly_Systems.Poly_Sys(q'range);
          evabq : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(q'range);
        begin
          q := Affine_Transformations.Make_Affine(p);
          abq := Standard_Mixed_Residuals.AbsVal(q);
          evabq := Standard_Complex_Poly_SysFun.Create(abq);
          new_line(file);
          put(file,natural32(q'last),q);
          Reporting_Root_Refiner
            (file,q,evabq,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
        end;
      else
        Multi_Projective_Transformations.Make_Affine(sols,mhom,idz.all);
        declare
          im : constant integer32 := integer32(mhom);
          q : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last-im);
          abq : Standard_Complex_Poly_Systems.Poly_Sys(q'range);
          evabq : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(q'range);
        begin
          q := Affine_Transformations.Make_Affine(p,mhom);
          abq := Standard_Mixed_Residuals.AbsVal(q);
          evabq := Standard_Complex_Poly_SysFun.Create(abq);
          new_line(file);
          put(file,natural32(q'last),q);
          Reporting_Root_Refiner
            (file,q,evabq,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
        end;
      end if;
    end if;
  end Refine_Roots;

  procedure Refine_Roots
              ( file : in file_type;
                abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
      := DoblDobl_Homotopy.Target_System;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout : boolean;

    use Root_Refining_Parameters,DoblDobl_Root_Refiners;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_to_series_trackers.Refine_Roots 5 ...");
    end if;
    DoblDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false;
    if mhom = 0 then
      Reporting_Root_Refiner
        (file,p,abh,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
    else -- dehomogenize p and the solutions
      if mhom = 1 then
        Projective_Transformations.Affine_Transformation(sols);
        declare
          q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last-1);
          abq : DoblDobl_Complex_Poly_Systems.Poly_Sys(q'range);
          evabq : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(q'range);
        begin
          q := Affine_Transformations.Make_Affine(p);
          abq := DoblDobl_Mixed_Residuals.AbsVal(q);
          evabq := DoblDobl_Complex_Poly_SysFun.Create(abq);
          new_line(file);
          put(file,natural32(q'last),q);
          Reporting_Root_Refiner
            (file,q,evabq,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
        end;
      else
        Multi_Projective_Transformations.Make_Affine(sols,mhom,idz.all);
        declare
          im : constant integer32 := integer32(mhom);
          q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last-im);
          abq : DoblDobl_Complex_Poly_Systems.Poly_Sys(q'range);
          evabq : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(q'range);
        begin
          q := Affine_Transformations.Make_Affine(p,mhom);
          abq := DoblDobl_Mixed_Residuals.AbsVal(q);
          evabq := DoblDobl_Complex_Poly_SysFun.Create(abq);
          new_line(file);
          put(file,natural32(q'last),q);
          Reporting_Root_Refiner
            (file,q,evabq,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
        end;
      end if;
    end if;
  end Refine_Roots;

  procedure Refine_Roots
              ( file : in file_type;
                abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use Root_Refining_Parameters,QuadDobl_Root_Refiners;

    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
      := QuadDobl_Homotopy.Target_System;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout : boolean;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_to_series_trackers.Refine_Roots 6 ...");
    end if;
    QuadDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false;
    if mhom = 0 then
      Reporting_Root_Refiner
        (file,p,abh,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
    else -- dehomogenize p and the solutions
      if mhom = 1 then
        Projective_Transformations.Affine_Transformation(sols);
        declare
          q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last-1);
          abq : QuadDobl_Complex_Poly_Systems.Poly_Sys(q'range);
          evabq : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(q'range);
        begin
          q := Affine_Transformations.Make_Affine(p);
          abq := QuadDobl_Mixed_Residuals.AbsVal(q);
          evabq := QuadDobl_Complex_Poly_SysFun.Create(abq);
          new_line(file);
          put(file,natural32(q'last),q);
          Reporting_Root_Refiner
            (file,q,evabq,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
        end;
      else
        Multi_Projective_Transformations.Make_Affine(sols,mhom,idz.all);
        declare
          im : constant integer32 := integer32(mhom);
          q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last-im);
          abq : QuadDobl_Complex_Poly_Systems.Poly_Sys(q'range);
          evabq : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(q'range);
        begin
          q := Affine_Transformations.Make_Affine(p,mhom);
          abq := QuadDobl_Mixed_Residuals.AbsVal(q);
          evabq := QuadDobl_Complex_Poly_SysFun.Create(abq);
          new_line(file);
          put(file,natural32(q'last),q);
          Reporting_Root_Refiner
            (file,q,evabq,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
        end;
      end if;
    end if;
  end Refine_Roots;

end Drivers_to_Series_trackers;
