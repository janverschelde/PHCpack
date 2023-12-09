with Ada.Calendar;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;
with TripDobl_Complex_Numbers;
with TripDobl_Complex_Numbers_cv;        use TripDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;        use QuadDobl_Complex_Numbers_cv;
with PentDobl_Complex_Numbers;
with PentDobl_Complex_Numbers_cv;        use PentDobl_Complex_Numbers_cv;
with OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers_cv;        use OctoDobl_Complex_Numbers_cv;
with DecaDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers_cv;        use DecaDobl_Complex_Numbers_cv;
with HexaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers_cv;        use HexaDobl_Complex_Numbers_cv;
with Characters_and_Numbers;
with Numbers_io;                         use Numbers_io;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Hessians;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Hessians;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with TripDobl_Complex_Poly_Systems;
with TripDobl_Complex_Poly_Systems_io;   use TripDobl_Complex_Poly_Systems_io;
with TripDobl_System_and_Solutions_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Hessians;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with PentDobl_Complex_Poly_Systems;
with PentDobl_Complex_Poly_Systems_io;   use PentDobl_Complex_Poly_Systems_io;
with PentDobl_System_and_Solutions_io;
with OctoDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Poly_Systems_io;   use OctoDobl_Complex_Poly_Systems_io;
with OctoDobl_System_and_Solutions_io;
with DecaDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Poly_Systems_io;   use DecaDobl_Complex_Poly_Systems_io;
with DecaDobl_System_and_Solutions_io;
with HexaDobl_Complex_Poly_Systems;
with HexaDobl_Complex_Poly_Systems_io;   use HexaDobl_Complex_Poly_Systems_io;
with HexaDobl_System_and_Solutions_io;
with Solution_Drops;
with Standard_Homotopy;
with Standard_Coefficient_Homotopy;
with DoblDobl_Homotopy;
with DoblDobl_Coefficient_Homotopy;
with TripDobl_Homotopy;
with TripDobl_Coefficient_Homotopy;
with QuadDobl_Homotopy;
with QuadDobl_Coefficient_Homotopy;
with PentDobl_Homotopy;
with PentDobl_Coefficient_Homotopy;
with OctoDobl_Homotopy;
with OctoDobl_Coefficient_Homotopy;
with DecaDobl_Homotopy;
with DecaDobl_Coefficient_Homotopy;
with HexaDobl_Homotopy;
with HexaDobl_Coefficient_Homotopy;
with Standard_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_SysFun;
with Standard_CSeries_Jaco_Matrices;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_SysFun;
with DoblDobl_CSeries_Jaco_Matrices;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_SysFun;
with QuadDobl_CSeries_Jaco_Matrices;
with Partitions_of_Sets_of_Unknowns_io;
with Singular_Values_of_Hessians;
with Series_and_Homotopies;
with Series_and_Trackers;
with Homotopy_Mixed_Residuals;
with Homotopy_Series_Readers;
with Homotopy_Continuation_Parameters_io;
with Standard_Pade_Trackers;
with DoblDobl_Pade_Trackers;
with QuadDobl_Pade_Trackers;
with Drivers_to_Series_Trackers;         use Drivers_to_Series_Trackers;

package body Series_Path_Trackers is

  procedure Standard_Write 
              ( file : in file_type; nq,nv : in natural32;
                idxpar : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is
  begin
    if idxpar /= 0 then
      put(file,nq,nv,Standard_Homotopy.Homotopy_System);
    else
      put(file,nq,Standard_Homotopy.Target_System);
      new_line(file);
      put_line(file,"THE START SYSTEM :");
      put(file,nq,Standard_Homotopy.Start_System);
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,Standard_Complex_Solutions.Length_Of(sols),
        natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars);
    new_line(file); flush(file);
  end Standard_Write;

  procedure DoblDobl_Write 
              ( file : in file_type; nq,nv : in natural32;
                idxpar : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is
  begin
    if idxpar /= 0 then
      put(file,nq,nv,DoblDobl_Homotopy.Homotopy_System);
    else
      put(file,nq,DoblDobl_Homotopy.Target_System);
      new_line(file);
      put_line(file,"THE START SYSTEM :");
      put(file,nq,DoblDobl_Homotopy.Start_System);
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
        natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars);
    new_line(file); flush(file);
  end DoblDobl_Write;

  procedure QuadDobl_Write 
              ( file : in file_type; nq,nv : in natural32;
                idxpar : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is
  begin
    if idxpar /= 0 then
      put(file,nq,nv,QuadDobl_Homotopy.Homotopy_System);
    else
      put(file,nq,QuadDobl_Homotopy.Target_System);
      new_line(file);
      put_line(file,"THE START SYSTEM :");
      put(file,nq,QuadDobl_Homotopy.Start_System);
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
        natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars);
    new_line(file); flush(file);
  end QuadDobl_Write;

  function Set_Dimension ( nvr,idxpar : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns nvr if idxpar is zero, otherwise returns nvr - 1.

  begin
    if idxpar = 0
     then return nvr;
     else return nvr-1;
    end if;
  end Set_Dimension;

  procedure Standard_Run
              ( monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq);
    abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(1..nq)
        := Homotopy_Mixed_Residuals.Standard_AbsVal_Homotopy;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq);
    fhm : Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(1..nq);
    fcf : Standard_Complex_Series_VecVecs.VecVec(1..nq);
    dim : constant integer32 := Set_Dimension(nvr,idxpar);
    jm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : Standard_Complex_Hessians.Link_to_Array_of_Hessians;
    ejm : Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(h'range,1..dim);
    mlt : Standard_CSeries_Jaco_Matrices.Mult_Factors(h'range,1..dim);
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    ans : character;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs,cntcut,cntfail : natural32;
    minsize,maxsize,smallest,largest : double_float;
    cntsstp,cntdstp,cntpstp : natural32;
    ratsstp,ratdstp,ratpstp : double_float := 0.0;
    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

    use Singular_Values_of_Hessians;

  begin
    if vrb > 0
     then put_line("-> in series_path_trackers.Standard_Run ...");
    end if;
    if idxpar /= 0
     then Standard_Jacobian_Hessians_of_Homotopy(idxpar,jm,hs);
     else Standard_Jacobian_Hessians_of_Homotopy(jm,hs);
    end if;
    h := Standard_Homotopy.Homotopy_System;
    if idxpar /= 0
     then s := Series_and_Homotopies.Create(h,idxpar,false);
     else s := Series_and_Homotopies.Create(h,nq+1,false);
    end if;
    fhm := Standard_CSeries_Poly_SysFun.Create(s);
    fcf := Standard_CSeries_Poly_SysFun.Coeff(s);
    Standard_CSeries_Jaco_Matrices.Create(s,ejm,mlt);
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (pars.maxsteps+1)*pars.corsteps+1; maxnbrcorrs := 0;
    smallest := pars.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put("Tracking path "); put(i,1); put_line(" ...");
      end if;
      Standard_Pade_Trackers.Track_One_Path
        (standard_output,abh,jm,hs,fhm,fcf,ejm,mlt,ls.all,pars,mhom,idz,
         nbrsteps,nbrcorrs,cntcut,cntfail,minsize,maxsize,
         cntsstp,cntdstp,cntpstp,verbose,vrb-1);
      if verbose then
        Series_and_Trackers.Write_Path_Statistics
          (standard_output,nbrsteps,nbrcorrs,cntcut,cntfail,minsize,maxsize,
           cntdstp,cntpstp);
      end if;
      put("Solution "); put(i,1); put_line(" :"); put(ls.all); new_line;
      put("Continue to the next path ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Series_and_Trackers.Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Series_and_Trackers.Update_MinMax(smallest,largest,minsize,maxsize);
      Series_and_Trackers.Update_Ratio_Sums(ratsstp,ratdstp,ratpstp,
        cntsstp,cntdstp,cntpstp,nbrsteps*natural32(len));
    end loop;
    tstop(timer);
    Series_and_Trackers.Write_Total_Path_Statistics
      (standard_output,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,
       smallest,largest,ratdstp,ratpstp);
    new_line;
    put_line("THE SOLUTIONS :");
    put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    Write_Timer(standard_output,pars.numdeg,pars.dendeg,0,timer);
    if idxpar = 0
     then Refine_Roots(standard_output,abh,mhom,idz,sols,vrb);
    end if;
    Write_Conclusion(standard_output,start_moment);
    Standard_Complex_Poly_Systems.Clear(h);
    Standard_Complex_Poly_SysFun.Clear(abh);
    Standard_CSeries_Poly_Systems.Clear(s);
    Standard_CSeries_Poly_SysFun.Clear(fhm);
    Standard_Complex_Series_VecVecs.Clear(fcf);
    Standard_CSeries_Jaco_Matrices.Clear(ejm);
    Standard_Cseries_Jaco_Matrices.Clear(mlt);
  end Standard_Run;

  procedure DoblDobl_Run
              ( monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq);
    abh : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(1..nq)
        := Homotopy_Mixed_Residuals.DoblDobl_AbsVal_Homotopy;
    s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq);
    jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
    fhm : DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(1..nq);
    fcf : DoblDobl_Complex_Series_VecVecs.VecVec(1..nq);
    dim : constant integer32 := Set_Dimension(nvr,idxpar);
    ejm : DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(h'range,1..dim);
    mlt : DoblDobl_CSeries_Jaco_Matrices.Mult_Factors(h'range,1..dim);
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    ans : character;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs,cntcut,cntfail : natural32;
    minsize,maxsize,smallest,largest : double_float;
    cntsstp,cntdstp,cntpstp : natural32;
    ratsstp,ratdstp,ratpstp : double_float := 0.0;
    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

    use Singular_Values_of_Hessians;

  begin
    if vrb > 0
     then put_line("-> in series_path_trackers.DoblDobl_Run ...");
    end if;
    if idxpar /= 0
     then DoblDobl_Jacobian_Hessians_of_Homotopy(idxpar,jm,hs);
     else DoblDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
    end if;
    h := DoblDobl_Homotopy.Homotopy_System;
    if idxpar /= 0
     then s := Series_and_Homotopies.Create(h,idxpar,false);
     else s := Series_and_Homotopies.Create(h,nq+1,false);
    end if;
    fhm := DoblDobl_CSeries_Poly_SysFun.Create(s);
    fcf := DoblDobl_CSeries_Poly_SysFun.Coeff(s);
    DoblDobl_CSeries_Jaco_Matrices.Create(s,ejm,mlt);
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (pars.maxsteps+1)*pars.corsteps+1; maxnbrcorrs := 0;
    smallest := pars.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put("Tracking path "); put(i,1); put_line(" ...");
      end if;
      DoblDobl_Pade_Trackers.Track_One_Path
        (standard_output,abh,jm,hs,fhm,fcf,ejm,mlt,ls.all,pars,mhom,idz,
         nbrsteps,nbrcorrs,cntcut,cntfail,minsize,maxsize,
         cntsstp,cntdstp,cntpstp,verbose,vrb-1);
      if verbose then
        Series_and_Trackers.Write_Path_Statistics
          (standard_output,nbrsteps,nbrcorrs,cntcut,cntfail,minsize,maxsize,
           cntdstp,cntpstp);
      end if;
      put("Solution "); put(i,1); put_line(" :"); put(ls.all); new_line;
      put("Continue to the next path ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Series_and_Trackers.Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Series_and_Trackers.Update_MinMax(smallest,largest,minsize,maxsize);
      Series_and_Trackers.Update_Ratio_Sums(ratsstp,ratdstp,ratpstp,
        cntsstp,cntdstp,cntpstp,nbrsteps*natural32(len));
    end loop;
    tstop(timer);
    Series_and_Trackers.Write_Total_Path_Statistics
      (standard_output,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,
       smallest,largest,ratdstp,ratpstp);
    new_line;
    put_line("THE SOLUTIONS :");
    put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    Write_Timer(standard_output,pars.numdeg,pars.dendeg,1,timer);
    if idxpar = 0
     then Refine_Roots(standard_output,abh,mhom,idz,sols,vrb);
    end if;
    Write_Conclusion(standard_output,start_moment);
    DoblDobl_Complex_Poly_Systems.Clear(h);
    DoblDobl_Complex_Poly_SysFun.Clear(abh);
    DoblDobl_CSeries_Poly_Systems.Clear(s);
    DoblDobl_CSeries_Poly_SysFun.Clear(fhm);
    DoblDobl_Complex_Series_VecVecs.Clear(fcf);
    DoblDobl_CSeries_Jaco_Matrices.Clear(ejm);
    DoblDobl_Cseries_Jaco_Matrices.Clear(mlt);
  end DoblDobl_Run;

  procedure QuadDobl_Run
              ( monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq);
    abh : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(1..nq)
        := Homotopy_Mixed_Residuals.QuadDobl_AbsVal_Homotopy;
    s : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq);
    jm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
    fhm : QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(1..nq);
    fcf : QuadDobl_Complex_Series_VecVecs.VecVec(1..nq);
    dim : constant integer32 := Set_Dimension(nvr,idxpar);
    ejm : QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(h'range,1..dim);
    mlt : QuadDobl_CSeries_Jaco_Matrices.Mult_Factors(h'range,1..dim);
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    ans : character;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs,cntcut,cntfail : natural32;
    minsize,maxsize,smallest,largest : double_float;
    cntsstp,cntdstp,cntpstp : natural32;
    ratsstp,ratdstp,ratpstp : double_float := 0.0;
    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

    use Singular_Values_of_Hessians;

  begin
    if vrb > 0
     then put_line("-> in series_path_trackers.QuadDobl_Run ...");
    end if;
    if idxpar /= 0
     then QuadDobl_Jacobian_Hessians_of_Homotopy(idxpar,jm,hs);
     else QuadDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
    end if;
    h := QuadDobl_Homotopy.Homotopy_System;
    if idxpar /= 0
     then s := Series_and_Homotopies.Create(h,idxpar,false);
     else s := Series_and_Homotopies.Create(h,nq+1,false);
    end if;
    fhm := QuadDobl_CSeries_Poly_SysFun.Create(s);
    fcf := QuadDobl_CSeries_Poly_SysFun.Coeff(s);
    QuadDobl_CSeries_Jaco_Matrices.Create(s,ejm,mlt);
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (pars.maxsteps+1)*pars.corsteps+1; maxnbrcorrs := 0;
    smallest := pars.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put("Tracking path "); put(i,1); put_line(" ...");
      end if;
      QuadDobl_Pade_Trackers.Track_One_Path
        (standard_output,abh,jm,hs,fhm,fcf,ejm,mlt,ls.all,pars,mhom,idz,
         nbrsteps,nbrcorrs,cntcut,cntfail,minsize,maxsize,
         cntsstp,cntdstp,cntpstp,verbose,vrb-1);
      if verbose then
        Series_and_Trackers.Write_Path_Statistics
          (standard_output,nbrsteps,nbrcorrs,cntcut,cntfail,minsize,maxsize,
           cntdstp,cntpstp);
      end if;
      put("Solution "); put(i,1); put_line(" :"); put(ls.all); new_line;
      put("Continue to the next path ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Series_and_Trackers.Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Series_and_Trackers.Update_MinMax(smallest,largest,minsize,maxsize);
      Series_and_Trackers.Update_Ratio_Sums(ratsstp,ratdstp,ratpstp,
        cntsstp,cntdstp,cntpstp,nbrsteps*natural32(len));
    end loop;
    tstop(timer);
    Series_and_Trackers.Write_Total_Path_Statistics
      (standard_output,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,
       smallest,largest,ratdstp,ratpstp);
    new_line;
    put_line("THE SOLUTIONS :");
    put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    Write_Timer(standard_output,pars.numdeg,pars.dendeg,2,timer);
    if idxpar = 0
     then Refine_Roots(standard_output,abh,mhom,idz,sols,vrb);
    end if;
    Write_Conclusion(standard_output,start_moment);
    QuadDobl_Complex_Poly_Systems.Clear(h);
    QuadDobl_Complex_Poly_SysFun.Clear(abh);
    QuadDobl_CSeries_Poly_Systems.Clear(s);
    QuadDobl_CSeries_Poly_SysFun.Clear(fhm);
    QuadDobl_Complex_Series_VecVecs.Clear(fcf);
    QuadDobl_CSeries_Jaco_Matrices.Clear(ejm);
    QuadDobl_Cseries_Jaco_Matrices.Clear(mlt);
  end QuadDobl_Run;

  procedure Standard_Run
              ( file : in file_type; monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out Standard_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq);
    abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(1..nq)
        := Homotopy_Mixed_Residuals.Standard_AbsVal_Homotopy;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq);
    fhm : Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(1..nq);
    fcf : Standard_Complex_Series_VecVecs.VecVec(1..nq);
    dim : constant integer32 := Set_Dimension(nvr,idxpar);
    jm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : Standard_Complex_Hessians.Link_to_Array_of_Hessians;
    ejm : Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(h'range,1..dim);
    mlt : Standard_CSeries_Jaco_Matrices.Mult_Factors(h'range,1..dim);
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs,cntcut,cntfail : natural32;
    minsize,maxsize,smallest,largest : double_float;
    cntsstp,cntdstp,cntpstp : natural32;
    ratsstp,ratdstp,ratpstp : double_float := 0.0;
    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

    use Singular_Values_of_Hessians;

  begin
    if vrb > 0
     then put_line("-> in series_path_trackers.Standard_Run ...");
    end if;
    if idxpar /= 0
     then Standard_Jacobian_Hessians_of_Homotopy(idxpar,jm,hs);
     else Standard_Jacobian_Hessians_of_Homotopy(jm,hs);
    end if;
    h := Standard_Homotopy.Homotopy_System;
    if idxpar /= 0
     then s := Series_and_Homotopies.Create(h,idxpar,false);
     else s := Series_and_Homotopies.Create(h,nq+1,false);
    end if;
    fhm := Standard_CSeries_Poly_SysFun.Create(s);
    fcf := Standard_CSeries_Poly_SysFun.Coeff(s);
    Standard_CSeries_Jaco_Matrices.Create(s,ejm,mlt);
    Standard_Write(file,natural32(nq),natural32(nvr),idxpar,sols,pars);
    if mhom > 1 then
      Write_Partition(file,natural32(nvr)-mhom,mhom,idz);
      new_line(file);
    end if;
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (pars.maxsteps+1)*pars.corsteps+1; maxnbrcorrs := 0;
    smallest := pars.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put("Tracking path "); put(i,1); put_line(" ...");
      end if;
      Standard_Pade_Trackers.Track_One_Path
        (file,abh,jm,hs,fhm,fcf,ejm,mlt,ls.all,pars,mhom,idz,nbrsteps,nbrcorrs,
         cntcut,cntfail,minsize,maxsize,cntsstp,cntdstp,cntpstp,verbose,vrb-1);
      if verbose then
        Series_and_Trackers.Write_Path_Statistics
          (file,nbrsteps,nbrcorrs,cntcut,cntfail,minsize,maxsize,
           cntdstp,cntpstp);
      end if;
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      put(file,ls.all); new_line(file);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Series_and_Trackers.Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Series_and_Trackers.Update_MinMax(smallest,largest,minsize,maxsize);
      Series_and_Trackers.Update_Ratio_Sums(ratsstp,ratdstp,ratpstp,
        cntsstp,cntdstp,cntpstp,nbrsteps*natural32(len));
    end loop;
    tstop(timer);
    Series_and_Trackers.Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,
       smallest,largest,ratdstp,ratpstp);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    Write_Timer(file,pars.numdeg,pars.dendeg,0,timer);
    if idxpar = 0
     then Refine_Roots(file,abh,mhom,idz,sols,vrb);
    end if;
    Write_Conclusion(file,start_moment);
    Standard_Complex_Poly_Systems.Clear(h);
    Standard_Complex_Poly_SysFun.Clear(abh);
    Standard_CSeries_Poly_Systems.Clear(s);
    Standard_CSeries_Poly_SysFun.Clear(fhm);
    Standard_Complex_Series_VecVecs.Clear(fcf);
    Standard_CSeries_Jaco_Matrices.Clear(ejm);
    Standard_Cseries_Jaco_Matrices.Clear(mlt);
  end Standard_Run;

  procedure DoblDobl_Run
              ( file : in file_type; monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq);
    abh : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(1..nq)
        := Homotopy_Mixed_Residuals.DoblDobl_AbsVal_Homotopy;
    s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq);
    jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
    fhm : DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(1..nq);
    fcf : DoblDobl_Complex_Series_VecVecs.VecVec(1..nq);
    dim : constant integer32 := Set_Dimension(nvr,idxpar);
    ejm : DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(h'range,1..dim);
    mlt : DoblDobl_CSeries_Jaco_Matrices.Mult_Factors(h'range,1..dim);
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs,cntcut,cntfail : natural32;
    minsize,maxsize,smallest,largest : double_float;
    cntsstp,cntdstp,cntpstp : natural32;
    ratsstp,ratdstp,ratpstp : double_float := 0.0;
    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

    use Singular_Values_of_Hessians;

  begin
    if vrb > 0
     then put_line("-> in series_path_trackers.DoblDobl_Run ...");
    end if;
    if idxpar /= 0
     then DoblDobl_Jacobian_Hessians_of_Homotopy(idxpar,jm,hs);
     else DoblDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
    end if;
    h := DoblDobl_Homotopy.Homotopy_System;
    if idxpar /= 0
     then s := Series_and_Homotopies.Create(h,idxpar,false);
     else s := Series_and_Homotopies.Create(h,nq+1,false);
    end if;
    fhm := DoblDobl_CSeries_Poly_SysFun.Create(s);
    fcf := DoblDobl_CSeries_Poly_SysFun.Coeff(s);
    DoblDobl_CSeries_Jaco_Matrices.Create(s,ejm,mlt);
    DoblDobl_Write(file,natural32(nq),natural32(nvr),idxpar,sols,pars);
    if mhom > 1 then
      Write_Partition(file,natural32(nvr)-mhom,mhom,idz);
      new_line(file);
    end if;
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (pars.maxsteps+1)*pars.corsteps+1; maxnbrcorrs := 0;
    smallest := pars.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put("Tracking path "); put(i,1); put_line(" ...");
      end if;
      DoblDobl_Pade_Trackers.Track_One_Path
        (file,abh,jm,hs,fhm,fcf,ejm,mlt,ls.all,pars,mhom,idz,nbrsteps,nbrcorrs,
         cntcut,cntfail,minsize,maxsize,cntsstp,cntdstp,cntpstp,verbose,vrb-1);
      if verbose then
        Series_and_Trackers.Write_Path_Statistics
          (file,nbrsteps,nbrcorrs,cntcut,cntfail,minsize,maxsize,
           cntdstp,cntpstp);
      end if;
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      put(file,ls.all); new_line(file);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Series_and_Trackers.Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Series_and_Trackers.Update_MinMax(smallest,largest,minsize,maxsize);
      Series_and_Trackers.Update_Ratio_Sums(ratsstp,ratdstp,ratpstp,
        cntsstp,cntdstp,cntpstp,nbrsteps*natural32(len));
    end loop;
    tstop(timer);
    Series_and_Trackers.Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,
       smallest,largest,ratdstp,ratpstp);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    Write_Timer(file,pars.numdeg,pars.dendeg,1,timer);
    if idxpar = 0
     then Refine_Roots(file,abh,mhom,idz,sols,vrb);
    end if;
    Write_Conclusion(file,start_moment);
    DoblDobl_Complex_Poly_Systems.Clear(h);
    DoblDobl_Complex_Poly_SysFun.Clear(abh);
    DoblDobl_CSeries_Poly_Systems.Clear(s);
    DoblDobl_CSeries_Poly_SysFun.Clear(fhm);
    DoblDobl_Complex_Series_VecVecs.Clear(fcf);
    DoblDobl_CSeries_Jaco_Matrices.Clear(ejm);
    DoblDobl_Cseries_Jaco_Matrices.Clear(mlt);
  end DoblDobl_Run;

  procedure QuadDobl_Run
              ( file : in file_type; monitor,verbose : in boolean;
                nq,nvr,idxpar : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                vrb : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq);
    abh : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(1..nq)
        := Homotopy_Mixed_Residuals.QuadDobl_AbsVal_Homotopy;
    s : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq);
    jm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
    fhm : QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(1..nq);
    fcf : QuadDobl_Complex_Series_VecVecs.VecVec(1..nq);
    dim : constant integer32 := Set_Dimension(nvr,idxpar);
    ejm : QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(h'range,1..dim);
    mlt : QuadDobl_CSeries_Jaco_Matrices.Mult_Factors(h'range,1..dim);
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs,cntcut,cntfail : natural32;
    minsize,maxsize,smallest,largest : double_float;
    cntsstp,cntdstp,cntpstp : natural32;
    ratsstp,ratdstp,ratpstp : double_float := 0.0;
    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

    use Singular_Values_of_Hessians;

  begin
    if vrb > 0
     then put_line("-> in series_path_trackers.QuadDobl_Run ...");
    end if;
    if idxpar /= 0
     then QuadDobl_Jacobian_Hessians_of_Homotopy(idxpar,jm,hs);
     else QuadDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
    end if;
    h := QuadDobl_Homotopy.Homotopy_System;
    if idxpar /= 0
     then s := Series_and_Homotopies.Create(h,idxpar,false);
     else s := Series_and_Homotopies.Create(h,nq+1,false);
    end if;
    fhm := QuadDobl_CSeries_Poly_SysFun.Create(s);
    fcf := QuadDobl_CSeries_Poly_SysFun.Coeff(s);
    QuadDobl_CSeries_Jaco_Matrices.Create(s,ejm,mlt);
    QuadDobl_Write(file,natural32(nq),natural32(nvr),idxpar,sols,pars);
    if mhom > 1 then
      Write_Partition(file,natural32(nvr)-mhom,mhom,idz);
      new_line(file);
    end if;
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (pars.maxsteps+1)*pars.corsteps+1; maxnbrcorrs := 0;
    smallest := pars.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put("Tracking path "); put(i,1); put_line(" ...");
      end if;
      QuadDobl_Pade_Trackers.Track_One_Path
        (file,abh,jm,hs,fhm,fcf,ejm,mlt,ls.all,pars,mhom,idz,nbrsteps,nbrcorrs,
         cntcut,cntfail,minsize,maxsize,cntsstp,cntdstp,cntpstp,verbose,vrb-1);
      if verbose then
        Series_and_Trackers.Write_Path_Statistics
          (file,nbrsteps,nbrcorrs,cntcut,cntfail,minsize,maxsize,
           cntdstp,cntpstp);
      end if;
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      put(file,ls.all); new_line(file);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Series_and_Trackers.Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Series_and_Trackers.Update_MinMax(smallest,largest,minsize,maxsize);
      Series_and_Trackers.Update_Ratio_Sums(ratsstp,ratdstp,ratpstp,
        cntsstp,cntdstp,cntpstp,nbrsteps*natural32(len));
    end loop;
    tstop(timer);
    Series_and_Trackers.Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,
       smallest,largest,ratdstp,ratpstp);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    Write_Timer(file,pars.numdeg,pars.dendeg,2,timer);
    if idxpar = 0
     then Refine_Roots(file,abh,mhom,idz,sols,vrb);
    end if;
    Write_Conclusion(file,start_moment);
    QuadDobl_Complex_Poly_Systems.Clear(h);
    QuadDobl_Complex_Poly_SysFun.Clear(abh);
    QuadDobl_CSeries_Poly_Systems.Clear(s);
    QuadDobl_CSeries_Poly_SysFun.Clear(fhm);
    QuadDobl_Complex_Series_VecVecs.Clear(fcf);
    QuadDobl_CSeries_Jaco_Matrices.Clear(ejm);
    QuadDobl_Cseries_Jaco_Matrices.Clear(mlt);
  end QuadDobl_Run;

  function Prompt_for_Artificial return boolean is

    ans : character;

  begin
    new_line;
    put_line("Either a homotopy has a parameter among its variables,");
    put_line("or the parameter is artificial and the homotopy connects");
    put_line("a target system to a start system with known solutions.");
    put("Is the homotopy an artificial parameter homotopy ? (y/n) ");
    Ask_Yes_or_No(ans);
    return (ans = 'y');
  end Prompt_for_Artificial;

  function Prompt_for_Homogenization
             ( dim : in natural32 ) return natural32 is

    res : natural32 := 0;

  begin
    new_line;
    put_line
      ("MENU for affine, homogeneous or multi-homogeneous coordinates :");
    put_line("  0 : in affine coordinates, in the original variables;");
    put_line("  1 : in 1-homogeous coordinates, in projective space;");
    put_line("  2 or higher : in multi-homogeous coordinates, in a multi-");
    put_line("  projective space defined by a partition of the variables.");
    loop
      put("Type a number between 0 and "); put(dim,1); put(" : ");
      Read_Natural(res);
      exit when res <= dim;
      put("-> your number is too high, as ");
      put(res,1); put(" > "); put(dim,1);
      put_line("; please try again.");
    end loop;
    return res;
  end Prompt_for_Homogenization;

  function Prompt_for_Partition
             ( nvr,mhom : in natural32 ) 
             return Standard_Natural_Vectors.Vector is
  begin
    new_line;
    put("Let us define a partition of "); put(nvr,1);
    put(" variables, of size "); put(mhom,1); put_line(" ...");
    declare
      res : constant Standard_Natural_Vectors.Vector(1..integer32(nvr))
          := Partitions_of_Sets_of_Unknowns_io.iget(mhom);
    begin
      return res;
    end;
  end Prompt_for_Partition;

  procedure Define_Partition
              ( n : in natural32; m : in out natural32;
                idx : out Standard_Natural_Vectors.Link_to_Vector;
                z : out Link_to_Partition ) is

    ans : character;

  begin
    loop
      declare
        wix : Standard_Natural_Vectors.Vector(1..integer32(n));
        wz : Partition(1..m);
      begin
        wix := Prompt_for_Partition(n,m);
        wz := Partitions_of_Sets_of_Unknowns_io.Make_Partition(n,m,wix);
        put("-> your partition : ");
        Partitions_of_Sets_of_Unknowns_io.put(wz);
        put(" Is this okay ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          idx := new Standard_Natural_Vectors.Vector'(wix);
          z := new Partition'(wz);
          return;
        end if;
      end;
      loop
        new_line;
        put("Give the number of sets in the partition : ");
        Read_Natural(m);
        exit when (m > 1 and m <= n);
        if m < 2 then
          put_line("-> Please enter a number larger than 1.");
        elsif m > n then
          put_line("-> Please enter a number not larger than ");
          put(n,1); put_line(".");
        end if;
      end loop;
    end loop;
  end Define_Partition;

  procedure Add_Multihomogeneous_Symbols
              ( m : in natural32; prefix : in string := "Z" ) is
  begin
    Symbol_Table.Enlarge(m);
    if m = 1 then
      Symbol_Table.Add_String(prefix & "0");
    else
      for i in 1..m loop
        declare
          sv : constant string := prefix & Characters_and_Numbers.nConvert(i);
        begin
          Symbol_Table.Add_String(sv);
        end;
      end loop;
    end if;
  end Add_Multihomogeneous_Symbols;

  procedure Write_Partition
              ( file : in file_type; n,m : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector ) is

    z : constant Partition(1..m)
      := Partitions_of_Sets_of_Unknowns_io.Make_Partition(n,m,idz.all);

  begin
    put(file,m,1);
    put_line(file,"-homogeneous coordinates are defined by the partition");
    Partitions_of_Sets_of_Unknowns_io.put(file,z);
    put(file," of the set of "); put(file,n,1); put_line(file," variables.");
  end Write_Partition;

  procedure Standard_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

    use Homotopy_Series_Readers;

  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
    nvr := Standard_Complex_Solutions.Head_Of(sols).n;
    nbq := target'last;
    mhom := Prompt_for_Homogenization(natural32(nvr));
    if mhom = 0 then
      Standard_Homotopy.Create(target.all,start.all,2,gamma);
    else
      if mhom = 1 then
        Standard_Projective_Transformation(target,start,sols);
        Add_Multihomogeneous_Symbols(1);
        nvr := nvr + 1; nbq := nbq + 1;
      else
        Define_Partition(natural32(nvr),mhom,idz,z);
        Standard_Multi_Projective_Transformation(target,start,sols,mhom,z.all);
        Add_Multihomogeneous_Symbols(mhom);
        nvr := nvr + integer32(mhom); nbq := nbq + integer32(mhom);
      end if;
      Standard_Homotopy.Create(target.all,start.all,1,gamma);
      Standard_Coefficient_Homotopy.Create(start.all,target.all,1,gamma);
    end if;
  end Standard_Define_Homotopy;

  procedure DoblDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    target,start : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    dd_gamma : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(gamma);

    use Homotopy_Series_Readers;

  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    nvr := DoblDobl_Complex_Solutions.Head_Of(sols).n;
    nbq := target'last;
    mhom := Prompt_for_Homogenization(natural32(nvr));
    if mhom = 0 then
      DoblDobl_Homotopy.Create(target.all,start.all,2,dd_gamma);
    else
      if mhom = 1 then
        DoblDobl_Projective_Transformation(target,start,sols);
        Add_Multihomogeneous_Symbols(1);
        nvr := nvr + 1; nbq := nbq + 1;
      else
        Define_Partition(natural32(nvr),mhom,idz,z);
        DoblDobl_Multi_Projective_Transformation(target,start,sols,mhom,z.all);
        Add_Multihomogeneous_Symbols(mhom);
        nvr := nvr + integer32(mhom); nbq := nbq + integer32(mhom);
      end if;
      DoblDobl_Homotopy.Create(target.all,start.all,1,dd_gamma);
      DoblDobl_Coefficient_Homotopy.Create(start.all,target.all,1,dd_gamma);
    end if;
  end DoblDobl_Define_Homotopy;

  procedure TripDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out TripDobl_Complex_Solutions.Solution_List ) is

    target,start : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    td_gamma : constant TripDobl_Complex_Numbers.Complex_Number
             := Standard_to_TripDobl_Complex(gamma);

    use Homotopy_Series_Readers;

  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    TripDobl_System_and_Solutions_io.get(start,sols);
    nvr := TripDobl_Complex_Solutions.Head_Of(sols).n;
    nbq := target'last;
    mhom := Prompt_for_Homogenization(natural32(nvr));
    if mhom = 0 then
      TripDobl_Homotopy.Create(target.all,start.all,2,td_gamma);
    else
      if mhom = 1 then
        TripDobl_Projective_Transformation(target,start,sols);
        Add_Multihomogeneous_Symbols(1);
        nvr := nvr + 1; nbq := nbq + 1;
      else
        Define_Partition(natural32(nvr),mhom,idz,z);
        TripDobl_Multi_Projective_Transformation(target,start,sols,mhom,z.all);
        Add_Multihomogeneous_Symbols(mhom);
        nvr := nvr + integer32(mhom); nbq := nbq + integer32(mhom);
      end if;
      TripDobl_Homotopy.Create(target.all,start.all,1,td_gamma);
      TripDobl_Coefficient_Homotopy.Create(start.all,target.all,1,td_gamma);
    end if;
  end TripDobl_Define_Homotopy;

  procedure QuadDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    target,start : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    qd_gamma : constant QuadDobl_Complex_Numbers.Complex_Number
             := Standard_to_QuadDobl_Complex(gamma);

    use Homotopy_Series_Readers;

  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    nvr := QuadDobl_Complex_Solutions.Head_Of(sols).n;
    nbq := target'last;
    mhom := Prompt_for_Homogenization(natural32(nvr));
    if mhom = 0 then
      QuadDobl_Homotopy.Create(target.all,start.all,2,qd_gamma);
    else
      if mhom = 1 then
        QuadDobl_Projective_Transformation(target,start,sols);
        Add_Multihomogeneous_Symbols(1);
        nvr := nvr + 1; nbq := nbq + 1;
      else
        Define_Partition(natural32(nvr),mhom,idz,z);
        QuadDobl_Multi_Projective_Transformation(target,start,sols,mhom,z.all);
        Add_Multihomogeneous_Symbols(mhom);
        nvr := nvr + integer32(mhom); nbq := nbq + integer32(mhom);
      end if;
      QuadDobl_Homotopy.Create(target.all,start.all,1,qd_gamma);
      QuadDobl_Coefficient_Homotopy.Create(start.all,target.all,1,qd_gamma);
    end if;
  end QuadDobl_Define_Homotopy;

  procedure PentDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out PentDobl_Complex_Solutions.Solution_List ) is

    target,start : PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    pd_gamma : constant PentDobl_Complex_Numbers.Complex_Number
             := Standard_to_PentDobl_Complex(gamma);

    use Homotopy_Series_Readers;

  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    PentDobl_System_and_Solutions_io.get(start,sols);
    nvr := PentDobl_Complex_Solutions.Head_Of(sols).n;
    nbq := target'last;
    mhom := Prompt_for_Homogenization(natural32(nvr));
    if mhom = 0 then
      PentDobl_Homotopy.Create(target.all,start.all,2,pd_gamma);
    else
      if mhom = 1 then
        PentDobl_Projective_Transformation(target,start,sols);
        Add_Multihomogeneous_Symbols(1);
        nvr := nvr + 1; nbq := nbq + 1;
      else
        Define_Partition(natural32(nvr),mhom,idz,z);
        PentDobl_Multi_Projective_Transformation(target,start,sols,mhom,z.all);
        Add_Multihomogeneous_Symbols(mhom);
        nvr := nvr + integer32(mhom); nbq := nbq + integer32(mhom);
      end if;
      PentDobl_Homotopy.Create(target.all,start.all,1,pd_gamma);
      PentDobl_Coefficient_Homotopy.Create(start.all,target.all,1,pd_gamma);
    end if;
  end PentDobl_Define_Homotopy;

  procedure OctoDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out OctoDobl_Complex_Solutions.Solution_List ) is

    target,start : OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    od_gamma : constant OctoDobl_Complex_Numbers.Complex_Number
             := Standard_to_OctoDobl_Complex(gamma);

    use Homotopy_Series_Readers;

  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    OctoDobl_System_and_Solutions_io.get(start,sols);
    nvr := OctoDobl_Complex_Solutions.Head_Of(sols).n;
    nbq := target'last;
    mhom := Prompt_for_Homogenization(natural32(nvr));
    if mhom = 0 then
      OctoDobl_Homotopy.Create(target.all,start.all,2,od_gamma);
    else
      if mhom = 1 then
        OctoDobl_Projective_Transformation(target,start,sols);
        Add_Multihomogeneous_Symbols(1);
        nvr := nvr + 1; nbq := nbq + 1;
      else
        Define_Partition(natural32(nvr),mhom,idz,z);
        OctoDobl_Multi_Projective_Transformation(target,start,sols,mhom,z.all);
        Add_Multihomogeneous_Symbols(mhom);
        nvr := nvr + integer32(mhom); nbq := nbq + integer32(mhom);
      end if;
      OctoDobl_Homotopy.Create(target.all,start.all,1,od_gamma);
      OctoDobl_Coefficient_Homotopy.Create(start.all,target.all,1,od_gamma);
    end if;
  end OctoDobl_Define_Homotopy;

  procedure DecaDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out DecaDobl_Complex_Solutions.Solution_List ) is

    target,start : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    da_gamma : constant DecaDobl_Complex_Numbers.Complex_Number
             := Standard_to_DecaDobl_Complex(gamma);

    use Homotopy_Series_Readers;

  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    DecaDobl_System_and_Solutions_io.get(start,sols);
    nvr := DecaDobl_Complex_Solutions.Head_Of(sols).n;
    nbq := target'last;
    mhom := Prompt_for_Homogenization(natural32(nvr));
    if mhom = 0 then
      DecaDobl_Homotopy.Create(target.all,start.all,2,da_gamma);
    else
      if mhom = 1 then
        DecaDobl_Projective_Transformation(target,start,sols);
        Add_Multihomogeneous_Symbols(1);
        nvr := nvr + 1; nbq := nbq + 1;
      else
        Define_Partition(natural32(nvr),mhom,idz,z);
        DecaDobl_Multi_Projective_Transformation(target,start,sols,mhom,z.all);
        Add_Multihomogeneous_Symbols(mhom);
        nvr := nvr + integer32(mhom); nbq := nbq + integer32(mhom);
      end if;
      DecaDobl_Homotopy.Create(target.all,start.all,1,da_gamma);
      DecaDobl_Coefficient_Homotopy.Create(start.all,target.all,1,da_gamma);
    end if;
  end DecaDobl_Define_Homotopy;

  procedure HexaDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out HexaDobl_Complex_Solutions.Solution_List ) is

    target,start : HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    da_gamma : constant HexaDobl_Complex_Numbers.Complex_Number
             := Standard_to_HexaDobl_Complex(gamma);

    use Homotopy_Series_Readers;

  begin
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    HexaDobl_System_and_Solutions_io.get(start,sols);
    nvr := HexaDobl_Complex_Solutions.Head_Of(sols).n;
    nbq := target'last;
    mhom := Prompt_for_Homogenization(natural32(nvr));
    if mhom = 0 then
      HexaDobl_Homotopy.Create(target.all,start.all,2,da_gamma);
    else
      if mhom = 1 then
        HexaDobl_Projective_Transformation(target,start,sols);
        Add_Multihomogeneous_Symbols(1);
        nvr := nvr + 1; nbq := nbq + 1;
      else
        Define_Partition(natural32(nvr),mhom,idz,z);
        HexaDobl_Multi_Projective_Transformation(target,start,sols,mhom,z.all);
        Add_Multihomogeneous_Symbols(mhom);
        nvr := nvr + integer32(mhom); nbq := nbq + integer32(mhom);
      end if;
      HexaDobl_Homotopy.Create(target.all,start.all,1,da_gamma);
      HexaDobl_Coefficient_Homotopy.Create(start.all,target.all,1,da_gamma);
    end if;
  end HexaDobl_Define_Homotopy;

  procedure Standard_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out Standard_Complex_Solutions.Solution_List ) is
  begin
    Standard_Define_Homotopy(nbq,nvr,pars.gamma,mhom,z,idz,sols);
  end Standard_Define_Homotopy;

  procedure DoblDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is
  begin
    DoblDobl_Define_Homotopy(nbq,nvr,pars.gamma,mhom,z,idz,sols);
  end DoblDobl_Define_Homotopy;

  procedure TripDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out TripDobl_Complex_Solutions.Solution_List ) is
  begin
    TripDobl_Define_Homotopy(nbq,nvr,pars.gamma,mhom,z,idz,sols);
  end TripDobl_Define_Homotopy;

  procedure QuadDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is
  begin
    QuadDobl_Define_Homotopy(nbq,nvr,pars.gamma,mhom,z,idz,sols);
  end QuadDobl_Define_Homotopy;

  procedure PentDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out PentDobl_Complex_Solutions.Solution_List ) is
  begin
    PentDobl_Define_Homotopy(nbq,nvr,pars.gamma,mhom,z,idz,sols);
  end PentDobl_Define_Homotopy;

  procedure OctoDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out OctoDobl_Complex_Solutions.Solution_List ) is
  begin
    OctoDobl_Define_Homotopy(nbq,nvr,pars.gamma,mhom,z,idz,sols);
  end OctoDobl_Define_Homotopy;

  procedure DecaDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out DecaDobl_Complex_Solutions.Solution_List ) is
  begin
    DecaDobl_Define_Homotopy(nbq,nvr,pars.gamma,mhom,z,idz,sols);
  end DecaDobl_Define_Homotopy;

  procedure HexaDobl_Define_Homotopy
              ( nbq,nvr : out integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : out natural32; z : out Link_to_Partition;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                sols : out HexaDobl_Complex_Solutions.Solution_List ) is
  begin
    HexaDobl_Define_Homotopy(nbq,nvr,pars.gamma,mhom,z,idz,sols);
  end HexaDobl_Define_Homotopy;

  procedure Standard_Main ( vrb : in integer32 := 0 ) is

    nbq,nvr,idx : integer32;
    sols,dropsols : Standard_Complex_Solutions.Solution_List;
    arth : constant boolean := Prompt_for_Artificial;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    mhom : natural32 := 0; -- by default, in affine coordinates
    z : Link_to_Partition;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    monitor,verbose,tofile : boolean;
    file : file_type;

  begin
    if vrb > 0
     then put_line("-> in series_path_trackers.Standard_Main ...");
    end if;
    Set_Output(file,monitor,verbose,tofile);
    new_line;
    if arth then
      Homotopy_Continuation_Parameters_io.Tune(pars);
      Standard_Define_Homotopy(nbq,nvr,pars,mhom,z,idz,sols);
      idx := 0;
      if tofile then
        Standard_Run(file,monitor,verbose,nbq,nvr,idx,pars,
                     mhom,idz,sols,vrb-1);
      else
        Standard_Run(monitor,verbose,nbq,nvr,idx,pars,mhom,idz,sols,vrb-1);
      end if;
    else
      pars.gamma := Standard_Complex_Numbers.Create(1.0);
      Homotopy_Continuation_Parameters_io.Tune(pars);
      Homotopy_Series_Readers.Standard_Parameter_Reader(nbq,nvr,idx,sols);
      dropsols := Solution_Drops.Drop(sols,natural32(idx));
      if tofile then
        Standard_Run(file,monitor,verbose,nbq,nvr,idx,pars,
                     mhom,idz,dropsols,vrb-1);
      else
        Standard_Run(monitor,verbose,nbq,nvr,idx,pars,mhom,idz,dropsols,vrb-1);
      end if;
    end if;
  end Standard_Main;

  procedure DoblDobl_Main ( vrb : in integer32 := 0 ) is

    nbq,nvr,idx : integer32;
    sols,dropsols : DoblDobl_Complex_Solutions.Solution_List;
    arth : constant boolean := Prompt_for_Artificial;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    mhom : natural32 := 0; -- by default, in affine coordinates
    z : Link_to_Partition;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    monitor,verbose,tofile : boolean;
    file : file_type;

  begin
    if vrb > 0
     then put_line("-> in series_path_trackers.DoblDobl_Main ...");
    end if;
    Set_Output(file,monitor,verbose,tofile);
    new_line;
    if arth then
      Homotopy_Continuation_Parameters_io.Tune(pars);
      DoblDobl_Define_Homotopy(nbq,nvr,pars,mhom,z,idz,sols);
      idx := 0;
      if tofile then
        DoblDobl_Run(file,monitor,verbose,nbq,nvr,idx,pars,
                     mhom,idz,sols,vrb-1);
      else
        DoblDobl_Run(monitor,verbose,nbq,nvr,idx,pars,mhom,idz,sols,vrb-1);
      end if;
    else
      pars.gamma := Standard_Complex_Numbers.Create(1.0);
      Homotopy_Continuation_Parameters_io.Tune(pars);
      Homotopy_Series_Readers.DoblDobl_Parameter_Reader(nbq,nvr,idx,sols);
      dropsols := Solution_Drops.Drop(sols,natural32(idx));
      if tofile then
        DoblDobl_Run(file,monitor,verbose,nbq,nvr,idx,pars,
                     mhom,idz,dropsols,vrb-1);
      else
        DoblDobl_Run(monitor,verbose,nbq,nvr,idx,pars,mhom,idz,dropsols,vrb-1);
      end if;
    end if;
  end DoblDobl_Main;

  procedure QuadDobl_Main ( vrb : in integer32 := 0 ) is

    nbq,nvr,idx : integer32;
    sols,dropsols : QuadDobl_Complex_Solutions.Solution_List;
    arth : constant boolean := Prompt_for_Artificial;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    mhom : natural32 := 0; -- by default, in affine coordinates
    z : Link_to_Partition;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    monitor,verbose,tofile : boolean;
    file : file_type;

  begin
    if vrb > 0
     then put_line("-> in series_path_trackers.QuadDobl_Main ...");
    end if;
    Set_Output(file,monitor,verbose,tofile);
    new_line;
    if arth then
      Homotopy_Continuation_Parameters_io.Tune(pars);
      QuadDobl_Define_Homotopy(nbq,nvr,pars,mhom,z,idz,sols);
      idx := 0;
      if tofile then
        QuadDobl_Run(file,monitor,verbose,nbq,nvr,idx,pars,
                     mhom,idz,sols,vrb-1);
      else
        QuadDobl_Run(monitor,verbose,nbq,nvr,idx,pars,mhom,idz,sols,vrb-1);
      end if;
    else
      pars.gamma := Standard_Complex_Numbers.Create(1.0);
      Homotopy_Continuation_Parameters_io.Tune(pars);
      Homotopy_Series_Readers.QuadDobl_Parameter_Reader(nbq,nvr,idx,sols);
      dropsols := Solution_Drops.Drop(sols,natural32(idx));
      if tofile then
        QuadDobl_Run(file,monitor,verbose,nbq,nvr,idx,pars,
                     mhom,idz,dropsols,vrb-1);
      else
        QuadDobl_Run(monitor,verbose,nbq,nvr,idx,pars,mhom,idz,dropsols,vrb-1);
      end if;
    end if;
  end QuadDobl_Main;

end Series_Path_Trackers;
