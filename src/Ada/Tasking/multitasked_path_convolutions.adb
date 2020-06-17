with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;                        use Time_Stamps;
with Write_Seed_Number;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;    use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;    use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;    use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;       use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;       use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;       use QuadDobl_Complex_Solutions_io;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;
with Homotopy_Continuation_Parameters_io;
with Residual_Convolution_Circuits;
with Corrector_Convolutions;             use Corrector_Convolutions;
with Predictor_Corrector_Trackers;       use Predictor_Corrector_Trackers;
with Projective_Transformations;
with Multi_Projective_Transformations;
with Track_Path_Convolutions;
with Standard_Solutions_Queue;
with DoblDobl_Solutions_Queue;
with QuadDobl_Solutions_Queue;
with Multitasking;
with Greeting_Banners;

package body Multitasked_Path_Convolutions is

  procedure Allocate ( v : in out Standard_Integer_VecVecs.VecVec;
                       n : in integer32 ) is
  begin
    for k in v'range loop
      v(k) := new Standard_Integer_Vectors.Vector(1..n);
    end loop;
  end Allocate;

  procedure Allocate ( v : in out Standard_Complex_VecVecs.VecVec;
                       n : in integer32 ) is
  begin
    for k in v'range loop
      v(k) := new Standard_Complex_Vectors.Vector(1..n);
    end loop;
  end Allocate;

  procedure Allocate ( v : in out DoblDobl_Complex_VecVecs.VecVec;
                       n : in integer32 ) is
  begin
    for k in v'range loop
      v(k) := new DoblDobl_Complex_Vectors.Vector(1..n);
    end loop;
  end Allocate;

  procedure Allocate ( v : in out QuadDobl_Complex_VecVecs.VecVec;
                       n : in integer32 ) is
  begin
    for k in v'range loop
      v(k) := new QuadDobl_Complex_Vectors.Vector(1..n);
    end loop;
  end Allocate;

  procedure Standard_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true ) is

    use Standard_Complex_Solutions,Standard_Speelpenning_Convolutions;
    use Standard_Predictor_Convolutions;

    maxit : constant integer32 := integer32(pars.numdeg + pars.dendeg + 2)/2;
    homsa,abhsa : System_Array(1..nbtasks);
    homlead,abhlead : VecVecVec(1..nbtasks);
    homcff : VecVecVec_Array(1..nbtasks);
    dx : Standard_Complex_VecVecs.VecVec(1..nbtasks);
    wrk : Standard_Complex_VecVecs.VecVec(1..nbtasks);
    ipvt : Standard_Integer_VecVecs.VecVec(1..nbtasks);
    first : constant Link_to_Solution := Head_Of(sols);
    prd : Predictor_Array(1..nbtasks)
        := Create(nbtasks,first.v,hom.neq,hom.deg,
                  integer32(pars.numdeg),integer32(pars.dendeg),SVD);
    psv : Predictor_Vectors_Array(1..nbtasks)
        := Create(nbtasks,hom.dim,hom.neq);
    svh : SVD_Hessians_Array(1..nbtasks) := Create(nbtasks,hom.dim);
    xtr : constant integer32 := 1;

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      t,mixres,minstpz,maxstpz : double_float := 0.0;
      tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
      nbrit : integer32 := 0;
      fail : boolean;

    begin
      loop
        myptr := Standard_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := 0.0;
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,maxit,
          mhom,idz,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,wrk(i),
          t,mixres,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz,fail);
        Restore_Coefficients(homcff(i),homsa(i).crc);
        Residual_Convolution_Circuits.Update_Radii_of_Constants
          (abhsa(i),homsa(i));
        Step_Coefficient(homsa(i),t);
        LU_Newton_Steps(homsa(i),abhsa(i),psv(i).all,integer32(pars.corsteps),
          nbrit,pars.tolres,ls.err,ls.res,dx(i).all,ipvt(i).all,ls.rco,
          fail,xtr);
        ls.v := psv(i).sol; ls.t := Standard_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
        Restore_Leading_Coefficients(abhlead(i),abhsa(i).crc);
        Restore_Coefficients(homcff(i),homsa(i).crc);
      end loop;
    end Silent_Track;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Track);

    procedure Report_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks with intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      cnt : integer32;
      t,mixres,minstpz,maxstpz : double_float := 0.0;
      tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
      nbrit : integer32 := 0;
      fail : boolean;

    begin
      loop
        myptr := Standard_Solutions_Queue.Next;
        cnt := Standard_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        put_line("Task " & Multitasking.to_string(i)
                         & " tracks path "
                         & Multitasking.to_string(cnt));
        psv(i).sol := ls.v; t := 0.0;
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,maxit,
          mhom,idz,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,wrk(i),
          t,mixres,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz,fail);
        Restore_Coefficients(homcff(i),homsa(i).crc);
        Residual_Convolution_Circuits.Update_Radii_of_Constants
          (abhsa(i),homsa(i));
        Step_Coefficient(homsa(i),t);
        LU_Newton_Steps(homsa(i),abhsa(i),psv(i).all,integer32(pars.corsteps),
          nbrit,pars.tolres,ls.err,ls.res,dx(i).all,ipvt(i).all,ls.rco,
          fail,xtr);
        ls.v := psv(i).sol; ls.t := Standard_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
        Restore_Leading_Coefficients(abhlead(i),abhsa(i).crc);
        Restore_Coefficients(homcff(i),homsa(i).crc);
      end loop;
    end Report_Track;
    procedure report_do_jobs is
      new Multitasking.Reporting_Workers(Report_Track);

  begin
    for k in 1..nbtasks loop
      Copy(hom,homsa(k)); Copy(abh,abhsa(k));
      Allocate_Leading_Coefficients(hom.crc,homlead(k));
      Allocate_Leading_Coefficients(abh.crc,abhlead(k));
      Store_Leading_Coefficients(hom.crc,homlead(k));
      Store_Leading_Coefficients(abh.crc,abhlead(k));
      Allocate_Coefficients(hom.crc,homcff(k));
      Store_Coefficients(hom.crc,homcff(k));
    end loop;
    Allocate(ipvt,hom.dim); Allocate(dx,hom.dim);
    wrk := Allocate_Coefficients(nbtasks,hom.deg);
    Standard_Solutions_Queue.Initialize(sols);
    if verbose
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
    Clear(homsa); Clear(abhsa);
    Clear(homlead); Clear(abhlead); Clear(homcff);
    Standard_Complex_VecVecs.Clear(dx);
    Standard_Complex_VecVecs.Clear(wrk);
    Standard_Integer_VecVecs.Clear(ipvt);
    Clear(prd); Clear(psv); Clear(svh);
  end Standard_Multitasked_Tracker;

  procedure DoblDobl_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true ) is

    use DoblDobl_Complex_Solutions,DoblDobl_Speelpenning_Convolutions;
    use DoblDobl_Predictor_Convolutions;

    maxit : constant integer32 := integer32(pars.numdeg + pars.dendeg + 2)/2;
    homsa,abhsa : System_Array(1..nbtasks);
    homlead,abhlead : VecVecVec(1..nbtasks);
    homcff : VecVecVec_Array(1..nbtasks);
    dx : DoblDobl_Complex_VecVecs.VecVec(1..nbtasks);
    wrk : DoblDobl_Complex_VecVecs.VecVec(1..nbtasks);
    ipvt : Standard_Integer_VecVecs.VecVec(1..nbtasks);
    first : constant Link_to_Solution := Head_Of(sols);
    prd : Predictor_Array(1..nbtasks)
        := Create(nbtasks,first.v,hom.neq,hom.deg,
                  integer32(pars.numdeg),integer32(pars.dendeg),SVD);
    psv : Predictor_Vectors_Array(1..nbtasks)
        := Create(nbtasks,hom.dim,hom.neq);
    svh : SVD_Hessians_Array(1..nbtasks) := Create(nbtasks,hom.dim);
    xtr : constant integer32 := 1;

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      t,mixres : double_double;
      tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
      minstpz,maxstpz : double_float;
      nbrit : integer32 := 0;
      fail : boolean;

    begin
      loop
        myptr := DoblDobl_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := Create(0.0);
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,maxit,
          mhom,idz,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,wrk(i),
          t,mixres,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz,fail);
        Restore_Coefficients(homcff(i),homsa(i).crc);
        Residual_Convolution_Circuits.Update_Radii_of_Constants
          (abhsa(i),homsa(i));
        Step_Coefficient(homsa(i),t);
        LU_Newton_Steps(homsa(i),abhsa(i),psv(i).all,integer32(pars.corsteps),
          nbrit,pars.tolres,ls.err,ls.res,dx(i).all,ipvt(i).all,ls.rco,
          fail,xtr);
        ls.v := psv(i).sol; ls.t := DoblDobl_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
        Restore_Leading_Coefficients(abhlead(i),abhsa(i).crc);
        Restore_Coefficients(homcff(i),homsa(i).crc);
      end loop;
    end Silent_Track;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Track);

    procedure Report_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks with intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      cnt : integer32;
      t,mixres : double_double;
      tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
      minstpz,maxstpz : double_float;
      nbrit : integer32 := 0;
      fail : boolean;

    begin
      loop
        myptr := DoblDobl_Solutions_Queue.Next;
        cnt := DoblDobl_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        put_line("Task " & Multitasking.to_string(i)
                         & " tracks path "
                         & Multitasking.to_string(cnt));
        psv(i).sol := ls.v; t := Create(0.0);
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,maxit,
          mhom,idz,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,wrk(i),
          t,mixres,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz,fail);
        Restore_Coefficients(homcff(i),homsa(i).crc);
        Residual_Convolution_Circuits.Update_Radii_of_Constants
          (abhsa(i),homsa(i));
        Step_Coefficient(homsa(i),t);
        LU_Newton_Steps(homsa(i),abhsa(i),psv(i).all,integer32(pars.corsteps),
          nbrit,pars.tolres,ls.err,ls.res,dx(i).all,ipvt(i).all,ls.rco,
          fail,xtr);
        ls.v := psv(i).sol; ls.t := DoblDobl_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
        Restore_Leading_Coefficients(abhlead(i),abhsa(i).crc);
        Restore_Coefficients(homcff(i),homsa(i).crc);
      end loop;
    end Report_Track;
    procedure report_do_jobs is
      new Multitasking.Reporting_Workers(Report_Track);

  begin
    for k in homsa'range loop
      Copy(hom,homsa(k)); Copy(abh,abhsa(k));
      Allocate_Leading_Coefficients(hom.crc,homlead(k));
      Allocate_Leading_Coefficients(abh.crc,abhlead(k));
      Store_Leading_Coefficients(hom.crc,homlead(k));
      Store_Leading_Coefficients(abh.crc,abhlead(k));
      Allocate_Coefficients(hom.crc,homcff(k));
      Store_Coefficients(hom.crc,homcff(k));
    end loop;
    Allocate(ipvt,hom.dim); Allocate(dx,hom.dim);
    wrk := Allocate_Coefficients(nbtasks,hom.deg);
    DoblDobl_Solutions_Queue.Initialize(sols);
    if verbose
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
    Clear(homsa); Clear(abhsa);
    Clear(homlead); Clear(abhlead); Clear(homcff);
    DoblDobl_Complex_VecVecs.Clear(dx);
    DoblDobl_Complex_VecVecs.Clear(wrk);
    Standard_Integer_VecVecs.Clear(ipvt);
    Clear(prd); Clear(psv); Clear(svh);
  end DoblDobl_Multitasked_Tracker;

  procedure QuadDobl_Multitasked_Tracker
              ( nbtasks : in integer32;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true ) is

    use QuadDobl_Complex_Solutions,QuadDobl_Speelpenning_Convolutions;
    use QuadDobl_Predictor_Convolutions;

    maxit : constant integer32 := integer32(pars.numdeg + pars.dendeg + 2)/2;
    homsa,abhsa : System_Array(1..nbtasks);
    homlead,abhlead : VecVecVec(1..nbtasks);
    homcff : VecVecVec_Array(1..nbtasks);
    dx : QuadDobl_Complex_VecVecs.VecVec(1..nbtasks);
    wrk : QuadDobl_Complex_VecVecs.VecVec(1..nbtasks);
    ipvt : Standard_Integer_VecVecs.VecVec(1..nbtasks);
    first : constant Link_to_Solution := Head_Of(sols);
    prd : Predictor_Array(1..nbtasks)
        := Create(nbtasks,first.v,hom.neq,hom.deg,
                  integer32(pars.numdeg),integer32(pars.dendeg),SVD);
    psv : Predictor_Vectors_Array(1..nbtasks)
        := Create(nbtasks,hom.dim,hom.neq);
    svh : SVD_Hessians_Array(1..nbtasks) := Create(nbtasks,hom.dim);
    xtr : constant integer32 := 1;

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      t,mixres : quad_double;
      tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
      minstpz,maxstpz : double_float;
      nbrit : integer32 := 0;
      fail : boolean;

    begin
      loop
        myptr := QuadDobl_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := Create(0.0);
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,maxit,
          mhom,idz,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,wrk(i),
          t,mixres,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz,fail);
        Restore_Coefficients(homcff(i),homsa(i).crc);
        Residual_Convolution_Circuits.Update_Radii_of_Constants
          (abhsa(i),homsa(i));
        Step_Coefficient(homsa(i),t);
        LU_Newton_Steps(homsa(i),abhsa(i),psv(i).all,integer32(pars.corsteps),
          nbrit,pars.tolres,ls.err,ls.res,dx(i).all,ipvt(i).all,ls.rco,
          fail,xtr);
        ls.v := psv(i).sol; ls.t := QuadDobl_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
        Restore_Leading_Coefficients(abhlead(i),abhsa(i).crc);
        Restore_Coefficients(homcff(i),homsa(i).crc);
      end loop;
    end Silent_Track;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Track);

    procedure Report_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks with intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      cnt : integer32;
      t,mixres : quad_double;
      tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
      minstpz,maxstpz : double_float;
      nbrit : integer32 := 0;
      fail : boolean;

    begin
      loop
        myptr := QuadDobl_Solutions_Queue.Next;
        cnt := QuadDobl_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        put_line("Task " & Multitasking.to_string(i)
                         & " tracks path "
                         & Multitasking.to_string(cnt));
        psv(i).sol := ls.v; t := Create(0.0);
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,maxit,
          mhom,idz,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,wrk(i),
          t,mixres,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz,fail);
        Restore_Coefficients(homcff(i),homsa(i).crc);
        Residual_Convolution_Circuits.Update_Radii_of_Constants
          (abhsa(i),homsa(i));
        Step_Coefficient(homsa(i),t);
        LU_Newton_Steps(homsa(i),abhsa(i),psv(i).all,integer32(pars.corsteps),
          nbrit,pars.tolres,ls.err,ls.res,dx(i).all,ipvt(i).all,ls.rco,
          fail,xtr);
        ls.v := psv(i).sol; ls.t := QuadDobl_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
        Restore_Leading_Coefficients(abhlead(i),abhsa(i).crc);
        Restore_Coefficients(homcff(i),homsa(i).crc);
      end loop;
    end Report_Track;
    procedure report_do_jobs is
      new Multitasking.Reporting_Workers(Report_Track);

  begin
    for k in homsa'range loop
      Copy(hom,homsa(k)); Copy(abh,abhsa(k));
      Allocate_Leading_Coefficients(hom.crc,homlead(k));
      Allocate_Leading_Coefficients(abh.crc,abhlead(k));
      Store_Leading_Coefficients(hom.crc,homlead(k));
      Store_Leading_Coefficients(abh.crc,abhlead(k));
      Allocate_Coefficients(hom.crc,homcff(k));
      Store_Coefficients(hom.crc,homcff(k));
    end loop;
    Allocate(ipvt,hom.dim); Allocate(dx,hom.dim);
    wrk := Allocate_Coefficients(nbtasks,hom.deg);
    QuadDobl_Solutions_Queue.Initialize(sols);
    if verbose
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
    Clear(homsa); Clear(abhsa);
    Clear(homlead); Clear(abhlead); Clear(homcff);
    QuadDobl_Complex_VecVecs.Clear(dx);
    QuadDobl_Complex_VecVecs.Clear(wrk);
    Standard_Integer_VecVecs.Clear(ipvt);
    Clear(prd); Clear(psv); Clear(svh);
  end QuadDobl_Multitasked_Tracker;

  procedure Track
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbt,mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean ) is

    ans : character;
    hcrd : constant boolean := (mhom > 0);
    verbose : boolean;
    startmoment,stopmoment : Ada.Calendar.Time;
  
  begin
    if not arth then
      put(file,natural32(hom.neq),natural32(hom.neq+1),
               Standard_Homotopy.Homotopy_System);
    else
      declare
        p : constant Standard_Complex_Poly_Systems.Poly_Sys
          := Standard_Homotopy.Target_System;
        q : constant Standard_Complex_Poly_Systems.Poly_Sys
          := Standard_Homotopy.Start_System;
      begin
        put(file,p'last,1); new_line(file); put(file,p);
        new_line(file);
        put_line(file,"THE START SYSTEM :");
        put(file,q'last,1); new_line(file); put(file,q);
      end;
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,Standard_Complex_Solutions.Length_Of(sols),
             natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars); flush(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
    startmoment := Ada.Calendar.Clock;
    Standard_Multitasked_Tracker(nbt,hom,abh,sols,pars,mhom,idz,verbose);
    stopmoment := Ada.Calendar.Clock;
    new_line(file);
    if arth and hcrd then
      if mhom = 1 then
        put_line(file,"THE 1-HOMOGENEOUS SOLUTIONS :");
      else
        put(file,"THE "); put(file,mhom,1);
        put_line(file,"-HOMOGENEOUS SOLUTIONS :");
      end if;  
    else
      put_line(file,"THE SOLUTIONS :");
    end if;
    put(file,Standard_Complex_Solutions.Length_Of(sols),
        natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
    if arth and hcrd then
      if mhom = 1 then
        Projective_Transformations.Affine_Transformation(sols);
      else
        Multi_Projective_Transformations.Make_Affine
          (sols,natural32(mhom),idz.all);
      end if;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Standard_Complex_Solutions.Length_Of(sols),
               natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
    new_line(file);
    put(file,"Elapsed wall clock time with ");
    put(file,nbt,1); put_line(file," tasks :");
    Time_Stamps.Write_Elapsed_Time(file,startmoment,stopmoment);
  end Track;

  procedure Track
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbt,mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean ) is

    ans : character;
    hcrd : constant boolean := (mhom > 0);
    verbose : boolean;
    startmoment,stopmoment : Ada.Calendar.Time;

  begin
    if not arth then
      put(file,natural32(hom.neq),natural32(hom.neq+1),
               DoblDobl_Homotopy.Homotopy_System);
    else
      declare
        p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
          := DoblDobl_Homotopy.Target_System;
        q : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
          := DoblDobl_Homotopy.Start_System;
      begin
        put(file,p'last,1); new_line(file); put(file,p);
        new_line(file);
        put_line(file,"THE START SYSTEM :");
        put(file,q'last,1); new_line(file); put(file,q);
      end;
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
             natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars); flush(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
    startmoment := Ada.Calendar.Clock;
    DoblDobl_Multitasked_Tracker(nbt,hom,abh,sols,pars,mhom,idz,verbose);
    stopmoment := Ada.Calendar.Clock;
    new_line(file);
    if arth and hcrd then
      if mhom = 1 then
        put_line(file,"THE 1-HOMOGENEOUS SOLUTIONS :");
      else
        put(file,"THE "); put(file,mhom,1);
        put_line(file,"-HOMOGENEOUS SOLUTIONS :");
      end if;
    else
      put_line(file,"THE SOLUTIONS :");
    end if;
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
             natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
    if arth and hcrd then
      if mhom = 1 then
        Projective_Transformations.Affine_Transformation(sols);
      else
        Multi_Projective_Transformations.Make_Affine
          (sols,natural32(mhom),idz.all);
      end if;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
               natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
    new_line(file);
    put(file,"Elapsed wall clock time with ");
    put(file,nbt,1); put_line(file," tasks :");
    Time_Stamps.Write_Elapsed_Time(file,startmoment,stopmoment);
  end Track;

  procedure Track
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbt,mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean ) is

    ans : character;
    hcrd : constant boolean := (mhom > 0);
    verbose : boolean;
    startmoment,stopmoment : Ada.Calendar.Time;

  begin
    if not arth then
      put(file,natural32(hom.neq),natural32(hom.neq+1),
               QuadDobl_Homotopy.Homotopy_System);
    else
      declare
        p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
          := QuadDobl_Homotopy.Target_System;
        q : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
          := QuadDobl_Homotopy.Start_System;
      begin
        put(file,p'last,1); new_line(file); put(file,p);
        new_line(file);
        put_line(file,"THE START SYSTEM :");
        put(file,q'last,1); new_line(file); put(file,q);
      end;
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
             natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars); flush(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
    startmoment := Ada.Calendar.Clock;
    QuadDobl_Multitasked_Tracker(nbt,hom,abh,sols,pars,mhom,idz,verbose);
    stopmoment := Ada.Calendar.Clock;
    new_line(file);
    if arth and hcrd then
      if mhom = 1 then
        put_line(file,"THE 1-HOMOGENEOUS SOLUTIONS :");
      else
        put(file,"THE "); put(file,mhom,1);
        put_line(file,"-HOMOGENEOUS SOLUTIONS :");
      end if;
    else
      put_line(file,"THE SOLUTIONS :");
    end if;
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
             natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
    if arth and hcrd then
      if mhom = 1 then
        Projective_Transformations.Affine_Transformation(sols);
      else
        Multi_Projective_Transformations.Make_Affine
          (sols,natural32(mhom),idz.all);
      end if;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
               natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
    new_line(file);
    put(file,"Elapsed wall clock time with ");
    put(file,nbt,1); put_line(file," tasks :");
    Time_Stamps.Write_Elapsed_Time(file,startmoment,stopmoment);
  end Track;

  procedure Standard_Main ( vrb : in integer32 := 0 ) is

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    pars : Homotopy_Continuation_Parameters.Parameters;
    nbt : integer32 := 0;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    artificial : boolean;
    file : file_type;
    start_moment,ended_moment : Ada.Calendar.Time;

  begin
    start_moment := Ada.Calendar.Clock;
    if vrb > 0
     then put_line("-> in track_path_convolutions.Standard_Main ...");
    end if;
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    new_line;
    put("Give the number of tasks : "); get(nbt); skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    Track(file,cnvhom,abshom,sols,pars,nbt,integer32(mhom),idz,artificial);
    ended_moment := Ada.Calendar.Clock;
    put(file,"PHC ran from "); Write_Time_Stamp(file,start_moment);
    put(file," till "); Write_Time_Stamp(file,ended_moment);
    put_line(file,".");
    put(file,"Number of tasks used in this run : ");
    put(file,nbt,1); new_line(file);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end Standard_Main;

  procedure DoblDobl_Main ( vrb : in integer32 := 0 ) is

    sols : DoblDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    pars : Homotopy_Continuation_Parameters.Parameters;
    nbt : integer32 := 0;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    artificial : boolean;
    file : file_type;
    start_moment,ended_moment : Ada.Calendar.Time;

  begin
    start_moment := Ada.Calendar.Clock;
    if vrb > 0
     then put_line("-> in track_path_convolutions.DoblDobl_Main ...");
    end if;
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    new_line;
    put("Give the number of tasks : "); get(nbt); skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    Track(file,cnvhom,abshom,sols,pars,nbt,integer32(mhom),idz,artificial);
    ended_moment := Ada.Calendar.Clock;
    put(file,"PHC ran from "); Write_Time_Stamp(file,start_moment);
    put(file," till "); Write_Time_Stamp(file,ended_moment);
    put_line(file,".");
    put(file,"Number of tasks used in this run : ");
    put(file,nbt,1); new_line(file);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end DoblDobl_Main;

  procedure QuadDobl_Main ( vrb : in integer32 := 0 ) is

    sols : QuadDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    pars : Homotopy_Continuation_Parameters.Parameters;
    nbt : integer32 := 0;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    artificial : boolean;
    file : file_type;
    start_moment,ended_moment : Ada.Calendar.Time;

  begin
    start_moment := Ada.Calendar.Clock;
    if vrb > 0
     then put_line("-> in track_path_convolutions.QuadDobl_Main ...");
    end if;
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    new_line;
    put("Give the number of tasks : "); get(nbt); skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    Track(file,cnvhom,abshom,sols,pars,nbt,integer32(mhom),idz,artificial);
    ended_moment := Ada.Calendar.Clock;
    put(file,"PHC ran from "); Write_Time_Stamp(file,start_moment);
    put(file," till "); Write_Time_Stamp(file,ended_moment);
    put_line(file,".");
    put(file,"Number of tasks used in this run : ");
    put(file,nbt,1); new_line(file);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end QuadDobl_Main;

end Multitasked_Path_Convolutions;
