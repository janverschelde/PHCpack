with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;                        use Time_Stamps;
with Write_Seed_Number;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Vector_Splitters;
with Standard_Complex_VecMats;
with Standard_Complex_Circuits;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;
with Residual_Convolution_Circuits;
with Standard_Circuit_Makers;
with Standard_Convolution_Splitters;
with Standard_Newton_Circuits;
with Standard_Coefficient_Storage;
with Corrector_Convolutions;             use Corrector_Convolutions;
with Predictor_Corrector_Trackers;       use Predictor_Corrector_Trackers;
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

  procedure Allocate ( v : in out Standard_Floating_VecVecs.VecVec;
                       n : in integer32 ) is
  begin
    for k in v'range loop
      v(k) := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
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
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                verbose : in boolean := true ) is

    use Standard_Complex_Solutions;
    use Standard_Predictor_Convolutions;

    maxit : constant integer32 := integer32(pars.numdeg + pars.dendeg + 2)/2;
    homsa : Standard_Coefficient_Convolutions.System_Array(1..nbtasks);
    rcf,icf : Standard_Floating_VecVecVecs.Link_to_VecVecVec;
    rx,ix : Standard_Floating_VecVecVecs.VecVecVec(1..nbtasks);
    cfhsa,abhsa : Standard_Coefficient_Circuits.System_Array(1..nbtasks);
    ipvt : Standard_Integer_VecVecs.VecVec(1..nbtasks);
    first : constant Link_to_Solution := Head_Of(sols);
    prd : Predictor_Array(1..nbtasks)
        := Create(nbtasks,first.v,hom.neq,hom.deg,
                  integer32(pars.numdeg),integer32(pars.dendeg),LU); --SVD);
    psv : Predictor_Vectors_Array(1..nbtasks)
        := Create(nbtasks,hom.dim,hom.neq);
    svh : SVD_Hessians_Array(1..nbtasks) := Create(nbtasks,hom.dim);
    xtr : constant natural32 := 1;
    xr,xi,pwt : Standard_Floating_VecVecs.VecVec(1..nbtasks);
    vh : Standard_Complex_VecMats.VecMat_Array(1..nbtasks);
    svls : Standard_Speelpenning_Convolutions.VecVecVec(1..nbtasks);

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      t,initres,mixres,minstpz,maxstpz : double_float := 0.0;
      tnbrit,nbpole,nbhess,nbmaxm,nbsteps,nbrit : natural32 := 0;
      fail : boolean;

    begin
      loop
        myptr := Standard_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := 0.0;
        Predictor_Corrector_Trackers.Track_One_Path
          (homsa(i),rcf,icf,cfhsa(i),abhsa(i),pars,maxit,mhom,idz,
           prd(i),psv(i).all,svh(i),rx(i),ix(i),xr(i),xi(i),
           vh(i).all,svls(i).all,ipvt(i).all,pwt(i),t,mixres,
           tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz,fail);
        Standard_Coefficient_Storage.Restore(rcf,icf,homsa(i).crc);
        EvalCffRad(homsa(i),cfhsa(i),abhsa(i),t);
        Standard_Newton_Circuits.LU_Newton_Steps
          (cfhsa(i),abhsa(i),psv(i).sol,psv(i).radsol,xr(i),xi(i),
           pars.corsteps,pars.tolres,pars.tolres,ipvt(i).all,initres,
           ls.res,ls.rco,ls.err,mixres,nbrit,fail,xtr);
        tnbrit := tnbrit + nbrit;
        ls.v := psv(i).sol; ls.t := Standard_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
        Standard_Coefficient_Storage.Restore(rcf,icf,homsa(i).crc);
      end loop;
    end Silent_Track;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Track);

    procedure Report_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks with intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      cnt : integer32;
      t,initres,mixres,minstpz,maxstpz : double_float := 0.0;
      tnbrit,nbpole,nbhess,nbmaxm,nbsteps,nbrit : natural32 := 0;
      fail : boolean;

    begin
      loop
        myptr := Standard_Solutions_Queue.Next;
        cnt := Standard_Solutions_Queue.Next_Counter;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := 0.0;
        put_line("Task " & Multitasking.to_string(i)
                         & " tracks path "
                         & Multitasking.to_string(cnt));
        Predictor_Corrector_Trackers.Track_One_Path
          (homsa(i),rcf,icf,cfhsa(i),abhsa(i),pars,maxit,mhom,idz,
           prd(i),psv(i).all,svh(i),rx(i),ix(i),xr(i),xi(i),
           vh(i).all,svls(i).all,ipvt(i).all,pwt(i),t,mixres,
           tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz,fail);
        Standard_Coefficient_Storage.Restore(rcf,icf,homsa(i).crc);
        EvalCffRad(homsa(i),cfhsa(i),abhsa(i),t);
        Standard_Newton_Circuits.LU_Newton_Steps
          (cfhsa(i),abhsa(i),psv(i).sol,psv(i).radsol,xr(i),xi(i),
           pars.corsteps,pars.tolres,pars.tolres,ipvt(i).all,initres,
           ls.res,ls.rco,ls.err,mixres,nbrit,fail,xtr);
        tnbrit := tnbrit + nbrit;
        ls.v := psv(i).sol; ls.t := Standard_Complex_Numbers.Create(t);
       -- Set_Head(myptr,ls);
        Standard_Coefficient_Storage.Restore(rcf,icf,homsa(i).crc);
      end loop;
    end Report_Track;
    procedure report_do_jobs is 
      new Multitasking.Reporting_Workers(Report_Track);

  begin
    for k in 1..nbtasks loop
      Standard_Coefficient_Convolutions.Copy(hom,homsa(k));
      cfhsa(k) := Standard_Coefficient_Circuits.Copy(cfh);
      abhsa(k) := Standard_Coefficient_Circuits.Copy(abh);
      declare
        vm : constant Standard_Complex_VecMats.VecMat(1..hom.neq)
           := Standard_Complex_Circuits.Allocate(hom.neq,hom.dim);
        rz : constant Standard_Floating_Vectors.Vector(0..hom.deg)
           := (0..hom.deg => 0.0);
        sv : constant Standard_Complex_VecVecs.VecVec(0..hom.dim)
           := Standard_Vector_Splitters.Allocate(hom.neq,hom.dim+1,0,1);
      begin
        pwt(k) := new Standard_Floating_Vectors.Vector'(rz);
        rx(k) := Standard_Vector_Splitters.Allocate_Floating_Coefficients
                   (hom.dim,hom.deg);
        ix(k) := Standard_Vector_Splitters.Allocate_Floating_Coefficients
                   (hom.dim,hom.deg);
        vh(k) := new Standard_Complex_VecMats.VecMat'(vm);
        svls(k) := new Standard_Complex_VecVecs.VecVec'(sv);
      end;
    end loop;
    Standard_Coefficient_Storage.Allocate_and_Store(hom.crc,rcf,icf);
    Allocate(ipvt,hom.dim);
    Allocate(xr,hom.dim); Allocate(xi,hom.dim);
    Standard_Solutions_Queue.Initialize(sols);
    if verbose
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
    Standard_Coefficient_Convolutions.Clear(homsa);
    Standard_Coefficient_Circuits.Clear(cfhsa);
    Standard_Coefficient_Circuits.Clear(abhsa);
    Clear(prd); Clear(psv); Clear(svh);
    Standard_Integer_VecVecs.Clear(ipvt);
    Standard_Floating_VecVecs.Clear(pwt);
    Standard_Floating_VecVecs.Clear(xr);
    Standard_Floating_VecVecs.Clear(xi);
    Standard_Floating_VecVecVecs.Clear(rcf);
    Standard_Floating_VecVecVecs.Clear(icf);
    Standard_Floating_VecVecVecs.Clear(rx);
    Standard_Floating_VecVecVecs.Clear(ix);
    Standard_Complex_VecMats.Clear(vh);
  end Standard_Multitasked_Tracker;

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
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbt,mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean ) is

    verbose : boolean;
    startmoment,stopmoment : Ada.Calendar.Time;
  
  begin
    Track_Path_Convolutions.Standard_Write_Homotopy
      (file,hom.neq,sols,pars,arth,verbose);
    startmoment := Ada.Calendar.Clock;
    Standard_Multitasked_Tracker(nbt,hom,cfh,abh,sols,pars,mhom,idz,verbose);
    stopmoment := Ada.Calendar.Clock;
    Track_Path_Convolutions.Standard_Write_Solutions(file,arth,mhom,idz,sols);
    new_line(file);
    put(file,"Elapsed wall clock time with ");
    put(file,nbt,1); put_line(file," tasks :");
    Time_Stamps.Write_Elapsed_Time(file,startmoment,stopmoment);
  end Track;

  procedure Track
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbt,mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean ) is

    verbose : boolean;
    startmoment,stopmoment : Ada.Calendar.Time;
  
  begin
    Track_Path_Convolutions.Standard_Write_Homotopy
      (file,hom.neq,sols,pars,arth,verbose);
    startmoment := Ada.Calendar.Clock;
    Standard_Multitasked_Tracker(nbt,hom,abh,sols,pars,mhom,idz,verbose);
    stopmoment := Ada.Calendar.Clock;
    Track_Path_Convolutions.Standard_Write_Solutions(file,arth,mhom,idz,sols);
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

    verbose : boolean;
    startmoment,stopmoment : Ada.Calendar.Time;

  begin
    Track_Path_Convolutions.DoblDobl_Write_Homotopy
      (file,hom.neq,sols,pars,arth,verbose);
    startmoment := Ada.Calendar.Clock;
    DoblDobl_Multitasked_Tracker(nbt,hom,abh,sols,pars,mhom,idz,verbose);
    stopmoment := Ada.Calendar.Clock;
    Track_Path_Convolutions.DoblDobl_Write_Solutions(file,arth,mhom,idz,sols);
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

    verbose : boolean;
    startmoment,stopmoment : Ada.Calendar.Time;

  begin
    Track_Path_Convolutions.QuadDobl_Write_Homotopy
      (file,hom.neq,sols,pars,arth,verbose);
    startmoment := Ada.Calendar.Clock;
    QuadDobl_Multitasked_Tracker(nbt,hom,abh,sols,pars,mhom,idz,verbose);
    stopmoment := Ada.Calendar.Clock;
    Track_Path_Convolutions.QuadDobl_Write_Solutions(file,arth,mhom,idz,sols);
    new_line(file);
    put(file,"Elapsed wall clock time with ");
    put(file,nbt,1); put_line(file," tasks :");
    Time_Stamps.Write_Elapsed_Time(file,startmoment,stopmoment);
  end Track;

  procedure Standard_Main ( nt : in natural32; vrb : in integer32 := 0 ) is

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    pars : Homotopy_Continuation_Parameters.Parameters;
    nbt : integer32 := 0;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    artificial : boolean;
    file : file_type;
    start_moment,ended_moment : Ada.Calendar.Time;
    ans : character;

  begin
    if vrb > 0
     then put_line("-> in multitasked_path_convolutions.Standard_Main ...");
    end if;
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    if nt > 0 then
      nbt := integer32(nt);
    else
      new_line;
      put("Give the number of tasks : "); get(nbt); skip_line;
    end if;
    if mhom > 0 then
      ans := 'n'; -- no support for homogeneous coefficient circuits yet
    else
     -- new_line;
     -- put("Running with coefficient convolution circuits ? (y/n) ");
     -- Ask_Yes_or_No(ans);
      ans := 'y'; -- make coefficient convolution circuits the default
    end if;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    start_moment := Ada.Calendar.Clock;
    if ans = 'n' then
      Track(file,cnvhom,abshom,sols,pars,nbt,integer32(mhom),idz,artificial);
    else
      declare
        cffhom : Standard_Coefficient_Convolutions.Link_to_System;
        cfs,abh : Standard_Coefficient_Circuits.Link_to_System;
      begin
        cffhom := Standard_Convolution_Splitters.Split(cnvhom);
        cfs := Standard_Circuit_Makers.Make_Coefficient_System(cffhom);
        abh := Standard_Coefficient_Circuits.Copy(cfs);
        Standard_Coefficient_Circuits.AbsVal(abh);
        Track(file,cffhom,cfs,abh,sols,pars,nbt,integer32(mhom),idz,artificial);
      end;
    end if;
    ended_moment := Ada.Calendar.Clock;
    put(file,"PHC ran from "); Write_Time_Stamp(file,start_moment);
    put(file," till "); Write_Time_Stamp(file,ended_moment);
    put_line(file,".");
    put(file,"Number of tasks used in this run : ");
    put(file,nbt,1); new_line(file);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end Standard_Main;

  procedure DoblDobl_Main ( nt : in natural32; vrb : in integer32 := 0 ) is

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
     then put_line("-> in multitasked_path_convolutions.DoblDobl_Main ...");
    end if;
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    if nt > 0 then
      nbt := integer32(nt);
    else
      new_line;
      put("Give the number of tasks : "); get(nbt); skip_line;
    end if;
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

  procedure QuadDobl_Main ( nt : in natural32; vrb : in integer32 := 0 ) is

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
     then put_line("-> in multitasked_path_convolutions.QuadDobl_Main ...");
    end if;
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    if nt > 0 then
      nbt := integer32(nt);
    else
      new_line;
      put("Give the number of tasks : "); get(nbt); skip_line;
    end if;
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
