with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vector_Norms;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;
with Corrector_Convolutions;             use Corrector_Convolutions;
with Predictor_Corrector_Loops;          use Predictor_Corrector_Loops;
with Standard_Solutions_Queue;
with DoblDobl_Solutions_Queue;
with QuadDobl_Solutions_Queue;
with Multitasking;

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

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      t,mixres : double_float := 0.0;
      nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
      fail : boolean;

    begin
      loop
        myptr := Standard_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := 0.0;
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
          maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
          wrk(i),t,mixres,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.err := Standard_Complex_Vector_Norms.Max_Norm(dx(i).all);
        ls.res := mixres;
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
      t,mixres : double_float := 0.0;
      nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
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
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
          maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
          wrk(i),t,mixres,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.err := Standard_Complex_Vector_Norms.Max_Norm(dx(i).all);
        ls.res := mixres;
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

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      t,mixres : double_double;
      nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
      fail : boolean;

    begin
      loop
        myptr := DoblDobl_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := Create(0.0);
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
          maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
          wrk(i),t,mixres,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.err := DoblDobl_Complex_Vector_Norms.Max_Norm(dx(i).all);
        ls.res := mixres;
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
      nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
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
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
          maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
          wrk(i),t,mixres,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.err := DoblDobl_Complex_Vector_Norms.Max_Norm(dx(i).all);
        ls.res := mixres;
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
      Copy(hom,homsa(k)); Copy(hom,abhsa(k));
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

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      t,mixres : quad_double;
      nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
      fail : boolean;

    begin
      loop
        myptr := QuadDobl_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        psv(i).sol := ls.v; t := Create(0.0);
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
          maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
          wrk(i),t,mixres,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.err := QuadDobl_Complex_Vector_Norms.Max_Norm(dx(i).all);
        ls.res := mixres;
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
      nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
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
        Track_One_Path(homsa(i),abhsa(i),homlead(i),abhlead(i),pars,
          maxit,prd(i),psv(i).all,svh(i),dx(i).all,ipvt(i).all,
          wrk(i),t,mixres,nbpole,nbhess,nbmaxm,nbsteps,fail);
        ls.err := QuadDobl_Complex_Vector_Norms.Max_Norm(dx(i).all);
        ls.res := mixres;
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
      Copy(hom,homsa(k)); Copy(hom,abhsa(k));
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

end Multitasked_Path_Convolutions;
