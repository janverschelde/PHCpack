with Timing_Package;                     use Timing_Package;
with Characters_and_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with Residual_Convolution_Circuits;      use Residual_Convolution_Circuits;
with Shift_Convolution_Circuits;
with Corrector_Convolutions;             use Corrector_Convolutions;
with Standard_Pade_Trackers;
with Series_and_Trackers;
with Hyperplane_Convolution_Scaling;

package body Predictor_Corrector_Loops is

  procedure Predictor_Corrector_Loop
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; hcrd : in boolean;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_float;
                step,mixres : out double_float; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean ) is

    use Standard_Predictor_Convolutions;

    info : integer32 := 0;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
    if hcrd then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    end if;
    SVD_Prediction(hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm);
    if pars.corsteps = 0 then -- no corrector for zero pars.corsteps
      mixres := 1.0; nbrit := 0;
    else
      loop
        Store_Leading_Coefficients(hom.crc,homlead);
        Store_Leading_Coefficients(abh.crc,abhlead);
        Step_Coefficient(hom,step);
        Update_Radii_of_Constants(abh,hom);
        LU_Newton_Steps(hom,abh,psv,integer32(pars.corsteps),nbrit,
                        pars.tolres,mixres,dx,ipvt,info,fail);
        Restore_Leading_Coefficients(homlead,hom.crc);
        Restore_Leading_Coefficients(abhlead,abh.crc);
        exit when not fail;   
        step := step/2.0;
        exit when (step < pars.minsize);
        Standard_Rational_Approximations.Evaluate
          (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
      end loop;
    end if;
  end Predictor_Corrector_Loop;

  procedure Predictor_Corrector_Loop
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; hcrd : in boolean;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_float;
                step,mixres : out double_float; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true ) is

    use Standard_Predictor_Convolutions;

    info : integer32 := 0;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
    if hcrd then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    end if;
    SVD_Prediction(file,hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm,false,verbose);
    if verbose then
      if fail
       then put(file,"Predictor failed to reach tolerance");
       else put(file,"Predictor reached tolerance");
      end if;
      put(file,pars.alpha,3);
      put(file," at t :"); put(file,acct,3); put_line(file,".");
    end if;
    if pars.corsteps = 0 then -- no corrector for zero pars.corsteps
      mixres := 1.0; nbrit := 0;
    else
      loop
        Store_Leading_Coefficients(hom.crc,homlead);
        Store_Leading_Coefficients(abh.crc,abhlead);
        Step_Coefficient(hom,step);
        Update_Radii_of_Constants(abh,hom);
        LU_Newton_Steps(file,hom,abh,psv,integer32(pars.corsteps),nbrit,
                        pars.tolres,mixres,dx,ipvt,info,fail,verbose);
        Restore_Leading_Coefficients(homlead,hom.crc);
        Restore_Leading_Coefficients(abhlead,abh.crc);
        exit when not fail;   
        step := step/2.0;
        if verbose then
          put(file,"Reduced step size to"); put(file,step,3);
          put_line(file,".");
        end if;
        exit when (step < pars.minsize);
        Standard_Rational_Approximations.Evaluate
          (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
      end loop;
    end if;
  end Predictor_Corrector_Loop;

  procedure Predictor_Corrector_Loop
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; hcrd : in boolean;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_double;
                step,mixres : out double_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean ) is

    use DoblDobl_Predictor_Convolutions;

    info : integer32 := 0;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
    if hcrd then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    end if;
    SVD_Prediction(hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm);
    if pars.corsteps = 0 then -- no corrector for zero pars.corsteps
      mixres := create(1.0); nbrit := 0;
    else
      loop
        Store_Leading_Coefficients(hom.crc,homlead);
        Store_Leading_Coefficients(abh.crc,abhlead);
        Step_Coefficient(hom,step);
        Update_Radii_of_Constants(abh,hom);
        LU_Newton_Steps(hom,abh,psv,integer32(pars.corsteps),nbrit,
                        pars.tolres,mixres,dx,ipvt,info,fail);
        Restore_Leading_Coefficients(homlead,hom.crc);
        Restore_Leading_Coefficients(abhlead,abh.crc);
        exit when not fail;
        step := step/2.0;
        exit when (step < pars.minsize);
        DoblDobl_Rational_Approximations.Evaluate
          (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
      end loop;
    end if;
  end Predictor_Corrector_Loop;

  procedure Predictor_Corrector_Loop
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; hcrd : in boolean;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_double;
                step,mixres : out double_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true ) is

    use DoblDobl_Predictor_Convolutions;

    info : integer32 := 0;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
    if hcrd then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    end if;
    SVD_Prediction(file,hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm,false,verbose);
    if verbose then
      if fail
       then put(file,"Predictor failed to reach tolerance");
       else put(file,"Predictor reached tolerance");
      end if;
      put(file,pars.alpha,3);
      put(file," at t : "); put(file,acct,3); put_line(file,".");
    end if;
    if pars.corsteps = 0 then -- no corrector for zero pars.corsteps
      mixres := create(1.0); nbrit := 0;
    else
      loop
        Store_Leading_Coefficients(hom.crc,homlead);
        Store_Leading_Coefficients(abh.crc,abhlead);
        Step_Coefficient(hom,step);
        Update_Radii_of_Constants(abh,hom);
        LU_Newton_Steps(file,hom,abh,psv,integer32(pars.corsteps),nbrit,
                        pars.tolres,mixres,dx,ipvt,info,fail,verbose);
        Restore_Leading_Coefficients(homlead,hom.crc);
        Restore_Leading_Coefficients(abhlead,abh.crc);
        exit when not fail;
        step := step/2.0;
        if verbose then
          put(file,"Reduced step size to "); put(file,step,3);
          put_line(file,".");
        end if;
        exit when (step < pars.minsize);
        DoblDobl_Rational_Approximations.Evaluate
          (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
      end loop;
    end if;
  end Predictor_Corrector_Loop;

  procedure Predictor_Corrector_Loop
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; hcrd : in boolean;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out quad_double;
                step,mixres : out quad_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean ) is

    use QuadDobl_Predictor_Convolutions;

    info : integer32 := 0;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
    if hcrd then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    end if;
    SVD_Prediction(hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm);
    if pars.corsteps = 0 then -- no corrector for zero pars.corsteps
      mixres := create(1.0); nbrit := 0;
    else
      loop
        Store_Leading_Coefficients(hom.crc,homlead);
        Store_Leading_Coefficients(abh.crc,abhlead);
        Step_Coefficient(hom,step);
        Update_Radii_of_Constants(abh,hom);
        LU_Newton_Steps(hom,abh,psv,integer32(pars.corsteps),nbrit,
                        pars.tolres,mixres,dx,ipvt,info,fail);
        Restore_Leading_Coefficients(homlead,hom.crc);
        Restore_Leading_Coefficients(abhlead,abh.crc);
        exit when not fail;
        step := step/2.0;
        exit when (step < pars.minsize);
        QuadDobl_Rational_Approximations.Evaluate
          (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
      end loop;
    end if;
  end Predictor_Corrector_Loop;

  procedure Predictor_Corrector_Loop
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; hcrd : in boolean;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out quad_double;
                step,mixres : out quad_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true ) is

    use QuadDobl_Predictor_Convolutions;

    info : integer32 := 0;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
    if hcrd then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    end if;
    SVD_Prediction(file,hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm,false,verbose);
    if verbose then
      if fail
       then put(file,"Predictor failed to reach tolerance");
       else put(file,"Predictor reached tolerance");
      end if;
      put(file,pars.alpha,3);
      put(file," at t : "); put(file,acct,3); put_line(file,".");
    end if;
    if pars.corsteps = 0 then -- no corrector for zero pars.corsteps
      mixres := create(1.0); nbrit := 0;
    else
      loop
        Store_Leading_Coefficients(hom.crc,homlead);
        Store_Leading_Coefficients(abh.crc,abhlead);
        Step_Coefficient(hom,step);
        Update_Radii_of_Constants(abh,hom);
        LU_Newton_Steps(file,hom,abh,psv,integer32(pars.corsteps),nbrit,
                        pars.tolres,mixres,dx,ipvt,info,fail,verbose);
        Restore_Leading_Coefficients(homlead,hom.crc);
        Restore_Leading_Coefficients(abhlead,abh.crc);
        exit when not fail;
        step := step/2.0;
        if verbose then
          put(file,"Reduced step size to "); put(file,step,3);
          put_line(file,".");
        end if;
        exit when (step < pars.minsize);
        QuadDobl_Rational_Approximations.Evaluate
          (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
      end loop;
    end if;
  end Predictor_Corrector_Loop;

  procedure Track_One_Path
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out double_float;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float;
                fail : out boolean ) is

    endt : constant double_float := 1.0;
    step : double_float := 0.0;
    togo : double_float; 
    nbrit : integer32;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    minstpz := 1.0; maxstpz := 0.0; tnbrit := 0;
    for k in 1..pars.maxsteps loop
      Predictor_Corrector_Loop(hom,abh,homlead,abhlead,pars,maxit,false,prd,
        psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,nbpole,nbhess,nbmaxm,fail);
      tnbrit := tnbrit + natural32(nbrit);
      Standard_Pade_Trackers.Update_Step_Sizes(minstpz,maxstpz,step);
      togo := endt - acct;
      if (abs(togo) < pars.epsilon)
       then nbsteps := k; exit;
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
    end loop;
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out double_float;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float;
                fail : out boolean; verbose : in boolean := true ) is

    endt : constant double_float := 1.0;
    step : double_float := 0.0;
    togo : double_float; 
    nbrit : integer32;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    minstpz := 1.0; maxstpz := 0.0; tnbrit := 0;
    for k in 1..pars.maxsteps loop
      if verbose
       then put(file,"t :"); put(file,acct,3); put_line(file," :");
      end if;
      Predictor_Corrector_Loop(file,hom,abh,homlead,abhlead,pars,maxit,false,
        prd,psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,nbpole,nbhess,nbmaxm,
        fail,verbose);
      tnbrit := tnbrit + natural32(nbrit);
      Standard_Pade_Trackers.Update_Step_Sizes(minstpz,maxstpz,step);
      if verbose then
        if fail
         then put_line(file,"Predictor-Corrector loop failed.");
         else put_line(file,"Predictor-Corrector loop succeeded.");
        end if;
      end if;
      togo := endt - acct;
      if (abs(togo) < pars.epsilon)
       then nbsteps := k; exit;
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
    end loop;
  end Track_One_Path;

  procedure Track_One_Path
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out double_double;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float;
                fail : out boolean ) is

    endt : constant double_float := 1.0;
    step : double_double := create(0.0);
    togo : double_double; 
    nbrit : integer32;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    minstpz := 1.0; maxstpz := 0.0; tnbrit := 0;
    for k in 1..pars.maxsteps loop
      Predictor_Corrector_Loop(hom,abh,homlead,abhlead,pars,maxit,false,prd,
        psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,nbpole,nbhess,nbmaxm,fail);
      tnbrit := tnbrit + natural32(nbrit);
      Standard_Pade_Trackers.Update_Step_Sizes(minstpz,maxstpz,hi_part(step));
      togo := endt - acct;
      if (abs(togo) < pars.epsilon)
       then nbsteps := k; exit;
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
    end loop;
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out double_double;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float;
                fail : out boolean; verbose : in boolean := true ) is

    endt : constant double_float := 1.0;
    step : double_double := create(0.0);
    togo : double_double; 
    nbrit : integer32;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    minstpz := 1.0; maxstpz := 0.0; tnbrit := 0;
    for k in 1..pars.maxsteps loop
      if verbose
       then put(file,"t : "); put(file,acct,3); put_line(file," :");
      end if;
      Predictor_Corrector_Loop(file,hom,abh,homlead,abhlead,pars,maxit,false,
        prd,psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,nbpole,nbhess,nbmaxm,
        fail,verbose);
      tnbrit := tnbrit + natural32(nbrit);
      Standard_Pade_Trackers.Update_Step_Sizes(minstpz,maxstpz,hi_part(step));
      if verbose then
        if fail
         then put_line(file,"Predictor-Corrector loop failed.");
         else put_line(file,"Predictor-Corrector loop succeeded.");
        end if;
      end if;
      togo := endt - acct;
      if (abs(togo) < pars.epsilon)
       then nbsteps := k; exit;
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
    end loop;
  end Track_One_Path;

  procedure Track_One_Path
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out quad_double;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float; fail : out boolean ) is

    endt : constant double_float := 1.0;
    step : quad_double := create(0.0);
    togo : quad_double; 
    nbrit : integer32;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    minstpz := 1.0; maxstpz := 0.0; tnbrit := 0;
    for k in 1..pars.maxsteps loop
      Predictor_Corrector_Loop(hom,abh,homlead,abhlead,pars,maxit,false,prd,
        psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,nbpole,nbhess,nbmaxm,fail);
      tnbrit := tnbrit + natural32(nbrit);
      Standard_Pade_Trackers.Update_Step_Sizes(minstpz,maxstpz,hihi_part(step));
      togo := endt - acct;
      if (abs(togo) < pars.epsilon)
       then nbsteps := k; exit;
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
    end loop;
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                acct,mixres : in out quad_double;
                tnbrit,nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                minstpz,maxstpz : out double_float;
                fail : out boolean; verbose : in boolean := true ) is

    endt : constant double_float := 1.0;
    step : quad_double := create(0.0);
    togo : quad_double; 
    nbrit : integer32;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    minstpz := 1.0; maxstpz := 0.0; tnbrit := 0;
    for k in 1..pars.maxsteps loop
      if verbose
       then put(file,"t : "); put(file,acct,3); put_line(file," :");
      end if;
      Predictor_Corrector_Loop(file,hom,abh,homlead,abhlead,pars,maxit,false,
        prd,psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,nbpole,nbhess,nbmaxm,
        fail,verbose);
      tnbrit := tnbrit + natural32(nbrit);
      Standard_Pade_Trackers.Update_Step_Sizes(minstpz,maxstpz,hihi_part(step));
      if verbose then
        if fail
         then put_line(file,"Predictor-Corrector loop failed.");
         else put_line(file,"Predictor-Corrector loop succeeded.");
        end if;
      end if;
      togo := endt - acct;
      if (abs(togo) < pars.epsilon)
       then nbsteps := k; exit;
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
    end loop;
  end Track_One_Path;

  procedure Track_All_Paths
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := true ) is

    use Standard_Complex_Solutions,Standard_Predictor_Convolutions;

    timer : Timing_Widget;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxit : constant integer32 := (numdeg + dendeg + 2)/2;
    ls : Link_to_Solution := Head_Of(sols);
    prd : Predictor := Create(ls.v,hom.neq,hom.deg,numdeg,dendeg,SVD);
    psv : Predictor_Vectors(hom.dim,hom.neq);
    svh : Link_to_SVD_Hessians := Create(hom.dim);
    solsptr : Solution_List := sols;
    tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    nbrit : integer32;
    fail : boolean;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : Standard_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : Standard_Complex_VecVecs.Link_to_VecVec;
    wrk : Standard_Complex_Vectors.Link_to_Vector
        := Standard_Speelpenning_Convolutions.Allocate_Coefficients(hom.deg);
    homcff : Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
    acct,mixres : double_float := 0.0;
    pathno : natural32 := 0;
    ratpole,rathess,ratmaxm,minstpz,maxstpz : double_float := 0.0;
    lensols : constant natural32 := Length_Of(sols);
    minpastp,maxpastp : double_float;
    mincorsteps,maxcorsteps,minnbrsteps,maxnbrsteps : natural32;

  begin
    Allocate_Coefficients(hom.crc,homcff);
    Store_Coefficients(hom.crc,homcff);
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    minpastp := 1.0; maxpastp := 0.0;
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    mincorsteps := (pars.maxsteps+1)*pars.corsteps+1; maxcorsteps := 0;
    tstart(timer);
    while not Is_Null(solsptr) loop
      pathno := pathno + 1;
      put(file,"Path "); put(file,pathno,1); put_line(file," :");
      ls := Head_Of(solsptr); psv.sol := ls.v; acct := 0.0;
      Track_One_Path(file,hom,abh,homlead,abhlead,pars,maxit,prd,psv,svh,
        dx,ipvt,wrk,acct,mixres,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,
        minstpz,maxstpz,fail,verbose);
      Restore_Coefficients(homcff,hom.crc);
      Update_Radii_of_Constants(abh,hom);
      Step_Coefficient(hom,acct);
      LU_Newton_Steps(file,hom,abh,psv,1,nbrit,pars.tolres,mixres,dx,
                      ipvt,ls.rco,fail,verbose);
      ls.v := psv.sol; ls.res := mixres;
      ls.err := Standard_Complex_Vector_Norms.Max_Norm(dx);
      ls.t := Standard_Complex_Numbers.Create(acct); Set_Head(solsptr,ls);
      put(file,ls.all); new_line(file);
      Write_Path_Statistics
        (file,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz);
      Series_and_Trackers.Update_Ratio_Sums
        (ratpole,rathess,ratmaxm,nbpole,nbhess,nbmaxm,nbsteps*lensols);
      Series_and_Trackers.Update_MinMax(minpastp,maxpastp,minstpz,maxstpz);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbsteps);
      Series_and_Trackers.Update_Counters(mincorsteps,maxcorsteps,tnbrit);
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      Restore_Leading_Coefficients(abhlead,abh.crc);
      Restore_Coefficients(homcff,hom.crc);
    end loop;
    tstop(timer);
    Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,mincorsteps,maxcorsteps,
       ratpole,rathess,ratmaxm,minpastp,maxpastp);
    new_line(file);
    print_times(file,timer,"tracking "
                & Characters_and_Numbers.nConvert(lensols)
                & " paths in double precision");
    Clear(prd); Clear(svh);
    Standard_Complex_Vectors.Clear(wrk);
    Standard_Complex_VecVecs.Deep_Clear(homlead);
    Standard_Complex_VecVecs.Deep_Clear(abhlead);
    Standard_Speelpenning_Convolutions.Clear(homcff);
  end Track_All_Paths;

  procedure Track_All_Paths
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := true ) is

    use DoblDobl_Complex_Solutions,DoblDobl_Predictor_Convolutions;

    timer : Timing_Widget;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxit : constant integer32 := (numdeg + dendeg + 2)/2;
    ls : Link_to_Solution := Head_Of(sols);
    prd : Predictor := Create(ls.v,hom.neq,hom.deg,numdeg,dendeg,SVD);
    psv : Predictor_Vectors(hom.dim,hom.neq);
    svh : Link_to_SVD_Hessians := Create(hom.dim);
    solsptr : Solution_List := sols;
    tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    nbrit : integer32 := 0;
    fail : boolean;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : DoblDobl_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    wrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := DoblDobl_Speelpenning_Convolutions.Allocate_Coefficients(hom.deg);
    homcff : DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
    acct,mixres : double_double;
    pathno : natural32 := 0;
    ratpole,rathess,ratmaxm,minstpz,maxstpz : double_float := 0.0;
    lensols : constant natural32 := Length_Of(sols);
    minpastp,maxpastp : double_float;
    mincorsteps,maxcorsteps,minnbrsteps,maxnbrsteps : natural32;

  begin
    Allocate_Coefficients(hom.crc,homcff);
    Store_Coefficients(hom.crc,homcff);
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    minpastp := 1.0; maxpastp := 0.0;
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    mincorsteps := (pars.maxsteps+1)*pars.corsteps+1; maxcorsteps := 0;
    tstart(timer);
    while not Is_Null(solsptr) loop
      pathno := pathno + 1;
      put(file,"Path "); put(file,pathno,1); put_line(file," :");
      ls := Head_Of(solsptr); psv.sol := ls.v; acct := create(0.0);
      Track_One_Path(file,hom,abh,homlead,abhlead,pars,maxit,prd,psv,svh,
        dx,ipvt,wrk,acct,mixres,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,
        minstpz,maxstpz,fail,verbose);
      Restore_Coefficients(homcff,hom.crc);
      Update_Radii_of_Constants(abh,hom);
      Step_Coefficient(hom,acct);
      LU_Newton_Steps(file,hom,abh,psv,1,nbrit,pars.tolres,mixres,dx,
                      ipvt,ls.rco,fail,verbose);
      ls.v := psv.sol; ls.res := mixres;
      ls.err := DoblDobl_Complex_Vector_Norms.Max_Norm(dx);
      ls.t := DoblDobl_Complex_Numbers.Create(acct); Set_Head(solsptr,ls);
      put(file,ls.all); new_line(file);
      Write_Path_Statistics
        (file,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz);
      Series_and_Trackers.Update_Ratio_Sums
        (ratpole,rathess,ratmaxm,nbpole,nbhess,nbmaxm,nbsteps*lensols);
      Series_and_Trackers.Update_MinMax(minpastp,maxpastp,minstpz,maxstpz);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbsteps);
      Series_and_Trackers.Update_Counters(mincorsteps,maxcorsteps,tnbrit);
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      Restore_Leading_Coefficients(abhlead,abh.crc);
      Restore_Coefficients(homcff,hom.crc);
    end loop;
    tstop(timer);
    Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,mincorsteps,maxcorsteps,
       ratpole,rathess,ratmaxm,minpastp,maxpastp);
    new_line(file);
    print_times(file,timer,"tracking "
                & Characters_and_Numbers.nConvert(lensols)
                & " paths in double precision");
    Clear(prd); Clear(svh);
    DoblDobl_Complex_Vectors.Clear(wrk);
    DoblDobl_Complex_VecVecs.Deep_Clear(homlead);
    DoblDobl_Complex_VecVecs.Deep_Clear(abhlead);
    DoblDobl_Speelpenning_Convolutions.Clear(homcff);
  end Track_All_Paths;

  procedure Track_All_Paths
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := true ) is

    use QuadDobl_Complex_Solutions,QuadDobl_Predictor_Convolutions;

    timer : Timing_Widget;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxit : constant integer32 := (numdeg + dendeg + 2)/2;
    ls : Link_to_Solution := Head_Of(sols);
    prd : Predictor := Create(ls.v,hom.neq,hom.deg,numdeg,dendeg,SVD);
    psv : Predictor_Vectors(hom.dim,hom.neq);
    svh : Link_to_SVD_Hessians := Create(hom.dim);
    solsptr : Solution_List := sols;
    tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    nbrit : integer32 := 0;
    fail : boolean;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : QuadDobl_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    wrk : QuadDobl_Complex_Vectors.Link_to_Vector
        := QuadDobl_Speelpenning_Convolutions.Allocate_Coefficients(hom.deg);
    homcff : QuadDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
    acct,mixres : quad_double;
    pathno : natural32 := 0;
    ratpole,rathess,ratmaxm,minstpz,maxstpz : double_float := 0.0;
    lensols : constant natural32 := Length_Of(sols);
    minpastp,maxpastp : double_float;
    mincorsteps,maxcorsteps,minnbrsteps,maxnbrsteps : natural32;

  begin
    Allocate_Coefficients(hom.crc,homcff);
    Store_Coefficients(hom.crc,homcff);
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    minpastp := 1.0; maxpastp := 0.0;
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    mincorsteps := (pars.maxsteps+1)*pars.corsteps+1; maxcorsteps := 0;
    tstart(timer);
    while not Is_Null(solsptr) loop
      pathno := pathno + 1;
      put(file,"Path "); put(file,pathno,1); put_line(file," :");
      ls := Head_Of(solsptr); psv.sol := ls.v; acct := Create(0.0);
      Track_One_Path(file,hom,abh,homlead,abhlead,pars,maxit,prd,psv,svh,
        dx,ipvt,wrk,acct,mixres,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,
        minstpz,maxstpz,fail,verbose);
      Restore_Coefficients(homcff,hom.crc);
      Update_Radii_of_Constants(abh,hom);
      Step_Coefficient(hom,acct);
      LU_Newton_Steps(file,hom,abh,psv,1,nbrit,pars.tolres,mixres,dx,
                      ipvt,ls.rco,fail,verbose);
      ls.v := psv.sol; ls.res := mixres;
      ls.err := QuadDobl_Complex_Vector_Norms.Max_Norm(dx);
      ls.t := QuadDobl_Complex_Numbers.Create(acct); Set_Head(solsptr,ls);
      put(file,ls.all); new_line(file);
      Write_Path_Statistics
        (file,tnbrit,nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz);
      Series_and_Trackers.Update_Ratio_Sums
        (ratpole,rathess,ratmaxm,nbpole,nbhess,nbmaxm,nbsteps*lensols);
      Series_and_Trackers.Update_MinMax(minpastp,maxpastp,minstpz,maxstpz);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbsteps);
      Series_and_Trackers.Update_Counters(mincorsteps,maxcorsteps,tnbrit);
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      Restore_Leading_Coefficients(abhlead,abh.crc);
      Restore_Coefficients(homcff,hom.crc);
    end loop;
    tstop(timer);
    Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,mincorsteps,maxcorsteps,
       ratpole,rathess,ratmaxm,minpastp,maxpastp);
    new_line(file);
    print_times(file,timer,"tracking "
                & Characters_and_Numbers.nConvert(lensols)
                & " paths in double precision");
    Clear(prd); Clear(svh);
    QuadDobl_Complex_Vectors.Clear(wrk);
    QuadDobl_Complex_VecVecs.Deep_Clear(homlead);
    QuadDobl_Complex_VecVecs.Deep_Clear(abhlead);
    QuadDobl_Speelpenning_Convolutions.Clear(homcff);
  end Track_All_Paths;

  procedure Write_Path_Statistics
              ( file : in file_type;
                nbrit,nbpole,nbhess,nbmaxm,nbsteps : in natural32;
                minstpz,maxstpz : in double_float ) is
  begin
    put(file,"The smallest step size on the path :");
    put(file,minstpz,2); new_line(file);
    put(file,"The largest step size on the path  :");
    put(file,maxstpz,2); new_line(file);
    put(file,"The total number of steps on the path          : ");
    put(file,nbsteps,1); new_line(file);
    put(file,"The number of corrector iterations on the path : ");
    put(file,nbrit,1); new_line(file);
    put(file,"Number of times the pole step was minimal      : ");
    put(file,nbpole,1); new_line(file);
    put(file,"Number of times the curvature step was minimal : ");
    put(file,nbhess,1); new_line(file);
    put(file,"Number of times the maximum step was minimal   : ");
    put(file,nbmaxm,1); new_line(file);
  end Write_Path_Statistics;

  procedure Write_Total_Path_Statistics
              ( file : in file_type;
                minnbrsteps,maxnbrsteps : in natural32;
                mincorsteps,maxcorsteps : in natural32;
                ratpole,rathess,ratmaxm : in double_float;
                minstpz,maxstpz : in double_float ) is
  begin
    new_line(file);
    put(file,"The smallest number of total steps : ");
    put(file,minnbrsteps,1); new_line(file);
    put(file,"The largest number of total steps  : ");
    put(file,maxnbrsteps,1); new_line(file);
    put(file,"The smallest number of total corrector iterations : ");
    put(file,mincorsteps,1); new_line(file);
    put(file,"The largest number of total corrector iterations  : ");
    put(file,maxcorsteps,1); new_line(file);
    put(file,"The smallest step size on the path :");
    put(file,minstpz,2); new_line(file);
    put(file,"The largest step size on the path  :");
    put(file,maxstpz,2); new_line(file);
    put(file,"Average ratio of times pole step was minimal      : ");
    put(file,ratpole,1,4,0); new_line(file);
    put(file,"Average ratio of times curvature step was minimal : ");
    put(file,rathess,1,4,0); new_line(file);
    put(file,"Average ratio of times maximum step was minimal   : ");
    put(file,ratmaxm,1,4,0); new_line(file);
  end Write_Total_Path_Statistics;

end Predictor_Corrector_Loops;
