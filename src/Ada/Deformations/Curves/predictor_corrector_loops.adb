with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with Residual_Convolution_Circuits;      use Residual_Convolution_Circuits;
with Shift_Convolution_Circuits;
with Corrector_Convolutions;             use Corrector_Convolutions;

package body Predictor_Corrector_Loops is

  procedure Predictor_Corrector_Loop
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
                endt : in double_float; acct : in out double_float;
                step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean ) is

    use Standard_Predictor_Convolutions;

    info,nbrit : integer32 := 0;
    mixres : double_float;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
    SVD_Prediction(hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm);
    if pars.corsteps > 0 then -- no corrector for zero pars.corsteps
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
                maxit : in integer32;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_float;
                step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true ) is

    use Standard_Predictor_Convolutions;

    info,nbrit : integer32 := 0;
    mixres : double_float;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
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
    if pars.corsteps > 0 then -- no corrector for zero pars.corsteps
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
                maxit : in integer32;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_double;
                step : out double_double;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean ) is

    use DoblDobl_Predictor_Convolutions;

    info,nbrit : integer32 := 0;
    mixres : double_double;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
    SVD_Prediction(hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm);
    if pars.corsteps > 0 then -- no corrector for zero pars.corsteps
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
                maxit : in integer32;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_double;
                step : out double_double;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true ) is

    use DoblDobl_Predictor_Convolutions;

    info,nbrit : integer32 := 0;
    mixres : double_double;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
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
    if pars.corsteps > 0 then -- no corrector for zero pars.corsteps
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
                maxit : in integer32;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out quad_double;
                step : out quad_double;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean ) is

    use QuadDobl_Predictor_Convolutions;

    info,nbrit : integer32 := 0;
    mixres : quad_double;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
    SVD_Prediction(hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm);
    if pars.corsteps > 0 then -- no corrector for zero pars.corsteps
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
                maxit : in integer32;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out quad_double;
                step : out quad_double;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true ) is

    use QuadDobl_Predictor_Convolutions;

    info,nbrit : integer32 := 0;
    mixres : quad_double;

  begin
    Set_Lead_Coefficients(prd,psv.sol);
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
    if pars.corsteps > 0 then -- no corrector for zero pars.corsteps
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
                acct : in out double_float;
                nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                fail : out boolean ) is

    endt : constant double_float := 1.0;
    step : double_float := 0.0;
    togo : double_float; 

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    for k in 1..pars.maxsteps loop
      Predictor_Corrector_Loop(hom,abh,homlead,abhlead,pars,maxit,prd,
        psv,svh,dx,ipvt,endt,acct,step,nbpole,nbhess,nbmaxm,fail);
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
                acct : in out double_float;
                nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                fail : out boolean; verbose : in boolean := true ) is

    endt : constant double_float := 1.0;
    step : double_float := 0.0;
    togo : double_float; 

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    for k in 1..pars.maxsteps loop
      if verbose
       then put(file,"t :"); put(file,acct,3); put_line(file," :");
      end if;
      Predictor_Corrector_Loop(file,hom,abh,homlead,abhlead,pars,maxit,prd,
        psv,svh,dx,ipvt,endt,acct,step,nbpole,nbhess,nbmaxm,fail,verbose);
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
                acct : in out double_double;
                nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                fail : out boolean ) is

    endt : constant double_float := 1.0;
    step : double_double := create(0.0);
    togo : double_double; 

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    for k in 1..pars.maxsteps loop
      Predictor_Corrector_Loop(hom,abh,homlead,abhlead,pars,maxit,prd,
        psv,svh,dx,ipvt,endt,acct,step,nbpole,nbhess,nbmaxm,fail);
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
                acct : in out double_double;
                nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                fail : out boolean; verbose : in boolean := true ) is

    endt : constant double_float := 1.0;
    step : double_double := create(0.0);
    togo : double_double; 

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    for k in 1..pars.maxsteps loop
      if verbose
       then put(file,"t : "); put(file,acct,3); put_line(file," :");
      end if;
      Predictor_Corrector_Loop(file,hom,abh,homlead,abhlead,pars,maxit,prd,
        psv,svh,dx,ipvt,endt,acct,step,nbpole,nbhess,nbmaxm,fail,verbose);
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
                acct : in out quad_double;
                nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                fail : out boolean ) is

    endt : constant double_float := 1.0;
    step : quad_double := create(0.0);
    togo : quad_double; 

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    for k in 1..pars.maxsteps loop
      Predictor_Corrector_Loop(hom,abh,homlead,abhlead,pars,maxit,prd,
        psv,svh,dx,ipvt,endt,acct,step,nbpole,nbhess,nbmaxm,fail);
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
                acct : in out quad_double;
                nbpole,nbhess,nbmaxm,nbsteps : out natural32;
                fail : out boolean; verbose : in boolean := true ) is

    endt : constant double_float := 1.0;
    step : quad_double := create(0.0);
    togo : quad_double; 

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0; nbsteps := pars.maxsteps;
    for k in 1..pars.maxsteps loop
      if verbose
       then put(file,"t : "); put(file,acct,3); put_line(file," :");
      end if;
      Predictor_Corrector_Loop(file,hom,abh,homlead,abhlead,pars,maxit,prd,
        psv,svh,dx,ipvt,endt,acct,step,nbpole,nbhess,nbmaxm,fail,verbose);
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

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxit : constant integer32 := (numdeg + dendeg + 2)/2;
    ls : Link_to_Solution := Head_Of(sols);
    prd : Predictor := Create(ls.v,hom.neq,hom.deg,numdeg,dendeg,SVD);
    psv : Predictor_Vectors(hom.dim,hom.neq);
    svh : Link_to_SVD_Hessians := Create(hom.dim);
    solsptr : Solution_List := sols;
    nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    fail : boolean;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : Standard_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : Standard_Complex_VecVecs.Link_to_VecVec;
    wrk : Standard_Complex_Vectors.Link_to_Vector
        := Standard_Speelpenning_Convolutions.Allocate_Coefficients(hom.deg);
    homcff : Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
    acct : double_float;

  begin
    Allocate_Coefficients(hom.crc,homcff);
    Store_Coefficients(hom.crc,homcff);
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr); psv.sol := ls.v; acct := 0.0;
      Track_One_Path(file,hom,abh,homlead,abhlead,pars,maxit,prd,psv,svh,
        dx,ipvt,wrk,acct,nbpole,nbhess,nbmaxm,nbsteps,fail,verbose);
      ls.v := psv.sol;
      ls.t := Standard_Complex_Numbers.Create(acct); Set_Head(solsptr,ls);
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      Restore_Leading_Coefficients(abhlead,abh.crc);
      Restore_Coefficients(homcff,hom.crc);
    end loop;
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

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxit : constant integer32 := (numdeg + dendeg + 2)/2;
    ls : Link_to_Solution := Head_Of(sols);
    prd : Predictor := Create(ls.v,hom.neq,hom.deg,numdeg,dendeg,SVD);
    psv : Predictor_Vectors(hom.dim,hom.neq);
    svh : Link_to_SVD_Hessians := Create(hom.dim);
    solsptr : Solution_List := sols;
    nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    fail : boolean;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : DoblDobl_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    wrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := DoblDobl_Speelpenning_Convolutions.Allocate_Coefficients(hom.deg);
    homcff : DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
    acct : double_double;

  begin
    Allocate_Coefficients(hom.crc,homcff);
    Store_Coefficients(hom.crc,homcff);
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr); psv.sol := ls.v; acct := create(0.0);
      Track_One_Path(file,hom,abh,homlead,abhlead,pars,maxit,prd,psv,svh,
        dx,ipvt,wrk,acct,nbpole,nbhess,nbmaxm,nbsteps,fail,verbose);
      ls.v := psv.sol;
      ls.t := DoblDobl_Complex_Numbers.Create(acct); Set_Head(solsptr,ls);
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      Restore_Leading_Coefficients(abhlead,abh.crc);
      Restore_Coefficients(homcff,hom.crc);
    end loop;
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

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxit : constant integer32 := (numdeg + dendeg + 2)/2;
    ls : Link_to_Solution := Head_Of(sols);
    prd : Predictor := Create(ls.v,hom.neq,hom.deg,numdeg,dendeg,SVD);
    psv : Predictor_Vectors(hom.dim,hom.neq);
    svh : Link_to_SVD_Hessians := Create(hom.dim);
    solsptr : Solution_List := sols;
    nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    fail : boolean;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : QuadDobl_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    wrk : QuadDobl_Complex_Vectors.Link_to_Vector
        := QuadDobl_Speelpenning_Convolutions.Allocate_Coefficients(hom.deg);
    homcff : QuadDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
    acct : quad_double;

  begin
    Allocate_Coefficients(hom.crc,homcff);
    Store_Coefficients(hom.crc,homcff);
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr); psv.sol := ls.v; acct := Create(0.0);
      Track_One_Path(file,hom,abh,homlead,abhlead,pars,maxit,prd,psv,svh,
        dx,ipvt,wrk,acct,nbpole,nbhess,nbmaxm,nbsteps,fail,verbose);
      ls.v := psv.sol;
      ls.t := QuadDobl_Complex_Numbers.Create(acct); Set_Head(solsptr,ls);
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      Restore_Leading_Coefficients(abhlead,abh.crc);
      Restore_Coefficients(homcff,hom.crc);
    end loop;
    Clear(prd); Clear(svh);
    QuadDobl_Complex_Vectors.Clear(wrk);
    QuadDobl_Complex_VecVecs.Deep_Clear(homlead);
    QuadDobl_Complex_VecVecs.Deep_Clear(abhlead);
    QuadDobl_Speelpenning_Convolutions.Clear(homcff);
  end Track_All_Paths;

end Predictor_Corrector_Loops;
