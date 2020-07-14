with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Inlined_Newton_Circuits;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with Residual_Convolution_Circuits;      use Residual_Convolution_Circuits;
with Corrector_Convolutions;             use Corrector_Convolutions;
with Hyperplane_Convolution_Scaling;

package body Predictor_Corrector_Loops is

-- ON COEFFICIENT CONVOLUTION CIRCUITS :

  procedure Predictor_Corrector_Loop
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_float;
                step,mixres : out double_float; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    initres,res,err : double_float;
    info : integer32;
    nbr : natural32;
    xtr : constant natural32 := 1; -- one extra Newton step

    use Standard_Predictor_Convolutions;

  begin
    if vrblvl > 0 then
      put("-> in predictor_corrector_loops.");
      put_line("Predictor_Corrector_Loop 1 ...");
    end if;
    Standard_Predictor_Convolutions.Set_Lead_Coefficients(prd,psv.sol);
    case prd.kind is
      when LU =>
        Standard_Predictor_Convolutions.LU_Prediction
          (hom,cfh,prd.ludata,svh,rx,ix,xr,xi,vh,svls,psv,maxit,
           pars.tolres,pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,
           pars.minsize,endt,acct,fail,step,nbpole,nbhess,nbmaxm,vrblvl-1);
      when SVD =>
        Standard_Predictor_Convolutions.SVD_Prediction
          (hom,cfh,prd.svdata,svh,rx,ix,xr,xi,vh,svls,psv,maxit,
           pars.tolres,pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,
           pars.minsize,endt,acct,fail,step,nbpole,nbhess,nbmaxm,vrblvl-1);
    end case;
    if pars.corsteps = 0 then -- no corrector for zero pars.corsteps
      mixres := 1.0; nbrit := 0;
    else
      nbrit := 0;
      loop
        Standard_Predictor_Convolutions.EvalCffRad(hom,cfh,abh,step);
        Standard_Inlined_Newton_Circuits.LU_Newton_Steps
          (cfh,abh,psv.sol,psv.radsol,xr,xi,pars.corsteps,pars.tolres,
           pars.tolres,ipvt,info,initres,res,err,mixres,nbr,fail,xtr);
        nbrit := nbrit + integer32(nbr);
        exit when not fail;
        step := step/2.0;
        exit when (step < pars.minsize);
        case prd.kind is
          when LU =>
            Standard_Rational_Approximations.Evaluate
              (prd.ludata.numcff,prd.ludata.dencff,step,psv.sol);
          when SVD =>
            Standard_Rational_Approximations.Evaluate
              (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
        end case;
      end loop;
    end if;
  end Predictor_Corrector_Loop;

  procedure Predictor_Corrector_Loop
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_float;
                step,mixres : out double_float; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    initres,res,err : double_float;
    info : integer32;
    nbr : natural32;
    xtr : constant natural32 := 1; -- one extra Newton step

    use Standard_Predictor_Convolutions;

  begin
    if vrblvl > 0 then
      put("-> in predictor_corrector_loops.");
      put_line("Predictor_Corrector_Loop 2 ...");
    end if;
    Standard_Predictor_Convolutions.Set_Lead_Coefficients(prd,psv.sol);
    case prd.kind is
      when LU =>
        Standard_Predictor_Convolutions.LU_Prediction
          (file,hom,cfh,prd.ludata,svh,rx,ix,xr,xi,vh,svls,psv,maxit,
           pars.tolres,pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,
           pars.minsize,endt,acct,fail,step,nbpole,nbhess,nbmaxm,
           verbose,verbose,vrblvl-1); -- false,verbose);
      when SVD =>
        Standard_Predictor_Convolutions.SVD_Prediction
          (file,hom,cfh,prd.svdata,svh,rx,ix,xr,xi,vh,svls,psv,maxit,
           pars.tolres,pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,
           pars.minsize,endt,acct,fail,step,nbpole,nbhess,nbmaxm,
           verbose,verbose,vrblvl-1); -- false,verbose);
    end case;
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
      nbrit := 0;
      loop
        Standard_Predictor_Convolutions.EvalCffRad(hom,cfh,abh,step);
        Standard_Inlined_Newton_Circuits.LU_Newton_Steps
          (file,cfh,abh,psv.sol,psv.radsol,xr,xi,pars.corsteps,pars.tolres,
           pars.tolres,ipvt,info,initres,res,err,mixres,nbr,fail,xtr,verbose);
        nbrit := nbrit + integer32(nbr);
        exit when not fail;
        step := step/2.0;
        exit when (step < pars.minsize);
        if verbose then
          put(file,"Reduced step size to"); put(file,step,3);
          put_line(file,".");
        end if;
        case prd.kind is
          when LU =>
            Standard_Rational_Approximations.Evaluate
              (prd.ludata.numcff,prd.ludata.dencff,step,psv.sol);
          when SVD =>
            Standard_Rational_Approximations.Evaluate
              (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
        end case;
      end loop;
    end if;
  end Predictor_Corrector_Loop;

-- ON COMPLEX CONVOLUTION CIRCUITS :

  procedure Predictor_Corrector_Loop
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_float;
                step,mixres : out double_float; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use Standard_Predictor_Convolutions;

    info : integer32 := 0;
    maxdx : double_float;
    xtr : constant integer32 := 1; -- one extra Newton step

  begin
    if vrblvl > 0 then
      put("-> in predictor_corrector_loops.");
      put_line("Predictor_Corrector_Loop 3 ...");
    end if;
    Set_Lead_Coefficients(prd,psv.sol);
    if mhom = 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    elsif mhom > 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol,idz,mhom);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radii(hom,abh,mhom);
    end if;
    SVD_Prediction(hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm,vrblvl-1);
    if pars.corsteps = 0 then -- no corrector for zero pars.corsteps
      mixres := 1.0; nbrit := 0;
    else
      loop
        Store_Leading_Coefficients(hom.crc,homlead);
        Store_Leading_Coefficients(abh.crc,abhlead);
        Step_Coefficient(hom,step);
        Update_Radii_of_Constants(abh,hom);
        LU_Newton_Steps(hom,abh,psv,integer32(pars.corsteps),nbrit,
                        pars.tolres,maxdx,mixres,dx,ipvt,info,fail,xtr);
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
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_float;
                step,mixres : out double_float; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    use Standard_Predictor_Convolutions;

    info : integer32 := 0;
    maxdx : double_float;
    xtr : constant integer32 := 1; -- one extra Newton step

  begin
    if vrblvl > 0 then
      put("-> in predictor_corrector_loops.");
      put_line("Predictor_Corrector_Loop 4 ...");
    end if;
    Set_Lead_Coefficients(prd,psv.sol);
    if mhom = 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    elsif mhom > 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol,idz,mhom);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radii(hom,abh,mhom);
    end if;
    SVD_Prediction(file,hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,
      endt,acct,fail,step,nbpole,nbhess,nbmaxm,false,verbose,vrblvl-1);
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
        LU_Newton_Steps
          (file,hom,abh,psv,integer32(pars.corsteps),nbrit,
           pars.tolres,maxdx,mixres,dx,ipvt,info,fail,xtr,verbose=>verbose);
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
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_double;
                step,mixres : out double_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use DoblDobl_Predictor_Convolutions;

    info : integer32 := 0;
    maxdx : double_double;
    xtr : constant integer32 := 1; -- one extra Newton step

  begin
    if vrblvl > 0 then
      put("-> in predictor_corrector_loops.");
      put_line("Predictor_Corrector_Loop 5 ...");
    end if;
    Set_Lead_Coefficients(prd,psv.sol);
    if mhom = 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    elsif mhom > 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol,idz,mhom);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radii(hom,abh,mhom);
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
                        pars.tolres,maxdx,mixres,dx,ipvt,info,fail,xtr);
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
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out double_double;
                step,mixres : out double_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Predictor_Convolutions;

    info : integer32 := 0;
    maxdx : double_double;
    xtr : constant integer32 := 1; -- one extra Newton step

  begin
    if vrblvl > 0 then
      put("-> in predictor_corrector_loops.");
      put_line("Predictor_Corrector_Loop 6 ...");
    end if;
    Set_Lead_Coefficients(prd,psv.sol);
    if mhom = 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    elsif mhom > 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol,idz,mhom);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radii(hom,abh,mhom);
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
        LU_Newton_Steps
          (file,hom,abh,psv,integer32(pars.corsteps),nbrit,
           pars.tolres,maxdx,mixres,dx,ipvt,info,fail,xtr,verbose=>verbose);
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
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out quad_double;
                step,mixres : out quad_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use QuadDobl_Predictor_Convolutions;

    info : integer32 := 0;
    maxdx : quad_double;
    xtr : constant integer32 := 1; -- one extra Newton step

  begin
    if vrblvl > 0 then
      put("-> in predictor_corrector_loops.");
      put_line("Predictor_Corrector_Loop 7 ...");
    end if;
    Set_Lead_Coefficients(prd,psv.sol);
    if mhom = 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    elsif mhom > 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol,idz,mhom);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radii(hom,abh,mhom);
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
                        pars.tolres,maxdx,mixres,dx,ipvt,info,fail,xtr);
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
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                endt : in double_float; acct : in out quad_double;
                step,mixres : out quad_double; nbrit : out integer32;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Predictor_Convolutions;

    info : integer32 := 0;
    maxdx : quad_double;
    xtr : constant integer32 := 1; -- one extra Newton step

  begin
    if vrblvl > 0 then
      put("-> in predictor_corrector_loops.");
      put_line("Predictor_Corrector_Loop 8 ...");
    end if;
    Set_Lead_Coefficients(prd,psv.sol);
    if mhom = 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radius(hom,abh);
    elsif mhom > 1 then
      Hyperplane_Convolution_Scaling.Scale_and_Adjust(hom,psv.sol,idz,mhom);
      Hyperplane_Convolution_Scaling.Adjust_Last_Radii(hom,abh,mhom);
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
        LU_Newton_Steps
          (file,hom,abh,psv,integer32(pars.corsteps),nbrit,
           pars.tolres,maxdx,mixres,dx,ipvt,info,fail,xtr,verbose=>verbose);
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

end Predictor_Corrector_Loops;
