with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Homotopy;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_CSeries_Vector_Functions;
with Homotopy_Pade_Approximants;
with Singular_Values_of_Hessians;
with Homotopy_Mixed_Residuals;
with Homotopy_Newton_Steps;
with Series_and_Homotopies;
with Series_and_Predictors;
with Standard_Pade_Trackers;
with Homotopy_Coefficient_Scaling;

package body DoblDobl_Pade_Trackers is

  function Residual_Prediction
              ( sol : DoblDobl_Complex_Vectors.Vector;
                t : double_float ) return double_float is

    res : double_double;
    ddt : constant double_double := create(t);
    cmplxt : constant DoblDobl_Complex_Numbers.Complex_Number
           := DoblDobl_Complex_Numbers.Create(ddt);
    val : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_Homotopy.Eval(sol,cmplxt);

  begin
    res := DoblDobl_Complex_Vector_Norms.Max_Norm(val);
    return hi_part(res);
  end Residual_Prediction;

  function Residual_Prediction
              ( abh : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                sol : DoblDobl_Complex_Vectors.Vector;
                t : double_float ) return double_float is

    res : double_double;
    ddt : constant double_double := create(t);
    cmplxt : constant DoblDobl_Complex_Numbers.Complex_Number
           := DoblDobl_Complex_Numbers.Create(ddt);

  begin
    res := Homotopy_mixed_Residuals.Residual(abh,sol,cmplxt);
    return hi_part(res);
  end Residual_Prediction;

  function Residual_Prediction
              ( file : file_type;
                abh : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                sol : DoblDobl_Complex_Vectors.Vector;
                t : double_float ) return double_float is

    res : double_double;
    ddt : constant double_double := create(t);
    cmplxt : constant DoblDobl_Complex_Numbers.Complex_Number
           := DoblDobl_Complex_Numbers.Create(ddt);

  begin
    res := Homotopy_mixed_Residuals.Residual(file,abh,sol,cmplxt);
    return hi_part(res);
  end Residual_Prediction;

  procedure Predictor_Feedback
              ( abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector;
                sol : out DoblDobl_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep : in double_float;
                cntcut : in out natural32 ) is

    dd_step : double_double;

  begin
    loop
      dd_step := create(step);
      sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
      predres := Residual_Prediction(abh,sol,t);
      exit when (predres <= tolpres);
      t := t - step; step := step/2.0; t := t + step;
      cntcut := cntcut + 1;
      exit when (step <= minstep);
    end loop;
  end Predictor_Feedback;

  procedure Predictor_Feedback
              ( file : in file_type;
                abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector;
                sol : out DoblDobl_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep : in double_float;
                cntcut : in out natural32 ) is

    dd_step : double_double;

  begin
    loop
      dd_step := create(step);
      sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
      predres := Residual_Prediction(abh,sol,t);
      put(file,"  predictor residual : ");
      put(file,predres,3); new_line(file);
      exit when (predres <= tolpres);
      t := t - step; step := step/2.0; t := t + step;
      cntcut := cntcut + 1;
      put(file,"Cut step size : "); put(file,step,3);
      put(file," t = "); put(file,t,3);
      exit when (step < minstep);
    end loop;
  end Predictor_Feedback;

  procedure Predictor_Feedback
              ( file : in file_type; verbose : in boolean;
                abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector;
                sol : out DoblDobl_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep : in double_float;
                cntcut : in out natural32 ) is
  begin
    if verbose then
      Predictor_Feedback
        (file,abh,pv,sol,predres,t,step,tolpres,minstep,cntcut);
    else
      Predictor_Feedback(abh,pv,sol,predres,t,step,tolpres,minstep,cntcut);
    end if;
  end Predictor_Feedback;

  procedure Predictor_Corrector
              ( abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector;
                sol : out DoblDobl_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep,tolcres : in double_float;
                maxit,extra : in natural32; nbrcorrs : in out natural32;
                err,rco,res : out double_float;
                cntcut,cntfail : in out natural32; fail : out boolean ) is

    nbrit : natural32 := 0;

  begin
    fail := true;
    loop
      Predictor_Feedback(abh,pv,sol,predres,t,step,tolpres,tolpres,cntcut);
     -- exit when ((step < pars.minsize) and (predres > alpha));
      Homotopy_Newton_Steps.Correct
        (abh,t,tolcres,maxit,nbrit,sol,err,rco,res,fail,extra);
      nbrcorrs := nbrcorrs + nbrit;
      if fail
       then fail := false; -- ignore corrector failure
      end if;
      exit when (not fail);
      t := t - step; step := step/2.0; t := t + step;
      cntfail := cntfail + 1;
      exit when (step < minstep);
    end loop;
  end Predictor_Corrector;

  procedure Predictor_Corrector
              ( file : in file_type; verbose : in boolean;
                abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector;
                sol : out DoblDobl_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep,tolcres : in double_float;
                maxit,extra : in natural32; nbrcorrs : in out natural32;
                err,rco,res : out double_float;
                cntcut,cntfail : in out natural32; fail : out boolean ) is

    nbrit : natural32 := 0;

  begin
    fail := true;
    loop
      Predictor_Feedback
        (file,verbose,abh,pv,sol,predres,t,step,tolpres,tolpres,cntcut);
     -- exit when ((step < pars.minsize) and (predres > alpha));
      Homotopy_Newton_Steps.Correct
        (file,abh,t,tolcres,maxit,nbrit,sol,err,rco,res,fail,extra,verbose);
      nbrcorrs := nbrcorrs + nbrit;
      if verbose then
        if fail then
          put_line(file,"Warning: The correct stage failed, will ignore.");
          fail := false;
        else
          put_line(file,"The correct stage succeeded.");
        end if;
      end if;
      exit when (not fail);
      t := t - step; step := step/2.0; t := t + step;
      cntfail := cntfail + 1;
      exit when (step < minstep);
    end loop;
  end Predictor_Corrector;

  procedure Step_Control
              ( jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out DoblDobl_Pade_Approximants.Pade_Vector;
                poles : in out DoblDobl_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 ) is

    srv : DoblDobl_Complex_Series_Vectors.Vector(sol'range);
    eva : DoblDobl_Complex_Series_Vectors.Vector(hom'range);
    frp : double_double;
    cfp : DoblDobl_Complex_Numbers.Complex_Number;
    sstep,dstep,pstep : double_float;
    dd_t : double_double;
    onetarget : constant double_float := 1.0;
   -- alpha : constant double_float := pars.alpha;
   -- tolcff : constant double_float := pars.epsilon;

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Step_Control 1 ...");
    end if;
    Series_and_Predictors.Newton_Prediction
      (maxdeg,nit,hom,sol,srv,eva,vrblvl=>vrblvl-1);
   -- sstep := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
    sstep := 1.0; -- disable series step -- pars.sbeta*sstep;
    Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
    dd_t := Create(t);
    dstep := Series_and_Predictors.Step_Distance
               (maxdeg,pars.cbeta,dd_t,jm,hs,sol,srv,pv);
    pstep := pars.pbeta*hi_part(frp);
    Standard_Pade_Trackers.Minimum_Step_Size
      (sstep,dstep,pstep,step,cntsstp,cntdstp,cntpstp);
    Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
    DoblDobl_Complex_Series_Vectors.Clear(eva);
    DoblDobl_Complex_Series_Vectors.Clear(srv);
  end Step_Control;

  procedure Step_Control
              ( file : in file_type; verbose : in boolean;
                jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out DoblDobl_Pade_Approximants.Pade_Vector;
                poles : in out DoblDobl_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 ) is

    srv : DoblDobl_Complex_Series_Vectors.Vector(sol'range);
    eva : DoblDobl_Complex_Series_Vectors.Vector(hom'range);
    frp,eta,nrm : double_double;
    cfp : DoblDobl_Complex_Numbers.Complex_Number;
    sstep,dstep,pstep : double_float;
    dd_t : double_double;
    onetarget : constant double_float := 1.0;
   -- alpha : constant double_float := pars.alpha;
   -- tolcff : constant double_float := pars.epsilon;

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Step_Control 2 ...");
    end if;
    Series_and_Predictors.Newton_Prediction
      (file,maxdeg,nit,hom,sol,srv,eva,false,vrblvl=>vrblvl-1); -- verbose);
   -- sstep := Series_and_Predictors.Set_Step_Size
   --            (file,eva,tolcff,alpha,verbose);
    sstep := 1.0; -- disable series step -- pars.sbeta*sstep;
   -- if verbose
   --  then put(file,"series step : "); put(file,sstep,2); new_line(file);
   -- end if;
    Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
    pstep := pars.pbeta*hi_part(frp);
    if verbose then
      put(file,"pole step : "); put(file,pstep,2);
      put(file,"  smallest pole radius : "); put(file,frp,2); new_line(file);
      put(file,"Closest pole : "); put(file,cfp); new_line(file);
    end if;
    dd_t := Create(t);
    if not verbose then
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,pars.cbeta,dd_t,jm,hs,sol,srv,pv);
    else
      declare -- must extend the solution vector with the value of t
        solxt : DoblDobl_Complex_Vectors.Vector(sol'first..sol'last+1);
        use Singular_Values_of_Hessians;
      begin
        solxt(sol'range) := sol;
        solxt(solxt'last) := DoblDobl_Complex_Numbers.Create(dd_t);
        eta := DoblDobl_Distance(jm.all,hs.all,solxt);
      end;
      nrm := Homotopy_Pade_Approximants.Solution_Error_Norm(srv,pv);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,pars.cbeta,hi_part(eta),hi_part(nrm));
      put(file,"Hessian step : "); put(file,dstep,2);
      put(file,"  eta : "); put(file,eta,2);
      put(file,"  nrm : "); put(file,nrm,2); new_line(file);
    end if;
    Standard_Pade_Trackers.Minimum_Step_Size
      (file,sstep,dstep,pstep,step,cntsstp,cntdstp,cntpstp);
    Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
    if verbose then
      put(file,"Step size : "); put(file,step,3);
      put(file," t = "); put(file,t,3);
    end if;
    DoblDobl_Complex_Series_Vectors.Clear(eva);
    DoblDobl_Complex_Series_Vectors.Clear(srv);
  end Step_Control;

  procedure Step_Control
              ( jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in DoblDobl_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out DoblDobl_Pade_Approximants.Pade_Vector;
                poles : in out DoblDobl_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 ) is

    srv : DoblDobl_Complex_Series_Vectors.Vector(sol'range);
    eva : DoblDobl_Complex_Series_Vectors.Vector(fhm'range);
    frp : double_double;
    cfp : DoblDobl_Complex_Numbers.Complex_Number;
    sstep,dstep,pstep : double_float;
    dd_t : double_double;
    onetarget : constant double_float := 1.0;
   -- alpha : constant double_float := pars.alpha;
   -- tolcff : constant double_float := pars.epsilon;

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Step_Control 3 ...");
    end if;
    Series_and_Predictors.Newton_Prediction
      (maxdeg,nit,fhm,fcf,ejm,mlt,sol,srv,eva,vrblvl=>vrblvl-1);
   -- sstep := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
    sstep := 1.0; -- disable series step -- pars.sbeta*sstep;
    Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
    pstep := pars.pbeta*hi_part(frp);
    dd_t := Create(t);
    dstep := Series_and_Predictors.Step_Distance
                (maxdeg,pars.cbeta,dd_t,jm,hs,sol,srv,pv);
    Standard_Pade_Trackers.Minimum_Step_Size
      (sstep,dstep,pstep,step,cntsstp,cntdstp,cntpstp);
    Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
    DoblDobl_Complex_Series_Vectors.Clear(eva);
    DoblDobl_Complex_Series_Vectors.Clear(srv);
  end Step_Control;

  procedure Step_Control
              ( file : in file_type; verbose : in boolean;
                jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in DoblDobl_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out DoblDobl_Pade_Approximants.Pade_Vector;
                poles : in out DoblDobl_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 ) is

    srv : DoblDobl_Complex_Series_Vectors.Vector(sol'range);
    eva : DoblDobl_Complex_Series_Vectors.Vector(fhm'range);
    frp,eta,nrm : double_double;
    cfp : DoblDobl_Complex_Numbers.Complex_Number;
    sstep,dstep,pstep : double_float;
    dd_t : double_double;
    onetarget : constant double_float := 1.0;
   -- alpha : constant double_float := pars.alpha;
   -- tolcff : constant double_float := pars.epsilon;

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Step_Control 4 ...");
    end if;
    Series_and_Predictors.Newton_Prediction -- verbose flag set to false
      (file,maxdeg,nit,fhm,fcf,ejm,mlt,sol,srv,eva,false,vrblvl=>vrblvl-1);
   -- sstep := Series_and_Predictors.Set_Step_Size
   --            (file,eva,tolcff,alpha,verbose);
    sstep := 1.0; -- disable series step -- pars.sbeta*sstep;
   -- if verbose
   --  then put(file,"series step : "); put(file,step,2); new_line(file);
   -- end if;
    Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
    pstep := pars.pbeta*hi_part(frp);
    if verbose then
      put(file,"pole step : "); put(file,pstep,2);
      put(file,"  smallest pole radius : "); put(file,frp,2); new_line(file);
      put(file,"closest pole : "); put(file,cfp); new_line(file);
    end if;
    dd_t := Create(t);
    if not verbose then
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,pars.cbeta,dd_t,jm,hs,sol,srv,pv);
    else
      declare -- must extend the solution vector with the value of t
        solxt : DoblDobl_Complex_Vectors.Vector(sol'first..sol'last+1);
        use Singular_Values_of_Hessians;
      begin
        solxt(sol'range) := sol;
        solxt(solxt'last) := DoblDobl_Complex_Numbers.Create(dd_t);
        eta := DoblDobl_Distance(jm.all,hs.all,solxt);
      end;
      nrm := Homotopy_Pade_Approximants.Solution_Error_Norm(srv,pv);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,pars.cbeta,hi_part(eta),hi_part(nrm));
      put(file,"Hessian step : "); put(file,dstep,2);
      put(file,"  eta : "); put(file,eta,2);
      put(file,"  nrm : "); put(file,nrm,2); new_line(file);
    end if;
    Standard_Pade_Trackers.Minimum_Step_Size
      (file,sstep,dstep,pstep,step,cntsstp,cntdstp,cntpstp);
    Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
    if verbose then
      put(file,"Step size : "); put(file,step,3);
      put(file," t = "); put(file,t,3);
    end if;
    DoblDobl_Complex_Series_Vectors.Clear(eva);
    DoblDobl_Complex_Series_Vectors.Clear(srv);
  end Step_Control;

  procedure Track_One_Path
              ( abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntcut,cntfail : out natural32;
                minsize,maxsize : out double_float;
                cntsstp,cntdstp,cntpstp : out natural32;
                vrblvl : in integer32 := 0 ) is

    wrk : DoblDobl_CSeries_Poly_Systems.Poly_Sys(hom'range);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := Standard_Pade_Trackers.Maximum(5,maxdeg/2);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(1..sol.n)
       := DoblDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_DoblDobl_Poles(sol.n,dendeg);
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    maxit : constant natural32 := pars.corsteps;
    extra : constant natural32 := 1;
    fail : boolean;
    t,step : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    err,rco,res : double_float;
    predres : double_float;
    nbrit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Track_One_Path 1 ...");
    end if;
    minsize := 1.0; maxsize := 0.0; cntsstp := 0; cntdstp := 0; cntpstp := 0;
    DoblDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0; cntcut := 0; cntfail := 0; nbrsteps := max_steps;
    for k in 1..max_steps loop
      Step_Control(jm,hs,wrk,wrk_sol,maxdeg,nit,pars,pv,poles,t,step,
                   cntsstp,cntdstp,cntpstp,vrblvl-1);
      Predictor_Corrector
        (abh,pv,wrk_sol,predres,t,step,alpha,pars.minsize,tolres,
         maxit,extra,nbrcorrs,err,rco,res,cntcut,cntfail,fail);
      Standard_Pade_Trackers.Update_Step_Sizes(minsize,maxsize,step);
      dd_step := create(step);
      Series_and_Homotopies.Shift(wrk,-dd_step);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
    end loop;
    DoblDobl_Pade_Approximants.Clear(pv);
    DoblDobl_Complex_VecVecs.Clear(poles);
    DoblDobl_CSeries_Poly_Systems.Clear(wrk);
    dd_t := create(-1.0);
    wrk := Series_and_Homotopies.Shift(hom,dd_t);
    Homotopy_Newton_Steps.Correct
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,extra);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := DoblDobl_Complex_Numbers.Create(Double_Double_Numbers.Create(t));
    sol.v := wrk_sol;
    sol.err := Double_Double_Numbers.create(err);
    sol.rco := Double_Double_Numbers.create(rco);
    sol.res := Double_Double_Numbers.create(res);
    DoblDobl_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntcut,cntfail : out natural32;
                minsize,maxsize : out double_float;
                cntsstp,cntdstp,cntpstp : out natural32;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    wrk : DoblDobl_CSeries_Poly_Systems.Poly_Sys(hom'range);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := Standard_Pade_Trackers.Maximum(5,maxdeg/2);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(1..sol.n)
       := DoblDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_DoblDobl_Poles(sol.n,dendeg);
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    maxit : constant natural32 := pars.corsteps;
    extra : constant natural32 := 1;
    fail : boolean;
    t,step : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    err,rco,res : double_float;
    predres : double_float;
    nbrit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Track_One_Path 2 ...");
    end if;
    minsize := 1.0; maxsize := 0.0; cntsstp := 0; cntdstp := 0; cntpstp := 0;
    DoblDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0; cntcut := 0; cntfail := 0; nbrsteps := max_steps;
    for k in 1..max_steps loop
      if verbose then
        put(file,"Step "); put(file,k,1); put_line(file," : ");
      end if;
      Step_Control
        (file,verbose,jm,hs,wrk,wrk_sol,maxdeg,nit,pars,pv,poles,t,step,
         cntsstp,cntdstp,cntpstp,vrblvl-1);
      Predictor_Corrector
        (file,verbose,abh,pv,wrk_sol,predres,t,step,alpha,pars.minsize,tolres,
         maxit,extra,nbrcorrs,err,rco,res,cntcut,cntfail,fail);
      Standard_Pade_Trackers.Update_Step_Sizes(minsize,maxsize,step);
      dd_step := create(step);
      Series_and_Homotopies.Shift(wrk,-dd_step);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
    end loop;
    DoblDobl_Pade_Approximants.Clear(pv);
    DoblDobl_Complex_VecVecs.Clear(poles);
    DoblDobl_CSeries_Poly_Systems.Clear(wrk);
    dd_t := create(-1.0);
    wrk := Series_and_Homotopies.Shift(hom,dd_t);
    Homotopy_Newton_Steps.Correct
      (file,abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,
       extra,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := DoblDobl_Complex_Numbers.Create(Double_Double_Numbers.Create(t));
    sol.v := wrk_sol;
    sol.err := Double_Double_Numbers.create(err);
    sol.rco := Double_Double_Numbers.create(rco);
    sol.res := Double_Double_Numbers.create(res);
    DoblDobl_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

-- VERSIONS ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure Track_One_Path
              ( abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                nbrsteps,nbrcorrs,cntcut,cntfail : out natural32;
                minsize,maxsize : out double_float;
                cntsstp,cntdstp,cntpstp : out natural32;
                vrblvl : in integer32 := 0 ) is

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := Standard_Pade_Trackers.Maximum(5,maxdeg/2);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(1..sol.n)
       := DoblDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_DoblDobl_Poles(sol.n,dendeg);
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    maxit : constant natural32 := pars.corsteps;
    extra : constant natural32 := 1;
    fail : boolean;
    t,step : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    err,rco,res,predres : double_float;
    nbrit : natural32 := 0;
    wrk_fcf : DoblDobl_Complex_Series_VecVecs.VecVec(fcf'range);
    dd_gamma : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(pars.gamma);

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Track_One_Path 3 ...");
    end if;
    minsize := 1.0; maxsize := 0.0; cntsstp := 0; cntdstp := 0; cntpstp := 0;
    nbrcorrs := 0; cntcut := 0; cntfail := 0;
    nbrsteps := max_steps;
    wrk_fcf := DoblDobl_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
    for k in 1..max_steps loop
      if mhom > 0 then
        dd_t := Double_Double_Numbers.Create(t);
        if mhom = 1 then
          Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
            (wrk_fcf,wrk_sol,dd_t,dd_gamma);
        elsif mhom > 1 then
          Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
            (wrk_fcf,wrk_sol,dd_t,dd_gamma,mhom,idz.all);
        end if;
      end if;
      Step_Control
        (jm,hs,fhm,wrk_fcf,ejm,mlt,wrk_sol,maxdeg,nit,pars,pv,poles,t,step,
         cntsstp,cntdstp,cntpstp,vrblvl-1);
      Predictor_Corrector
        (abh,pv,wrk_sol,predres,t,step,alpha,pars.minsize,tolres,
         maxit,extra,nbrcorrs,err,rco,res,cntcut,cntfail,fail);
      Standard_Pade_Trackers.Update_Step_Sizes(minsize,maxsize,step);
      dd_step := create(step);
      DoblDobl_CSeries_Vector_Functions.Shift(wrk_fcf,-dd_step);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
    end loop;
    DoblDobl_Pade_Approximants.Clear(pv);
    DoblDobl_Complex_VecVecs.Clear(poles);
    Homotopy_Newton_Steps.Correct
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,extra);
    nbrcorrs := nbrcorrs + nbrit;
    dd_t := Double_Double_Numbers.Create(t);
    sol.t := DoblDobl_Complex_Numbers.Create(dd_t);
    sol.v := wrk_sol;
    sol.err := Double_Double_Numbers.Create(err);
    sol.rco := Double_Double_Numbers.Create(rco);
    sol.res := Double_Double_Numbers.Create(res);
    DoblDobl_CSeries_Vector_Functions.Deep_Clear(wrk_fcf);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                nbrsteps,nbrcorrs,cntcut,cntfail : out natural32;
                minsize,maxsize : out double_float;
                cntsstp,cntdstp,cntpstp : out natural32;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := Standard_Pade_Trackers.Maximum(5,maxdeg/2);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(1..sol.n)
       := DoblDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_DoblDobl_Poles(sol.n,dendeg);
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    maxit : constant natural32 := pars.corsteps;
    extra : constant natural32 := 1;
    fail : boolean;
    t,step : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    err,rco,res,predres : double_float;
    nbrit : natural32 := 0;
    wrk_fcf : DoblDobl_Complex_Series_VecVecs.VecVec(fcf'range);
    dd_gamma : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(pars.gamma);

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Track_One_Path 4 ...");
    end if;
    minsize := 1.0; maxsize := 0.0; cntsstp := 0; cntdstp := 0; cntpstp := 0;
    nbrcorrs := 0; cntcut := 0; cntfail := 0; nbrsteps := max_steps;
    wrk_fcf := DoblDobl_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
    for k in 1..max_steps loop
      dd_t := Double_Double_Numbers.Create(t);
      if verbose then
        put(file,"Step "); put(file,k,1); put_line(file," : ");
        if vrblvl > 0 then
          if mhom = 1 then
            Homotopy_Coefficient_Scaling.Last_Coefficients
              (file,wrk_fcf(wrk_fcf'last),dd_t,dd_gamma);
          elsif mhom > 1 then
            Homotopy_Coefficient_Scaling.Last_Coefficients
              (file,wrk_fcf,dd_t,dd_gamma,mhom);
          end if;
        end if;
      end if;
      if mhom = 1 then
        if vrblvl > 0 then
          Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
            (file,fhm,wrk_fcf,wrk_sol,dd_t,dd_gamma,true);
        else
          Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
            (file,fhm,wrk_fcf,wrk_sol,dd_t,dd_gamma);
        end if;
      elsif mhom > 1 then
        if vrblvl > 0 then
          Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
            (file,fhm,wrk_fcf,wrk_sol,dd_t,dd_gamma,mhom,idz.all,true);
        else
          Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
            (file,fhm,wrk_fcf,wrk_sol,dd_t,dd_gamma,mhom,idz.all);
        end if;
      end if;
      Step_Control
        (file,verbose,jm,hs,fhm,wrk_fcf,ejm,mlt,wrk_sol,maxdeg,nit,pars,
         pv,poles,t,step,cntsstp,cntdstp,cntpstp,vrblvl-1);
      Predictor_Corrector
        (file,verbose,abh,pv,wrk_sol,predres,t,step,alpha,pars.minsize,tolres,
         maxit,extra,nbrcorrs,err,rco,res,cntcut,cntfail,fail);
      Standard_Pade_Trackers.Update_Step_Sizes(minsize,maxsize,step);
      dd_step := create(step);
      DoblDobl_CSeries_Vector_Functions.Shift(wrk_fcf,-dd_step);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
    end loop;
    DoblDobl_Pade_Approximants.Clear(pv);
    DoblDobl_Complex_VecVecs.Clear(poles);
    Homotopy_Newton_Steps.Correct
      (file,abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,
       extra,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    dd_t := Double_Double_Numbers.Create(t);
    sol.t := DoblDobl_Complex_Numbers.Create(dd_t);
    sol.v := wrk_sol;
    sol.err := Double_Double_Numbers.Create(err);
    sol.rco := Double_Double_Numbers.Create(rco);
    sol.res := Double_Double_Numbers.Create(res);
    DoblDobl_CSeries_Vector_Functions.Deep_Clear(wrk_fcf);
  end Track_One_Path;

end DoblDobl_Pade_Trackers;
