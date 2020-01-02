with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vector_Norms;
with Standard_Homotopy;
with Standard_Complex_Series_Vectors;
with Standard_CSeries_Vector_Functions;
with Homotopy_Pade_Approximants;
with Singular_Values_of_Hessians;
with Homotopy_Mixed_Residuals;
with Homotopy_Newton_Steps;
with Series_and_Homotopies;
with Series_and_Predictors;
with Homotopy_Coefficient_Scaling;

package body Standard_Pade_Trackers is

  function Maximum ( a,b : integer32 ) return integer32 is
  begin
    if a > b
     then return a;
     else return b;
    end if;
  end Maximum;

  function Minimum ( a,b : double_float ) return double_float is
  begin
    if a < b
     then return a;
     else return b;
    end if;
  end Minimum;

  procedure Set_Step
              ( t,step : in out double_float;
                maxstep,target : in double_float ) is

    update : double_float;

  begin
    if step > maxstep
     then step := maxstep;
    end if;
    update := t + step;
    if update <= target then
      t := update;
    else
      step := target - t;
      t := target; 
    end if;
  end Set_Step;

  procedure Update_Step_Sizes
              ( minsize,maxsize : in out double_float;
                step : in double_float ) is
  begin
    if step < minsize then
      minsize := step;
    elsif step > maxsize then
      maxsize := step;
    end if;
  end Update_Step_Sizes;

  function Residual_Prediction
              ( sol : Standard_Complex_Vectors.Vector;
                t : double_float ) return double_float is

    res : double_float;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    val : constant Standard_Complex_Vectors.Vector
       := Standard_Homotopy.Eval(sol,cmplxt);

  begin
    res := Standard_Complex_Vector_Norms.Max_Norm(val);
    return res;
  end Residual_Prediction;

  function Residual_Prediction
              ( abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                sol : Standard_Complex_Vectors.Vector;
                t : double_float ) return double_float is

    res : double_float;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);

  begin
    res := Homotopy_mixed_Residuals.Residual(abh,sol,cmplxt);
    return res;
  end Residual_Prediction;

  function Residual_Prediction
              ( file : file_type;
                abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                sol : Standard_Complex_Vectors.Vector;
                t : double_float ) return double_float is

    res : double_float;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);

  begin
    res := Homotopy_mixed_Residuals.Residual(file,abh,sol,cmplxt);
    return res;
  end Residual_Prediction;

  procedure Predictor_Feedback
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                sol : out Standard_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep : in double_float;
                cntcut : in out natural32 ) is
  begin
    loop
      sol := Series_and_Predictors.Predicted_Solution(pv,step);
      predres := Residual_Prediction(abh,sol,t);
      exit when (predres <= tolpres);
      t := t - step; step := step/2.0; t := t + step;
      cntcut := cntcut + 1;
      exit when (step <= minstep);
    end loop;
  end Predictor_Feedback;

  procedure Predictor_Feedback
              ( file : in file_type;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                sol : out Standard_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep : in double_float;
                cntcut : in out natural32 ) is
  begin
    loop
      sol := Series_and_Predictors.Predicted_Solution(pv,step);
      predres := Residual_Prediction(abh,sol,t);
      put(file,"  predictor residual : ");
      put(file,predres,3); new_line(file);
      exit when (predres <= tolpres);
      t := t - step; step := step/2.0; t := t + step;
      cntcut := cntcut + 1;
      put(file,"Cut step size : "); put(file,step,3);
      put(file," t = "); put(file,t,3);
      exit when (step <= minstep);
    end loop;
  end Predictor_Feedback;

  procedure Predictor_Feedback
              ( file : in file_type; verbose : in boolean;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                sol : out Standard_Complex_Vectors.Vector;
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
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                sol : out Standard_Complex_Vectors.Vector;
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
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                sol : out Standard_Complex_Vectors.Vector;
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
          put_line(file,"Warning: the correct stage failed, will ignore...");
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

  procedure Minimum_Step_Size
              ( sstp,dstp,pstp : in double_float; minstp : out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32 ) is
  begin
    if sstp < dstp then
      if sstp < pstp then
        minstp := sstp; cntsstp := cntsstp + 1;
      else -- pstp <= sstp < dstp
        minstp := pstp; cntpstp := cntpstp + 1;
      end if;
    elsif pstp < dstp then
      minstp := pstp; cntpstp := cntpstp + 1;
    else -- dstp <= pstp <= sstp
      minstp := dstp; cntdstp := cntdstp + 1;
    end if;
  end Minimum_Step_Size;

  procedure Minimum_Step_Size
              ( file : in file_type;
                sstp,dstp,pstp : in double_float; minstp : out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32 ) is
  begin
    if sstp < dstp then
      if sstp < pstp then
        minstp := sstp; cntsstp := cntsstp + 1;
        put(file,"series step is mimimal, count = ");
        put(file,cntsstp,1); new_line(file);
      else -- pstp <= sstp < dstp
        minstp := pstp; cntpstp := cntpstp + 1;
        put(file,"pole step is mimimal, count = ");
        put(file,cntpstp,1); new_line(file);
      end if;
    elsif pstp < dstp then
      minstp := pstp; cntpstp := cntpstp + 1;
      put(file,"pole step is mimimal, count = ");
      put(file,cntpstp,1); new_line(file);
    else -- dstp <= pstp <= sstp
      minstp := dstp; cntdstp := cntdstp + 1;
      put(file,"curvature step is mimimal, count = ");
      put(file,cntdstp,1); new_line(file);
    end if;
  end Minimum_Step_Size;

  procedure Step_Control
              ( jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out Standard_Pade_Approximants.Pade_Vector;
                poles : in out Standard_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 ) is

    srv : Standard_Complex_Series_Vectors.Vector(sol'range);
    eva : Standard_Complex_Series_Vectors.Vector(hom'range);
    frp : double_float;
    cfp : Standard_Complex_Numbers.Complex_Number;
    sstep,dstep,pstep : double_float;
    onetarget : constant double_float := 1.0;
   -- alpha : constant double_float := pars.alpha;
   -- tolcff : constant double_float := pars.epsilon;

  begin
    if vrblvl > 0
     then put_line("-> in standard_pade_trackers.Step_Control 1 ...");
    end if;
    Series_and_Predictors.Newton_Prediction
      (maxdeg,nit,hom,sol,srv,eva,vrblvl=>vrblvl-1);
   -- sstep := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
    sstep := 1.0; -- disable series step --  pars.sbeta*sstep;
    Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
    dstep := Series_and_Predictors.Step_Distance
               (maxdeg,pars.cbeta,t,jm,hs,sol,srv,pv);
    pstep := pars.pbeta*frp;
    Minimum_Step_Size(sstep,dstep,pstep,step,cntsstp,cntdstp,cntpstp);
    Set_Step(t,step,pars.maxsize,onetarget);
    Standard_Complex_Series_Vectors.Clear(eva);
    Standard_Complex_Series_Vectors.Clear(srv);
  end Step_Control;

  procedure Step_Control
              ( jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in Standard_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out Standard_Pade_Approximants.Pade_Vector;
                poles : in out Standard_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 ) is

    srv : Standard_Complex_Series_Vectors.Vector(sol'range);
    eva : Standard_Complex_Series_Vectors.Vector(fhm'range);
    frp : double_float;
    cfp : Standard_Complex_Numbers.Complex_Number;
    sstep,dstep,pstep : double_float;
    onetarget : constant double_float := 1.0;
   -- alpha : constant double_float := pars.alpha;
   -- tolcff : constant double_float := pars.epsilon;

  begin
    if vrblvl > 0
     then put_line("-> in standard_pade_trackers.Step_Control 2 ...");
    end if;
    Series_and_Predictors.Newton_Prediction
      (maxdeg,nit,fhm,fcf,ejm,mlt,sol,srv,eva,vrblvl=>vrblvl-1);
   -- sstep := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
    sstep := 1.0; -- disable series step -- pars.sbeta*sstep;
    Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
    dstep := Series_and_Predictors.Step_Distance
               (maxdeg,pars.cbeta,t,jm,hs,sol,srv,pv);
    pstep := pars.pbeta*frp;
    Minimum_Step_Size(sstep,dstep,pstep,step,cntsstp,cntdstp,cntpstp);
    Set_Step(t,step,pars.maxsize,onetarget);
    Standard_Complex_Series_Vectors.Clear(eva);
    Standard_Complex_Series_Vectors.Clear(srv);
  end Step_Control;

  procedure Step_Control
              ( file : in file_type; verbose : in boolean;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out Standard_Pade_Approximants.Pade_Vector;
                poles : in out Standard_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 ) is

    srv : Standard_Complex_Series_Vectors.Vector(sol'range);
    eva : Standard_Complex_Series_Vectors.Vector(hom'range);
    frp,eta,nrm : double_float;
    cfp : Standard_Complex_Numbers.Complex_Number;
    sstep,dstep,pstep : double_float;
    onetarget : constant double_float := 1.0;
   -- alpha : constant double_float := pars.alpha;
   -- tolcff : constant double_float := pars.epsilon;

  begin
    if vrblvl > 0
     then put_line("-> in standard_pade_trackers.Step_Control 3 ...");
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
    pstep := pars.pbeta*frp;
    if verbose then
      put(file,"pole step : "); put(file,pstep,2); new_line(file);
      put(file,"  smallest pole radius : "); put(file,frp,2);
      put(file,"closest pole : "); put(file,cfp); new_line(file);
    end if;
    if not verbose then
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,pars.cbeta,t,jm,hs,sol,srv,pv);
    else
      declare -- must extend the solution vector with the value for t
        solxt : Standard_Complex_Vectors.Vector(sol'first..sol'last+1);
        use Singular_Values_of_Hessians;
      begin
        solxt(sol'range) := sol;
        solxt(solxt'last) := Standard_Complex_Numbers.Create(t);
        eta := Standard_Distance(jm.all,hs.all,solxt);
      end;
      nrm := Homotopy_Pade_Approximants.Solution_Error_Norm(srv,pv);
      dstep := Series_and_Predictors.Step_Distance(maxdeg,pars.cbeta,eta,nrm);
      put(file,"Hessian step : "); put(file,dstep,2);
      put(file,"  eta : "); put(file,eta,2);
      put(file,"  nrm : "); put(file,nrm,2); new_line(file);
    end if;
    if verbose then
      Minimum_Step_Size(file,sstep,dstep,pstep,step,cntsstp,cntdstp,cntpstp);
    else
      Minimum_Step_Size(sstep,dstep,pstep,step,cntsstp,cntdstp,cntpstp);
    end if;
    Set_Step(t,step,pars.maxsize,onetarget);
    if verbose then
      put(file,"Step size : "); put(file,step,3);
      put(file," t = "); put(file,t,3);
    end if;
    Standard_Complex_Series_Vectors.Clear(eva);
    Standard_Complex_Series_Vectors.Clear(srv);
  end Step_Control;

  procedure Step_Control
              ( file : in file_type; verbose : in boolean;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in Standard_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out Standard_Pade_Approximants.Pade_Vector;
                poles : in out Standard_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 ) is

    srv : Standard_Complex_Series_Vectors.Vector(sol'range);
    eva : Standard_Complex_Series_Vectors.Vector(fhm'range);
    frp,eta,nrm : double_float;
    cfp : Standard_Complex_Numbers.Complex_Number;
    sstep,dstep,pstep : double_float;
    onetarget : constant double_float := 1.0;
   -- alpha : constant double_float := pars.alpha;
   -- tolcff : constant double_float := pars.epsilon;

  begin
    if vrblvl > 0
     then put_line("-> in standard_pade_trackers.Step_Control 4 ...");
    end if;
    Series_and_Predictors.Newton_Prediction
      (file,maxdeg,nit,fhm,fcf,ejm,mlt,sol,srv,eva,false,vrblvl=>vrblvl-1);
   -- verbose);
   -- sstep := Series_and_Predictors.Set_Step_Size
   --            (file,eva,tolcff,alpha,verbose);
    sstep := 1.0; -- disable series step -- pars.sbeta*sstep;
   -- if verbose
   --  then put(file,"series step : "); put(file,sstep,2); new_line(file);
   -- end if;
    Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
    pstep := pars.pbeta*frp;
    if verbose then
      put(file,"pole step : "); put(file,pstep,2);
      put(file,"  smallest pole radius : "); put(file,frp,2); new_line(file);
      put(file,"closest pole : "); put(file,cfp); new_line(file);
    end if;
    if not verbose then
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,pars.cbeta,t,jm,hs,sol,srv,pv);
    else
      declare -- must extend the solution vector with the value for t
        solxt : Standard_Complex_Vectors.Vector(sol'first..sol'last+1);
        use Singular_Values_of_Hessians;
      begin
        solxt(sol'range) := sol;
        solxt(solxt'last) := Standard_Complex_Numbers.Create(t);
        eta := Standard_Distance(jm.all,hs.all,solxt);
      end;
      nrm := Homotopy_Pade_Approximants.Solution_Error_Norm(srv,pv);
      dstep := Series_and_Predictors.Step_Distance(maxdeg,pars.cbeta,eta,nrm);
      put(file,"Hessian step : "); put(file,dstep,2);
      put(file,"  eta : "); put(file,eta,2);
      put(file,"  nrm : "); put(file,nrm,2); new_line(file);
    end if;
    if verbose then
      Minimum_Step_Size(file,sstep,dstep,pstep,step,cntsstp,cntdstp,cntpstp);
    else
      Minimum_Step_Size(sstep,dstep,pstep,step,cntsstp,cntdstp,cntpstp);
    end if;
    Set_Step(t,step,pars.maxsize,onetarget);
    if verbose then
      put(file,"Step size : "); put(file,step,3);
      put(file,"  t = "); put(file,t,3);
    end if;
    Standard_Complex_Series_Vectors.Clear(eva);
    Standard_Complex_Series_Vectors.Clear(srv);
  end Step_Control;

  procedure Track_One_Path
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntcut,cntfail : out natural32;
                minsize,maxsize : out double_float;
                cntsstp,cntdstp,cntpstp : out natural32;
                vrblvl : in integer32 := 0 ) is

    wrk : Standard_CSeries_Poly_Systems.Poly_Sys(hom'range);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := Maximum(5,maxdeg/2);
    pv : Standard_Pade_Approximants.Pade_Vector(1..sol.n)
       := Standard_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : Standard_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_Standard_Poles(sol.n,dendeg);
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    maxit : constant natural32 := pars.corsteps;
    extra : constant natural32 := 1;
    fail : boolean;
    t,step : double_float := 0.0;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    err,rco,res,predres : double_float;
    nbrit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in standard_pade_trackers.Track_One_Path 1 ...");
    end if;
    minsize := 1.0; maxsize := 0.0; cntsstp := 0; cntdstp := 0; cntpstp := 0;
    Standard_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0; cntcut := 0; cntfail := 0; nbrsteps := max_steps;
    for k in 1..max_steps loop
      Step_Control(jm,hs,wrk,wrk_sol,maxdeg,nit,pars,pv,poles,t,step,
                   cntsstp,cntdstp,cntpstp,vrblvl-1);
      Predictor_Corrector
        (abh,pv,wrk_sol,predres,t,step,alpha,pars.minsize,tolres,
         maxit,extra,nbrcorrs,err,rco,res,cntcut,cntfail,fail);
      Update_Step_Sizes(minsize,maxsize,step);
      Series_and_Homotopies.Shift(wrk,-step);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
    end loop;
    Standard_Pade_Approximants.Clear(pv);
    Standard_Complex_VecVecs.Clear(poles);
    Standard_CSeries_Poly_Systems.Clear(wrk);
    wrk := Series_and_Homotopies.Shift(hom,-1.0);
    Homotopy_Newton_Steps.Correct
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,extra);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := Standard_Complex_Numbers.Create(t);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntcut,cntfail : out natural32;
                minsize,maxsize : out double_float;
                cntsstp,cntdstp,cntpstp : out natural32;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    wrk : Standard_CSeries_Poly_Systems.Poly_Sys(hom'range);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := Maximum(5,maxdeg/2);
    pv : Standard_Pade_Approximants.Pade_Vector(1..sol.n)
       := Standard_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : Standard_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_Standard_Poles(sol.n,dendeg);
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    maxit : constant natural32 := pars.corsteps;
    extra : constant natural32 := 1;
    fail : boolean;
    t,step : double_float := 0.0;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    err,rco,res,predres : double_float;
    nbrit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in standard_pade_trackers.Track_One_Path 2 ...");
    end if;
    minsize := 1.0; maxsize := 0.0; cntsstp := 0; cntdstp := 0; cntpstp := 0;
    Standard_CSeries_Poly_Systems.Copy(hom,wrk);
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
      Update_Step_Sizes(minsize,maxsize,step);
      Series_and_Homotopies.Shift(wrk,-step);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
    end loop;
    Standard_Pade_Approximants.Clear(pv);
    Standard_Complex_VecVecs.Clear(poles);
    Standard_CSeries_Poly_Systems.Clear(wrk);
    Homotopy_Newton_Steps.Correct
      (file,abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,
       extra,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := Standard_Complex_Numbers.Create(t);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
  end Track_One_Path;

-- VERSIONS ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure Track_One_Path
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out Standard_Complex_Solutions.Solution;
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
    nit : constant integer32 := Maximum(5,maxdeg/2) + 1; -- one extra
    pv : Standard_Pade_Approximants.Pade_Vector(1..sol.n)
       := Standard_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : Standard_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_Standard_Poles(sol.n,dendeg);
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    maxit : constant natural32 := pars.corsteps;
    extra : constant natural32 := 1;
    fail : boolean;
    t,step : double_float := 0.0;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    err,rco,res,predres : double_float;
    nbrit : natural32 := 0;
    wrk_fcf : Standard_Complex_Series_VecVecs.VecVec(fcf'range);

  begin
    if vrblvl > 0
     then put_line("-> in standard_pade_trackers.Track_One_Path 3 ...");
    end if;
    minsize := 1.0; maxsize := 0.0; cntsstp := 0; cntdstp := 0; cntpstp := 0;
    nbrcorrs := 0; cntcut := 0; cntfail := 0; nbrsteps := max_steps;
    wrk_fcf := Standard_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
    for k in 1..max_steps loop
      if mhom = 1 then
        Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
          (wrk_fcf,wrk_sol,t,pars.gamma);
      elsif mhom > 1 then
        Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
          (wrk_fcf,wrk_sol,t,pars.gamma,mhom,idz.all);
      end if;
      Step_Control
        (jm,hs,fhm,wrk_fcf,ejm,mlt,wrk_sol,maxdeg,nit,pars,pv,poles,t,step,
         cntsstp,cntdstp,cntpstp,vrblvl-1);
      Predictor_Corrector
        (abh,pv,wrk_sol,predres,t,step,alpha,pars.minsize,tolres,
         maxit,extra,nbrcorrs,err,rco,res,cntcut,cntfail,fail);
      Update_Step_Sizes(minsize,maxsize,step);
      Standard_CSeries_Vector_Functions.Shift(wrk_fcf,-step);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
    end loop;
    Standard_Pade_Approximants.Clear(pv);
    Standard_Complex_VecVecs.Clear(poles);
    Homotopy_Newton_Steps.Correct
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,extra);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := Standard_Complex_Numbers.Create(t);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_CSeries_Vector_Functions.Deep_Clear(wrk_fcf);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out Standard_Complex_Solutions.Solution;
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
    nit : constant integer32 := Maximum(5,maxdeg/2) + 1; -- one extra
    pv : Standard_Pade_Approximants.Pade_Vector(1..sol.n)
       := Standard_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : Standard_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_Standard_Poles(sol.n,dendeg);
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    maxit : constant natural32 := pars.corsteps;
    extra : constant natural32 := 1;
    fail : boolean;
    t,step : double_float := 0.0;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    err,rco,res,predres : double_float;
    nbrit : natural32 := 0;
    wrk_fcf : Standard_Complex_Series_VecVecs.VecVec(fcf'range);

  begin
    if vrblvl > 0
     then put_line("-> in standard_pade_trackers.Track_One_Path 4 ...");
    end if;
    minsize := 1.0; maxsize := 0.0; cntsstp := 0; cntdstp := 0; cntpstp := 0;
    nbrcorrs := 0; cntcut := 0; cntfail := 0; nbrsteps := max_steps;
    wrk_fcf := Standard_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
    for k in 1..max_steps loop
      if verbose then
        put(file,"Step "); put(file,k,1); put_line(file," : ");
        if vrblvl > 0 then
          if mhom = 1 then
            Homotopy_Coefficient_Scaling.Last_Coefficients
              (file,wrk_fcf(wrk_fcf'last),t,pars.gamma);
          elsif mhom > 1 then
            Homotopy_Coefficient_Scaling.Last_Coefficients
              (file,wrk_fcf,t,pars.gamma,mhom);
          end if;
        end if;
      end if;
      if mhom = 1 then
        if vrblvl > 0 then
          Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
            (file,fhm,wrk_fcf,wrk_sol,t,pars.gamma,true);
        else
          Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
            (file,fhm,wrk_fcf,wrk_sol,t,pars.gamma);
        end if;
      elsif mhom > 1 then
        if vrblvl > 0 then
          Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
            (file,fhm,wrk_fcf,wrk_sol,t,pars.gamma,mhom,idz.all,true);
        else
          Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
            (file,fhm,wrk_fcf,wrk_sol,t,pars.gamma,mhom,idz.all);
        end if;
      end if;
      Step_Control
        (file,verbose,jm,hs,fhm,wrk_fcf,ejm,mlt,wrk_sol,maxdeg,nit,pars,
         pv,poles,t,step,cntsstp,cntdstp,cntpstp,vrblvl-1);
      Predictor_Corrector
        (file,verbose,abh,pv,wrk_sol,predres,t,step,alpha,pars.minsize,tolres,
         maxit,extra,nbrcorrs,err,rco,res,cntcut,cntfail,fail);
      Update_Step_Sizes(minsize,maxsize,step);
      Standard_CSeries_Vector_Functions.Shift(wrk_fcf,-step);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
    end loop;
    Standard_Pade_Approximants.Clear(pv);
    Standard_Complex_VecVecs.Clear(poles);
    Homotopy_Newton_Steps.Correct
      (file,abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,
       extra,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := Standard_Complex_Numbers.Create(t); sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_CSeries_Vector_Functions.Deep_Clear(wrk_fcf);
  end Track_One_Path;

end Standard_Pade_Trackers;
