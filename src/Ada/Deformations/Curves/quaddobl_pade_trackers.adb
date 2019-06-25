with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Homotopy;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_CSeries_Vector_Functions;
with QuadDobl_Pade_Approximants;
with Homotopy_Pade_Approximants;
with Homotopy_Mixed_Residuals;
with Homotopy_Newton_Steps;
with Series_and_Homotopies;
with Series_and_Predictors;

with Standard_Pade_Trackers;

package body QuadDobl_Pade_Trackers is

  function Minimum ( a, b : in double_float ) return double_float is
  begin
    if a < b
     then return a;
     else return b;
    end if;
  end Minimum;

  function Residual_Prediction
              ( sol : QuadDobl_Complex_Vectors.Vector;
                t : double_float ) return double_float is

    res : quad_double;
    qdt : constant quad_double := create(t);
    cmplxt : constant QuadDobl_Complex_Numbers.Complex_Number
           := QuadDobl_Complex_Numbers.Create(qdt);
    val : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_Homotopy.Eval(sol,cmplxt);

  begin
    res := QuadDobl_Complex_Vector_Norms.Max_Norm(val);
    return hihi_part(res);
  end Residual_Prediction;

  function Residual_Prediction
              ( abh : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                sol : QuadDobl_Complex_Vectors.Vector;
                t : double_float ) return double_float is

    res : quad_double;
    qdt : constant quad_double := create(t);
    cmplxt : constant QuadDobl_Complex_Numbers.Complex_Number
           := QuadDobl_Complex_Numbers.Create(qdt);

  begin
    res := Homotopy_mixed_Residuals.Residual(abh,sol,cmplxt);
    return hihi_part(res);
  end Residual_Prediction;

  function Residual_Prediction
              ( file : file_type;
                abh : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                sol : QuadDobl_Complex_Vectors.Vector;
                t : double_float ) return double_float is

    res : quad_double;
    qdt : constant quad_double := create(t);
    cmplxt : constant QuadDobl_Complex_Numbers.Complex_Number
           := QuadDobl_Complex_Numbers.Create(qdt);

  begin
    res := Homotopy_mixed_Residuals.Residual(file,abh,sol,cmplxt);
    return hihi_part(res);
  end Residual_Prediction;

  procedure Track_One_Path
              ( abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                vrblvl : in integer32 := 0 ) is

    wrk : QuadDobl_CSeries_Poly_Systems.Poly_Sys(hom'range);
   -- nbq : constant integer32 := hom'last;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- + 2;
    nit : constant integer32 := integer32(pars.corsteps+2);
    srv : QuadDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : QuadDobl_Complex_Series_Vectors.Vector(hom'range);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range)
       := QuadDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_QuadDobl_Poles(sol.n,dendeg);
   -- tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 50;
    extra : constant natural32 := maxit/10;
    fail : boolean;
    t,step,dstep : double_float := 0.0;
    qd_t,qd_step : quad_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : QuadDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : double_float;
    frp : quad_double;
    cfp : QuadDobl_Complex_Numbers.Complex_Number;
    predres : double_float;
    nbrit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in series_and_trackers.Track_One_Path 3 ...");
    end if;
    minsize := 1.0; maxsize := 0.0;
    QuadDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction(maxdeg,nit,wrk,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
     -- step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
     -- step := pars.sbeta*step;
     -- step := Series_and_Predictors.Cap_Step_Size
     --           (step,hihi_part(frp),pars.pbeta); -- ignore series step
      QuadDobl_Complex_Series_Vectors.Clear(eva);
      qd_t := Create(t);
     -- dstep := Series_and_Predictors.Step_Distance
     --            (maxdeg,dbeta,qd_t,jm,hs,wrk_sol,srv,pv);
      dstep := pars.maxsize; -- ignore Hessian step
      step := Series_and_Predictors.Cap_Step_Size
                (dstep,hihi_part(frp),pars.pbeta);
     -- step := Minimum(step,dstep);
      Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
      loop
        loop
          qd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
         -- predres := Residual_Prediction(wrk_sol,t);
          predres := Residual_Prediction(abh,wrk_sol,t);
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
         -- exit when (step < pars.minsize);
          exit when (step <= alpha);
        end loop;
        Standard_Pade_Trackers.Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,extra);
         -- (nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,extra);
        nbrcorrs := nbrcorrs + nbrit;
        exit when (not fail);
        step := step/2.0; cntfail := cntfail + 1;
        exit when (step < pars.minsize);
      end loop;
      QuadDobl_Complex_Series_Vectors.Clear(srv);
      qd_step := create(step);
      Series_and_Homotopies.Shift(wrk,-qd_step);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
    end loop;
    QuadDobl_Pade_Approximants.Clear(pv);
    QuadDobl_Complex_VecVecs.Clear(poles);
    QuadDobl_CSeries_Poly_Systems.Clear(wrk);
    qd_t := create(-1.0);
    wrk := Series_and_Homotopies.Shift(hom,qd_t);
    Homotopy_Newton_Steps.Correct
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,extra);
     -- (nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := QuadDobl_Complex_Numbers.Create(Quad_Double_Numbers.Create(t));
    sol.v := wrk_sol;
    sol.err := Quad_Double_Numbers.create(err);
    sol.rco := Quad_Double_Numbers.create(rco);
    sol.res := Quad_Double_Numbers.create(res);
    QuadDobl_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    wrk : QuadDobl_CSeries_Poly_Systems.Poly_Sys(hom'range);
   -- nbq : constant integer32 := hom'last;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- + 2;
    nit : constant integer32 := integer32(pars.corsteps+2);
    srv : QuadDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : QuadDobl_Complex_Series_Vectors.Vector(hom'range);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range)
       := QuadDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_QuadDobl_Poles(sol.n,dendeg);
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 50;
    extra : constant natural32 := maxit/10;
    fail : boolean;
    t,step,dstep : double_float := 0.0;
    qd_t,qd_step : quad_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : QuadDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : double_float;
    frp : quad_double;
    cfp : QuadDobl_Complex_Numbers.Complex_Number;
    predres : double_float;
    nbrit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in series_and_trackers.Track_One_Path 6 ...");
    end if;
    minsize := 1.0; maxsize := 0.0;
    QuadDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      if verbose then
        put(file,"Step "); put(file,k,1); put(file," : ");
      end if;
      Series_and_Predictors.Newton_Prediction
        (file,maxdeg,nit,wrk,wrk_sol,srv,eva,false); -- verbose);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
      if verbose then
        put(file,"Smallest pole radius : ");
        put(file,frp,3); new_line(file);
        put(file,"Closest pole : "); put(file,cfp); new_line(file);
      end if;
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,alpha,verbose);
     -- step := pars.sbeta*step;
      step := pars.maxsize; -- ignore series step
      step := Series_and_Predictors.Cap_Step_Size
                (step,hihi_part(frp),pars.pbeta);
      QuadDobl_Complex_Series_Vectors.Clear(eva);
      qd_t := Create(t);
     -- dstep := Series_and_Predictors.Step_Distance
     --            (maxdeg,dbeta,qd_t,jm,hs,wrk_sol,srv,pv);
      dstep := pars.maxsize; -- ignore Hessian step
      step := Minimum(step,dstep);
      Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);
      end if;
      loop
        loop
          qd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
         -- predres := Residual_Prediction(wrk_sol,t);
          predres := Residual_Prediction(abh,wrk_sol,t);
          if verbose
           then put(file,"  residual : "); put(file,predres,3); new_line(file);
          end if;
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
          if verbose then
            put(file,"Step size : "); put(file,step,3);
            put(file," t = "); put(file,t,3);
          end if;
         -- exit when (step < pars.minsize);
          exit when (step < alpha);
        end loop;
        Standard_Pade_Trackers.Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (file,abh,t,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,
           extra,verbose);
         -- (file,nbq,t,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,
         --  extra,verbose);
        nbrcorrs := nbrcorrs + nbrit;
        exit when (not fail);
        step := step/2.0; cntfail := cntfail + 1;
        exit when (step < pars.minsize);
      end loop;
      QuadDobl_Complex_Series_Vectors.Clear(srv);
     -- QuadDobl_CSeries_Poly_Systems.Clear(wrk);
     -- qd_t := create(-t);
     -- wrk := Series_and_Homotopies.Shift(hom,qd_t);
      qd_step := create(step);
      Series_and_Homotopies.Shift(wrk,-qd_step);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
    end loop;
    QuadDobl_Pade_Approximants.Clear(pv);
    QuadDobl_Complex_VecVecs.Clear(poles);
    QuadDobl_CSeries_Poly_Systems.Clear(wrk);
    qd_t := create(-1.0);
    wrk := Series_and_Homotopies.Shift(hom,qd_t);
    Homotopy_Newton_Steps.Correct
      (file,abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,extra,verbose);
     -- (file,nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,
     --  extra,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := QuadDobl_Complex_Numbers.Create(Quad_Double_Numbers.Create(t));
    sol.v := wrk_sol;
    sol.err := Quad_Double_Numbers.Create(err);
    sol.rco := Quad_Double_Numbers.Create(rco);
    sol.res := Quad_Double_Numbers.Create(res);
    QuadDobl_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

-- VERSIONS ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure Track_One_Path
              ( abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out QuadDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                vrblvl : in integer32 := 0 ) is

   -- nbq : constant integer32 := fhm'last;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- + 2;
    nit : constant integer32 := integer32(pars.corsteps+2);
    srv : QuadDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : QuadDobl_Complex_Series_Vectors.Vector(fhm'range);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range)
       := QuadDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_QuadDobl_Poles(sol.n,dendeg);
   -- tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 50;
    extra : constant natural32 := maxit/10;
    fail : boolean;
    t,step,dstep : double_float := 0.0;
    qd_t,qd_step : quad_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : QuadDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res,predres : double_float;
    frp : quad_double;
    cfp : QuadDobl_Complex_Numbers.Complex_Number;
    nbrit : natural32 := 0;
    wrk_fcf : QuadDobl_Complex_Series_VecVecs.VecVec(fcf'range);

  begin
    if vrblvl > 0
     then put_line("-> in series_and_trackers.Track_One_Path 9 ...");
    end if;
    minsize := 1.0; maxsize := 0.0;
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    wrk_fcf := QuadDobl_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction
        (maxdeg,nit,fhm,wrk_fcf,ejm,mlt,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
     -- step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
     -- step := pars.sbeta*step;
      QuadDobl_Complex_Series_Vectors.Clear(eva);
     -- step := Series_and_Predictors.Cap_Step_Size
     --           (step,hihi_part(frp),pars.pbeta);
      qd_t := Create(t);
     -- dstep := Series_and_Predictors.Step_Distance
     --            (maxdeg,dbeta,qd_t,jm,hs,wrk_sol,srv,pv);
     -- step := Minimum(step,dstep); -- ignore series step
      dstep := pars.maxsize; -- ignore Hessian step
      step := Series_and_Predictors.Cap_Step_Size
                (dstep,hihi_part(frp),pars.pbeta);
      Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
      loop
        loop
          qd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
         -- predres := Residual_Prediction(wrk_sol,t);
          predres := Residual_Prediction(abh,wrk_sol,t);
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
         -- exit when (step < pars.minsize);
          exit when (step <= alpha);
        end loop;
        Standard_Pade_Trackers.Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,extra);
         -- (nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,extra);
        nbrcorrs := nbrcorrs + nbrit;
        exit when (not fail);
        step := step/2.0; cntfail := cntfail + 1;
        exit when (step < pars.minsize);
      end loop;
      QuadDobl_Complex_Series_Vectors.Clear(srv);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
      qd_step := create(step);
      QuadDobl_CSeries_Vector_Functions.Shift(wrk_fcf,-qd_step);
    end loop;
    QuadDobl_Pade_Approximants.Clear(pv);
    QuadDobl_Complex_VecVecs.Clear(poles);
    Homotopy_Newton_Steps.Correct
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,extra);
     -- (nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,extra);
    nbrcorrs := nbrcorrs + nbrit;
    qd_t := Quad_Double_Numbers.Create(t);
    sol.t := QuadDobl_Complex_Numbers.Create(qd_t);
    sol.v := wrk_sol;
    sol.err := Quad_Double_Numbers.Create(err);
    sol.rco := Quad_Double_Numbers.Create(rco);
    sol.res := Quad_Double_Numbers.Create(res);
    QuadDobl_CSeries_Vector_Functions.Deep_Clear(wrk_fcf);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out QuadDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

   -- nbq : constant integer32 := fhm'last;
    nit : constant integer32 := integer32(pars.corsteps+2);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- 2;
    srv : QuadDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : QuadDobl_Complex_Series_Vectors.Vector(fhm'range);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range)
       := QuadDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_QuadDobl_Poles(sol.n,dendeg);
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 50;
    extra : constant natural32 := maxit/10;
    fail : boolean;
    t,step,dstep,pstep : double_float := 0.0;
    qd_t,qd_step : quad_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : QuadDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res,predres : double_float;
    frp : quad_double;
    cfp : QuadDobl_Complex_Numbers.Complex_Number;
    nbrit : natural32 := 0;
    wrk_fcf : QuadDobl_Complex_Series_VecVecs.VecVec(fcf'range);

  begin
    if vrblvl > 0
     then put_line("-> in series_and_trackers.Track_One_Path 12 ...");
    end if;
    minsize := 1.0; maxsize := 0.0;
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    wrk_fcf := QuadDobl_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
    for k in 1..max_steps loop
      if verbose then
        put(file,"Step "); put(file,k,1); put(file," : ");
      end if;
      Series_and_Predictors.Newton_Prediction -- verbose flag set to false
        (file,maxdeg,nit,fhm,wrk_fcf,ejm,mlt,wrk_sol,srv,eva,false);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
      if verbose then
        put(file,"Smallest pole radius :");
        put(file,frp,3); new_line(file);
        put(file,"Closest pole :"); put(file,cfp); new_line(file);
      end if;
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,alpha,verbose);
      step := pars.sbeta*step;
      if verbose then
        put(file,"series step : "); put(file,step,2);
      end if;
     -- step := Series_and_Predictors.Cap_Step_Size
     --           (step,hihi_part(frp),pars.pbeta);
      QuadDobl_Complex_Series_Vectors.Clear(eva);
      qd_t := Create(t);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,qd_t,jm,hs,wrk_sol,srv,pv);
      if verbose then
        put(file,"  Hessian step : "); put(file,dstep,2);
        pstep := hihi_part(frp)*pars.pbeta;
        put(file,"  pole step : "); put(file,pstep,2); new_line(file);
      end if;
      dstep := pars.maxsize; -- ignore Hessian step
      step := Minimum(dstep,pstep);
      Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);
      end if;
      loop
        loop
          qd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
         -- predres := Residual_Prediction(wrk_sol,t);
          predres := Residual_Prediction(abh,wrk_sol,t);
          if verbose
           then put(file,"  residual : "); put(file,predres,3); new_line(file);
          end if;
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
          if verbose then
            put(file,"Step size : "); put(file,step,3);
            put(file," t = "); put(file,t,3);
          end if;
         -- exit when (step < pars.minsize);
          exit when (step <= alpha);
        end loop;
        Standard_Pade_Trackers.Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (file,abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,
           extra,verbose);
         -- (file,nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,
         --  extra,verbose);
        nbrcorrs := nbrcorrs + nbrit;
        exit when (not fail);
        step := step/2.0; cntfail := cntfail + 1;
        exit when (step < pars.minsize);
      end loop;
      QuadDobl_Complex_Series_Vectors.Clear(srv);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
      qd_step := create(step);
      QuadDobl_CSeries_Vector_Functions.Shift(wrk_fcf,-qd_step);
    end loop;
    QuadDobl_Pade_Approximants.Clear(pv);
    QuadDobl_Complex_VecVecs.Clear(poles);
    Homotopy_Newton_Steps.Correct
      (file,abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,
       extra,verbose);
     -- (file,nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,
     --  extra,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    qd_t := Quad_Double_Numbers.Create(t);
    sol.t := QuadDobl_Complex_Numbers.Create(qd_t);
    sol.v := wrk_sol;
    sol.err := Quad_Double_Numbers.Create(err);
    sol.rco := Quad_Double_Numbers.Create(rco);
    sol.res := Quad_Double_Numbers.Create(res);
    QuadDobl_CSeries_Vector_Functions.Deep_Clear(wrk_fcf);
  end Track_One_Path;

end QuadDobl_Pade_Trackers;
