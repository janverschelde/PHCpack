with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Homotopy;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;
with DoblDobl_Complex_Series_Vectors;
with Standard_CSeries_Jaco_Matrices;
with DoblDobl_CSeries_Vector_Functions;
with DoblDobl_Pade_Approximants;
with Homotopy_Pade_Approximants;
with Homotopy_Mixed_Residuals;
with Homotopy_Newton_Steps;
with Series_and_Homotopies;
with Series_and_Predictors;

with Standard_Pade_Trackers;

package body DoblDobl_Pade_Trackers is

  function Minimum ( a, b : in double_float ) return double_float is
  begin
    if a < b
     then return a;
     else return b;
    end if;
  end Minimum;

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

  procedure Track_One_Path
              ( abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                vrblvl : in integer32 := 0 ) is

    wrk : DoblDobl_CSeries_Poly_Systems.Poly_Sys(hom'range);
   -- nbq : constant integer32 := hom'last;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- + 2;
    nit : constant integer32 := integer32(pars.corsteps+2);
    srv : DoblDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : DoblDobl_Complex_Series_Vectors.Vector(hom'range);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range)
       := DoblDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_DoblDobl_Poles(sol.n,dendeg);
   -- tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 50;
    extra : constant natural32 := maxit/10;
    fail : boolean;
    t,step,dstep : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : double_float;
    frp : double_double;
    cfp : DoblDobl_Complex_Numbers.Complex_Number;
    predres : double_float;
    nbrit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Track_One_Path 1 ...");
    end if;
    minsize := 1.0; maxsize := 0.0;
    DoblDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction(maxdeg,nit,wrk,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
     -- step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
     -- step := pars.sbeta*step;
     -- step := Series_and_Predictors.Cap_Step_Size
     --            (step,hi_part(frp),pars.pbeta); -- ignore series step
      DoblDobl_Complex_Series_Vectors.Clear(eva);
      dd_t := Create(t);
     -- dstep := Series_and_Predictors.Step_Distance
     --            (maxdeg,dbeta,dd_t,jm,hs,wrk_sol,srv,pv);
      dstep := pars.maxsize; -- ignore Hessian step size
      step := Series_and_Predictors.Cap_Step_Size
                 (dstep,hi_part(frp),pars.pbeta);
     -- step := Minimum(step,dstep);
      Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
      loop
        loop
          dd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
         -- predres := Residual_Prediction(wrk_sol,t);
          predres := Residual_Prediction(abh,wrk_sol,t);
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
         -- exit when (step < pars.minsize);
          exit when (step <= alpha);
        end loop;
        Standard_Pade_Trackers.Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > pars.alpha));
        Homotopy_Newton_Steps.Correct
          (abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,extra);
         -- (nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,extra);
        nbrcorrs := nbrcorrs + nbrit;
        exit when (not fail);
        step := step/2.0; cntfail := cntfail + 1;
        exit when (step < pars.minsize);
      end loop;
      DoblDobl_Complex_Series_Vectors.Clear(srv);
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
     -- (nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,extra);
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
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    wrk : DoblDobl_CSeries_Poly_Systems.Poly_Sys(hom'range);
   -- nbq : constant integer32 := hom'last;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- 2;
    nit : constant integer32 := integer32(pars.corsteps+2);
    srv : DoblDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : DoblDobl_Complex_Series_Vectors.Vector(hom'range);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range)
       := DoblDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_DoblDobl_Poles(sol.n,dendeg);
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 50;
    extra : constant natural32 := maxit/10;
    fail : boolean;
    t,step,dstep : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : double_float;
    frp : double_double;
    cfp : DoblDobl_Complex_Numbers.Complex_Number;
    predres : double_float;
    nbrit : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Track_One_Path 2 ...");
    end if;
    minsize := 1.0; maxsize := 0.0;
    DoblDobl_CSeries_Poly_Systems.Copy(hom,wrk);
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
      step := pars.sbeta*step;
      step := Series_and_Predictors.Cap_Step_Size
                (step,hi_part(frp),pars.pbeta);
      DoblDobl_Complex_Series_Vectors.Clear(eva);
      dd_t := Create(t);
     -- dstep := Series_and_Predictors.Step_Distance
     --            (maxdeg,dbeta,dd_t,jm,hs,wrk_sol,srv,pv);
      dstep := pars.maxsize; -- ignore Hessian step
      step := Minimum(step,dstep);
      Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);
      end if;
      loop
        loop
          dd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
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
      DoblDobl_Complex_Series_Vectors.Clear(srv);
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
     -- (file,nbq,1.0,tolres,pars.corsteps,nbrit,
     --  wrk_sol,err,rco,res,fail,extra,verbose);
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
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                vrblvl : in integer32 := 0 ) is

   -- nbq : constant integer32 := fhm'last;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- 2;
    nit : constant integer32 := integer32(pars.corsteps+2);
    srv : DoblDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : DoblDobl_Complex_Series_Vectors.Vector(fhm'range);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range)
       := DoblDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_DoblDobl_Poles(sol.n,dendeg);
   -- tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 50;
    extra : constant natural32 := maxit/10;
    fail : boolean;
    t,step,dstep : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res,predres : double_float;
    frp : double_double;
    cfp : DoblDobl_Complex_Numbers.Complex_Number;
    nbrit : natural32 := 0;
    wrk_fcf : DoblDobl_Complex_Series_VecVecs.VecVec(fcf'range);

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Track_One_Path 3 ...");
    end if;
    minsize := 1.0; maxsize := 0.0;
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    wrk_fcf := DoblDobl_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction
        (maxdeg,nit,fhm,wrk_fcf,ejm,mlt,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
     -- step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
     -- step := pars.sbeta*step;
      DoblDobl_Complex_Series_Vectors.Clear(eva);
     -- step := Series_and_Predictors.Cap_Step_Size
     --           (step,hi_part(frp),pars.pbeta);
      dd_t := Create(t);
     -- dstep := Series_and_Predictors.Step_Distance
     --            (maxdeg,dbeta,dd_t,jm,hs,wrk_sol,srv,pv);
     -- step := Minimum(step,dstep); -- ignore series step
      dstep := pars.maxsize; -- ignore Hessian step
      step := Series_and_Predictors.Cap_Step_Size
                (dstep,hi_part(frp),pars.pbeta);
      Standard_Pade_Trackers.Set_Step(t,step,pars.maxsize,onetarget);
      loop
        loop
          dd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
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
      DoblDobl_Complex_Series_Vectors.Clear(srv);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
      dd_step := create(step);
      DoblDobl_CSeries_Vector_Functions.Shift(wrk_fcf,-dd_step);
    end loop;
    DoblDobl_Pade_Approximants.Clear(pv);
    DoblDobl_Complex_VecVecs.Clear(poles);
    Homotopy_Newton_Steps.Correct
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,extra);
     -- (nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,extra);
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
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

   -- nbq : constant integer32 := fhm'last;
    nit : constant integer32 := integer32(pars.corsteps+2);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- 2;
    srv : DoblDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : DoblDobl_Complex_Series_Vectors.Vector(fhm'range);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range)
       := DoblDobl_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_DoblDobl_Poles(sol.n,dendeg);
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 50;
    extra : constant natural32 := maxit/10;
    fail : boolean;
    t,step,dstep,pstep : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res,predres : double_float;
    frp : double_double;
    cfp : DoblDobl_Complex_Numbers.Complex_Number;
    nbrit : natural32 := 0;
    wrk_fcf : DoblDobl_Complex_Series_VecVecs.VecVec(fcf'range);

  begin
    if vrblvl > 0
     then put_line("-> in dobldobl_pade_trackers.Track_One_Path 4 ...");
    end if;
    minsize := 1.0; maxsize := 0.0;
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    wrk_fcf := DoblDobl_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
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
        put(file,"Closest pole : "); put(file,cfp); new_line(file);
      end if;
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,alpha,verbose);
      step := pars.sbeta*step;
      if verbose then
        put(file,"series step : "); put(file,step,2);
      end if;
     -- step := Series_and_Predictors.Cap_Step_Size
     --           (step,hi_part(frp),pars.pbeta);
      DoblDobl_Complex_Series_Vectors.Clear(eva);
      dd_t := Create(t);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,dd_t,jm,hs,wrk_sol,srv,pv);
      if verbose then
        put(file,"  Hessian step : "); put(file,dstep,2);
        pstep := hi_part(frp)*pars.pbeta;
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
          dd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
         -- predres := Residual_Prediction(wrk_sol,t);
          if not verbose then
            predres := Residual_Prediction(abh,wrk_sol,t);
          else
            predres := Residual_Prediction(file,abh,wrk_sol,t);
            put(file,"  residual : "); put(file,predres,3); new_line(file);
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
      DoblDobl_Complex_Series_Vectors.Clear(srv);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
      dd_step := create(step);
      DoblDobl_CSeries_Vector_Functions.Shift(wrk_fcf,-dd_step);
    end loop;
    DoblDobl_Pade_Approximants.Clear(pv);
    DoblDobl_Complex_VecVecs.Clear(poles);
    Homotopy_Newton_Steps.Correct
      (file,abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,
       extra,verbose);
     -- (file,nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail,
     --  extra,verbose);
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
