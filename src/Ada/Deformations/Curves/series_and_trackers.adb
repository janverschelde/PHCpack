with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;  
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Complex_VecVecs;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Vector_Norms;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;
with Standard_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with Standard_CSeries_Vector_Functions;
with DoblDobl_CSeries_Vector_Functions;
with QuadDobl_CSeries_Vector_Functions;
with Standard_Pade_Approximants;
with DoblDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants;
with Homotopy_Pade_Approximants;
with Homotopy_Mixed_Residuals;
with Homotopy_Newton_Steps;
with Series_and_Homotopies;
with Series_and_Predictors;

package body Series_and_Trackers is

  function Minimum ( a, b : in double_float ) return double_float is
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

  procedure Track_One_Path
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float ) is

    wrk : Standard_CSeries_Poly_Systems.Poly_Sys(hom'range);
   -- nbq : constant integer32 := hom'last;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- + 2;
    nit : constant integer32 := integer32(pars.corsteps+2);
    srv : Standard_Complex_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Complex_Series_Vectors.Vector(hom'range);
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range)
       := Standard_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : Standard_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_Standard_Poles(sol.n,dendeg);
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 500;
    fail : boolean;
    t,step,dstep : double_float := 0.0;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res,frp,predres : double_float;
    cfp : Standard_Complex_Numbers.Complex_Number;
    nbrit : natural32 := 0;

  begin
    minsize := 1.0; maxsize := 0.0;
    Standard_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction(maxdeg,nit,wrk,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
      step := pars.sbeta*step;
      Standard_Complex_Series_Vectors.Clear(eva);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,t,jm,hs,wrk_sol,srv,pv);
      step := Series_and_Predictors.Cap_Step_Size(step,frp,pars.pbeta);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      loop
        loop
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,step);
          predres := Residual_Prediction(wrk_sol,t);
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
         -- exit when (step < pars.minsize);
          exit when (step <= alpha);
        end loop;
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
         -- (nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
        nbrcorrs := nbrcorrs + nbrit;
        exit when (not fail);
        step := step/2.0; cntfail := cntfail + 1;
        exit when (step < pars.minsize);
      end loop;
      Standard_Complex_Series_Vectors.Clear(srv);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
      Series_and_Homotopies.Shift(wrk,-step);
    end loop;
    Standard_Pade_Approximants.Clear(pv);
    Standard_Complex_VecVecs.Clear(poles);
    Standard_CSeries_Poly_Systems.Clear(wrk);
    wrk := Series_and_Homotopies.Shift(hom,-1.0);
    Homotopy_Newton_Steps.Correct
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
     -- (nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := Standard_Complex_Numbers.Create(t);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float ) is

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
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 500;
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
    minsize := 1.0; maxsize := 0.0;
    DoblDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction(maxdeg,nit,wrk,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
      step := pars.sbeta*step;
      step := Series_and_Predictors.Cap_Step_Size
                 (step,hi_part(frp),pars.pbeta);
      DoblDobl_Complex_Series_Vectors.Clear(eva);
      dd_t := Create(t);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,dd_t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      loop
        loop
          dd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
          predres := Residual_Prediction(wrk_sol,t);
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
         -- exit when (step < pars.minsize);
          exit when (step <= alpha);
        end loop;
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > pars.alpha));
        Homotopy_Newton_Steps.Correct
          (abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
         -- (nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
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
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
     -- (nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := DoblDobl_Complex_Numbers.Create(Double_Double_Numbers.Create(t));
    sol.v := wrk_sol;
    sol.err := Double_Double_Numbers.create(err);
    sol.rco := Double_Double_Numbers.create(rco);
    sol.res := Double_Double_Numbers.create(res);
    DoblDobl_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float ) is

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
    maxit : constant natural32 := 500;
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
    minsize := 1.0; maxsize := 0.0;
    QuadDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction(maxdeg,nit,wrk,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
      step := pars.sbeta*step;
      step := Series_and_Predictors.Cap_Step_Size
                (step,hihi_part(frp),pars.pbeta);
      QuadDobl_Complex_Series_Vectors.Clear(eva);
      qd_t := Create(t);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,qd_t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      loop
        loop
          qd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
          predres := Residual_Prediction(wrk_sol,t);
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
         -- exit when (step < pars.minsize);
          exit when (step <= alpha);
        end loop;
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
         -- (nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
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
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
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
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false ) is

    wrk : Standard_CSeries_Poly_Systems.Poly_Sys(hom'range);
   -- nbq : constant integer32 := hom'last;
    nit : constant integer32 := integer32(pars.corsteps+2);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- + 2;
    srv : Standard_Complex_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Complex_Series_Vectors.Vector(hom'range);
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range)
       := Standard_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : Standard_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_Standard_Poles(sol.n,dendeg);
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 500;
    fail : boolean;
    t,step,dstep : double_float := 0.0;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res,frp,predres : double_float;
    cfp : Standard_Complex_Numbers.Complex_Number;
    nbrit : natural32 := 0;

  begin
    minsize := 1.0; maxsize := 0.0;
    Standard_CSeries_Poly_Systems.Copy(hom,wrk);
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
        put(file,"Smallest pole radius :");
        put(file,frp,3); new_line(file);
        if Standard_Complex_Numbers.REAL_PART(cfp) >= 0.0
         then put(file,"Closest pole :"); put(file,cfp); new_line(file);
        end if;
      end if;
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,alpha,verbose);
      step := pars.sbeta*step;
      step := Series_and_Predictors.Cap_Step_Size(step,frp,pars.pbeta);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Standard_Complex_Series_Vectors.Clear(eva);
      Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);
      end if;
      loop
        loop
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,step);
          predres := Residual_Prediction(wrk_sol,t);
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
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (file,abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
         -- (file,nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
        nbrcorrs := nbrcorrs + nbrit;
        exit when (not fail);
        step := step/2.0; cntfail := cntfail + 1;
        exit when (step < pars.minsize);
      end loop;
      Standard_Complex_Series_Vectors.Clear(srv);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
      Series_and_Homotopies.Shift(wrk,-step);
    end loop;
    Standard_Pade_Approximants.Clear(pv);
    Standard_Complex_VecVecs.Clear(poles);
    Standard_CSeries_Poly_Systems.Clear(wrk);
    Homotopy_Newton_Steps.Correct
      (file,abh,1.0,tolres,pars.corsteps,nbrit,
       wrk_sol,err,rco,res,fail,verbose);
     -- (file,nbq,1.0,tolres,pars.corsteps,nbrit,
     --  wrk_sol,err,rco,res,fail,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := Standard_Complex_Numbers.Create(t);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
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
                verbose : in boolean := false ) is

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
    maxit : constant natural32 := 500;
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
        if DoblDobl_Complex_Numbers.REAL_PART(cfp) >= 0.0
         then put(file,"Closest pole : "); put(file,cfp); new_line(file);
        end if;
      end if;
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,alpha,verbose);
      step := pars.sbeta*step;
      step := Series_and_Predictors.Cap_Step_Size
                (step,hi_part(frp),pars.pbeta);
      DoblDobl_Complex_Series_Vectors.Clear(eva);
      dd_t := Create(t);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,dd_t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);
      end if;
      loop
        loop
          dd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
          predres := Residual_Prediction(wrk_sol,t);
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
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (file,abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
         -- (file,nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
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
      (file,abh,1.0,tolres,pars.corsteps,nbrit,
       wrk_sol,err,rco,res,fail,verbose);
     -- (file,nbq,1.0,tolres,pars.corsteps,nbrit,
     --  wrk_sol,err,rco,res,fail,verbose);
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
                abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false ) is

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
    maxit : constant natural32 := 500;
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
        if QuadDobl_Complex_Numbers.REAL_PART(cfp) >= 0.0
         then put(file,"Closest pole : "); put(file,cfp); new_line(file);
        end if;
      end if;
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,alpha,verbose);
      step := pars.sbeta*step;
      step := Series_and_Predictors.Cap_Step_Size
                (step,hihi_part(frp),pars.pbeta);
      QuadDobl_Complex_Series_Vectors.Clear(eva);
      qd_t := Create(t);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,qd_t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);
      end if;
      loop
        loop
          qd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
          predres := Residual_Prediction(wrk_sol,t);
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
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (file,abh,t,tolres,pars.corsteps,nbrit,
           wrk_sol,err,rco,res,fail,verbose);
         -- (file,nbq,t,tolres,pars.corsteps,nbrit,
         --  wrk_sol,err,rco,res,fail,verbose);
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
      (file,abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
     -- (file,nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
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
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float ) is

   -- nbq : constant integer32 := fhm'last;
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- + 2;
    nit : constant integer32 := integer32(pars.corsteps+2);
    srv : Standard_Complex_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Complex_Series_Vectors.Vector(fhm'range);
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range)
       := Standard_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : Standard_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_Standard_Poles(sol.n,dendeg);
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 500;
    fail : boolean;
    t,step,dstep : double_float := 0.0;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res,frp,predres : double_float;
    cfp : Standard_Complex_Numbers.Complex_Number;
    nbrit : natural32 := 0;
    wrk_fcf : Standard_Complex_Series_VecVecs.VecVec(fcf'range);

  begin
    minsize := 1.0; maxsize := 0.0;
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    wrk_fcf := Standard_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction
        (maxdeg,nit,fhm,wrk_fcf,ejm,mlt,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
      step := pars.sbeta*step;
      Standard_Complex_Series_Vectors.Clear(eva);
      step := Series_and_Predictors.Cap_Step_Size(step,frp,pars.pbeta);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      loop
        loop
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,step);
          predres := Residual_Prediction(wrk_sol,t);
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
         -- exit when (step < pars.minsize);
          exit when (step <= alpha);
        end loop;
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
         -- (nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
        nbrcorrs := nbrcorrs + nbrit;
        exit when (not fail);
        step := step/2.0; cntfail := cntfail + 1;
        exit when (step < pars.minsize);
      end loop;
      Standard_Complex_Series_Vectors.Clear(srv);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
      Standard_CSeries_Vector_Functions.Shift(wrk_fcf,-step);
    end loop;
    Standard_Pade_Approximants.Clear(pv);
    Standard_Complex_VecVecs.Clear(poles);
    Homotopy_Newton_Steps.Correct
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
     -- (nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := Standard_Complex_Numbers.Create(t);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_CSeries_Vector_Functions.Deep_Clear(wrk_fcf);
  end Track_One_Path;

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
                minsize,maxsize : out double_float ) is

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
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 500;
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
    minsize := 1.0; maxsize := 0.0;
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    wrk_fcf := DoblDobl_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction
        (maxdeg,nit,fhm,wrk_fcf,ejm,mlt,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
      step := pars.sbeta*step;
      DoblDobl_Complex_Series_Vectors.Clear(eva);
      step := Series_and_Predictors.Cap_Step_Size
                (step,hi_part(frp),pars.pbeta);
      dd_t := Create(t);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,dd_t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      loop
        loop
          dd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
          predres := Residual_Prediction(wrk_sol,t);
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
         -- exit when (step < pars.minsize);
          exit when (step <= alpha);
        end loop;
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
         -- (nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
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
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
     -- (nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
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
                minsize,maxsize : out double_float ) is

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
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 500;
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
    minsize := 1.0; maxsize := 0.0;
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    wrk_fcf := QuadDobl_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction
        (maxdeg,nit,fhm,wrk_fcf,ejm,mlt,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(srv,pv,poles,frp,cfp);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
      step := pars.sbeta*step;
      QuadDobl_Complex_Series_Vectors.Clear(eva);
      step := Series_and_Predictors.Cap_Step_Size
                (step,hihi_part(frp),pars.pbeta);
      qd_t := Create(t);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,qd_t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      loop
        loop
          qd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
          predres := Residual_Prediction(wrk_sol,t);
          exit when (predres <= alpha);
          t := t - step; step := step/2.0; t := t + step;
         -- exit when (step < pars.minsize);
          exit when (step <= alpha);
        end loop;
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
         -- (nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail);
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
      (abh,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
     -- (nbq,1.0,tolres,pars.corsteps,nbrit,wrk_sol,err,rco,res,fail);
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
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntfail : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false ) is

   -- nbq : constant integer32 := fhm'last;
    nit : constant integer32 := integer32(pars.corsteps+2);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- + 2;
    srv : Standard_Complex_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Complex_Series_Vectors.Vector(fhm'range);
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range)
       := Standard_Pade_Approximants.Allocate(sol.n,numdeg,dendeg);
    poles : Standard_Complex_VecVecs.VecVec(pv'range)
          := Homotopy_Pade_Approximants.Allocate_Standard_Poles(sol.n,dendeg);
    tolcff : constant double_float := pars.epsilon;
    alpha : constant double_float := pars.alpha;
    tolres : constant double_float := pars.tolres;
    dbeta : constant double_float := 0.005;
    maxit : constant natural32 := 500;
    fail : boolean;
    t,step,dstep : double_float := 0.0;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res,frp,predres : double_float;
    cfp : Standard_Complex_Numbers.Complex_Number;
    nbrit : natural32 := 0;
    wrk_fcf : Standard_Complex_Series_VecVecs.VecVec(fcf'range);

  begin
    minsize := 1.0; maxsize := 0.0;
    nbrcorrs := 0; cntfail := 0;
    nbrsteps := max_steps;
    wrk_fcf := Standard_CSeries_Vector_Functions.Make_Deep_Copy(fcf);
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
        if Standard_Complex_Numbers.REAL_PART(cfp) >= 0.0
         then put(file,"Closest pole :"); put(file,cfp); new_line(file);
        end if;
      end if;
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,alpha,verbose);
      step := pars.sbeta*step;
      step := Series_and_Predictors.Cap_Step_Size(step,frp,pars.pbeta);
      Standard_Complex_Series_Vectors.Clear(eva);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);
      end if;
      loop
        loop
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,step);
          predres := Residual_Prediction(wrk_sol,t);
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
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (file,abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
         -- (file,nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
        nbrcorrs := nbrcorrs + nbrit;
        exit when (not fail);
        step := step/2.0; cntfail := cntfail + 1;
        exit when (step < pars.minsize);
      end loop;
      Standard_Complex_Series_Vectors.Clear(srv);
      if t = 1.0 then        -- converged and reached the end
        nbrsteps := k; exit;
      elsif (fail and (step < pars.minsize)) then -- diverged
        nbrsteps := k; exit;
      end if;
      Standard_CSeries_Vector_Functions.Shift(wrk_fcf,-step);
    end loop;
    Standard_Pade_Approximants.Clear(pv);
    Standard_Complex_VecVecs.Clear(poles);
    Homotopy_Newton_Steps.Correct
      (file,abh,1.0,tolres,pars.corsteps,nbrit,
       wrk_sol,err,rco,res,fail,verbose);
     -- (file,nbq,1.0,tolres,pars.corsteps,nbrit,
     --  wrk_sol,err,rco,res,fail,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := Standard_Complex_Numbers.Create(t);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_CSeries_Vector_Functions.Deep_Clear(wrk_fcf);
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
                verbose : in boolean := false ) is

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
    maxit : constant natural32 := 500;
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
        if DoblDobl_Complex_Numbers.REAL_PART(cfp) >= 0.0
         then put(file,"Closest pole :"); put(file,cfp); new_line(file);
        end if;
      end if;
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,alpha,verbose);
      step := pars.sbeta*step;
      step := Series_and_Predictors.Cap_Step_Size
                (step,hi_part(frp),pars.pbeta);
      DoblDobl_Complex_Series_Vectors.Clear(eva);
      dd_t := Create(t);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,dd_t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);
      end if;
      loop
        loop
          dd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
          predres := Residual_Prediction(wrk_sol,t);
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
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (file,abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
         -- (file,nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
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
      (file,abh,1.0,tolres,pars.corsteps,nbrit,
       wrk_sol,err,rco,res,fail,verbose);
     -- (file,nbq,1.0,tolres,pars.corsteps,nbrit,
     --  wrk_sol,err,rco,res,fail,verbose);
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
                verbose : in boolean := false ) is

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
    maxit : constant natural32 := 500;
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
        if QuadDobl_Complex_Numbers.REAL_PART(cfp) >= 0.0
         then put(file,"Closest pole :"); put(file,cfp); new_line(file);
        end if;
      end if;
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,alpha,verbose);
      step := pars.sbeta*step;
      step := Series_and_Predictors.Cap_Step_Size
                (step,hihi_part(frp),pars.pbeta);
      QuadDobl_Complex_Series_Vectors.Clear(eva);
      qd_t := Create(t);
      dstep := Series_and_Predictors.Step_Distance
                 (maxdeg,dbeta,qd_t,jm,hs,wrk_sol,srv,pv);
      step := Minimum(step,dstep);
      Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);
      end if;
      loop
        loop
          qd_step := create(step);
          wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
          predres := Residual_Prediction(wrk_sol,t);
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
        Update_Step_Sizes(minsize,maxsize,step);
        exit when ((step < pars.minsize) and (predres > alpha));
        Homotopy_Newton_Steps.Correct
          (file,abh,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
         -- (file,nbq,t,tolres,maxit,nbrit,wrk_sol,err,rco,res,fail,verbose);
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
      (file,abh,1.0,tolres,pars.corsteps,nbrit,
       wrk_sol,err,rco,res,fail,verbose);
     -- (file,nbq,1.0,tolres,pars.corsteps,nbrit,
     --  wrk_sol,err,rco,res,fail,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    qd_t := Quad_Double_Numbers.Create(t);
    sol.t := QuadDobl_Complex_Numbers.Create(qd_t);
    sol.v := wrk_sol;
    sol.err := Quad_Double_Numbers.Create(err);
    sol.rco := Quad_Double_Numbers.Create(rco);
    sol.res := Quad_Double_Numbers.Create(res);
    QuadDobl_CSeries_Vector_Functions.Deep_Clear(wrk_fcf);
  end Track_One_Path;

  procedure Update_Counters
              ( mincnt,maxcnt : in out natural32; cnt : in natural32 ) is
  begin
    if cnt < mincnt then
      mincnt := cnt;
    elsif cnt > maxcnt then
      maxcnt := cnt;
    end if;
  end Update_Counters;

  procedure Update_MinMax
              ( smallest,largest : in out double_float;
                minsize,maxsize : in double_float ) is
  begin
    if minsize < smallest
     then smallest := minsize;
    end if;
    if maxsize > largest
     then largest := maxsize;
    end if;
  end Update_MinMax;

  procedure Track_Many_Paths
              ( file : in file_type;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                monitor,verbose : in boolean := false ) is

    use Standard_Complex_Solutions;

    abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(hom'range)
        := Homotopy_Mixed_Residuals.Standard_AbsVal_Homotopy;
    nvr : constant integer32 := Head_Of(sols).n;
    fhm : Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(hom'range)
        := Standard_CSeries_Poly_SysFun.Create(hom);
    fcf : Standard_Complex_Series_VecVecs.VecVec(hom'range)
        := Standard_CSeries_Poly_SysFun.Coeff(hom);
    ejm : Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(hom'range,1..nvr);
    mlt : Standard_CSeries_Jaco_Matrices.Mult_Factors(hom'range,1..nvr);
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs,cntfail : natural32;
    minsize,maxsize,smallest,largest : double_float;

  begin
    Standard_CSeries_Jaco_Matrices.Create(hom,ejm,mlt);
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (pars.maxsteps+1)*pars.corsteps+1; maxnbrcorrs := 0;
    smallest := pars.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      end if;
      Track_One_Path(file,abh,jm,hs,fhm,fcf,ejm,mlt,ls.all,pars,
                     nbrsteps,nbrcorrs,cntfail,minsize,maxsize,verbose);
      if verbose then
        Write_Path_Statistics(file,nbrsteps,nbrcorrs,cntfail,minsize,maxsize);
      end if;
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      Standard_Complex_Solutions_io.put(file,ls.all); new_line(file);
      tmp := Tail_Of(tmp);
      Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Update_MinMax(smallest,largest,minsize,maxsize);
    end loop;
    tstop(timer);
    Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,smallest,largest);
    new_line(file);
    print_times(file,timer,"Tracking in double precision.");
    Standard_CSeries_Poly_SysFun.Clear(fhm);
    Standard_Complex_Series_VecVecs.Clear(fcf);
    Standard_CSeries_Jaco_Matrices.Clear(ejm);
    Standard_CSeries_Jaco_Matrices.Clear(mlt);
    Standard_Complex_Poly_SysFun.Clear(abh);
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( file : in file_type;
                jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                monitor,verbose : in boolean := false ) is

    use DoblDobl_Complex_Solutions;

    abh : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(hom'range)
        := Homotopy_Mixed_Residuals.DoblDobl_AbsVal_Homotopy;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs,cntfail : natural32;
    minsize,maxsize,smallest,largest : double_float;

  begin
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (pars.maxsteps+1)*pars.corsteps+1; maxnbrcorrs := 0;
    smallest := pars.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      end if;
      Track_One_Path(file,abh,jm,hs,hom,ls.all,pars,nbrsteps,nbrcorrs,cntfail,
                     minsize,maxsize,verbose);
      if verbose then
        Write_Path_Statistics(file,nbrsteps,nbrcorrs,cntfail,minsize,maxsize);
      end if;
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      DoblDobl_Complex_Solutions_io.put(file,ls.all); new_line(file);
      tmp := Tail_Of(tmp);
      Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Update_MinMax(smallest,largest,minsize,maxsize);
    end loop;
    tstop(timer);
    Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,smallest,largest);
    new_line(file);
    print_times(file,timer,"Tracking in double double precision.");
    DoblDobl_Complex_Poly_SysFun.Clear(abh);
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( file : in file_type;
                jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                monitor,verbose : in boolean := false ) is

    use QuadDobl_Complex_Solutions;

    abh : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(hom'range)
        := Homotopy_Mixed_Residuals.QuadDobl_AbsVal_Homotopy;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs,cntfail : natural32;
    minsize,maxsize,smallest,largest : double_float;

  begin
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (pars.maxsteps+1)*pars.corsteps+1; maxnbrcorrs := 0;
    smallest := pars.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      end if;
      Track_One_Path(file,abh,jm,hs,hom,ls.all,pars,nbrsteps,nbrcorrs,cntfail,
                     minsize,maxsize,verbose);
      if verbose then
        Write_Path_Statistics(file,nbrsteps,nbrcorrs,cntfail,minsize,maxsize);
      end if;
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      QuadDobl_Complex_Solutions_io.put(file,ls.all); new_line(file);
      tmp := Tail_Of(tmp);
      Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Update_MinMax(smallest,largest,minsize,maxsize);
    end loop;
    tstop(timer);
    Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,smallest,largest);
    new_line(file);
    print_times(file,timer,"Tracking in quad double precision.");
    QuadDobl_Complex_Poly_SysFun.Clear(abh);
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

    use Standard_Complex_Solutions;

    abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(hom'range)
        := Homotopy_Mixed_Residuals.Standard_AbsVal_Homotopy;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    nbrsteps,nbrcorrs,cntfail : natural32;
    minsize,maxsize : double_float;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Track_One_Path
        (abh,jm,hs,hom,ls.all,pars,nbrsteps,nbrcorrs,cntfail,minsize,maxsize);
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Complex_Poly_SysFun.Clear(abh);
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

    use DoblDobl_Complex_Solutions;

    abh : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(hom'range)
        := Homotopy_Mixed_Residuals.DoblDobl_AbsVal_Homotopy;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    nbrsteps,nbrcorrs,cntfail : natural32;
    minsize,maxsize : double_float;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Track_One_Path
        (abh,jm,hs,hom,ls.all,pars,nbrsteps,nbrcorrs,cntfail,minsize,maxsize);
      tmp := Tail_Of(tmp);
    end loop;
    DoblDobl_Complex_Poly_SysFun.Clear(abh);
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

    use QuadDobl_Complex_Solutions;

    abh : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(hom'range)
        := Homotopy_Mixed_Residuals.QuadDobl_AbsVal_Homotopy;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    nbrsteps,nbrcorrs,cntfail : natural32;
    minsize,maxsize : double_float;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Track_One_Path
        (abh,jm,hs,hom,ls.all,pars,nbrsteps,nbrcorrs,cntfail,minsize,maxsize);
      tmp := Tail_Of(tmp);
    end loop;
    QuadDobl_Complex_Poly_SysFun.Clear(abh);
  end Track_Many_Paths;

  procedure Write_Path_Statistics
              ( file : in file_type;
                nbrsteps,nbrcorrs,cntfail : in natural32;
                minsize,maxsize : in double_float ) is
  begin
    put(file,"The total number of steps on the path     : ");
    put(file,nbrsteps,1); new_line(file);
    put(file,"Total number of correct steps on the path : ");
    put(file,nbrcorrs,1); new_line(file);
    put(file,"Number of corrector failures on the path  : ");
    put(file,cntfail,1); new_line(file);
    put(file,"The smallest step size on the path        :");
    put(file,minsize,2); new_line(file);
    put(file,"The largest step size on the path         :");
    put(file,maxsize,2); new_line(file);
  end Write_Path_Statistics;

  procedure Write_Total_Path_Statistics
              ( file : in file_type;
                minnbrsteps,maxnbrsteps : in natural32;
                minnbrcorrs,maxnbrcorrs : in natural32;
                smallestsize,largestsize : in double_float ) is
  begin
    new_line(file);
    put(file,"The smallest number of total steps : ");
    put(file,minnbrsteps,1); new_line(file);
    put(file,"The largest number of total steps  : ");
    put(file,maxnbrsteps,1); new_line(file);
    put(file,"The smallest number of corrector iterations : ");
    put(file,minnbrcorrs,1); new_line(file);
    put(file,"The largest number of corrector iterations  : ");
    put(file,maxnbrcorrs,1); new_line(file);
    put(file,"The smallest step size on a path :");
    put(file,smallestsize,2); new_line(file);
    put(file,"The largest step size on a path  :");
    put(file,largestsize,2); new_line(file);
  end Write_Total_Path_Statistics;

end Series_and_Trackers;
