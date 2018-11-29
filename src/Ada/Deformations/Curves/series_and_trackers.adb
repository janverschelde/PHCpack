with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;
with Standard_Root_Refiners;
with DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;
with Standard_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with Standard_Pade_Approximants;
with DoblDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants;
with Series_and_Homotopies;
with Series_and_Predictors;

package body Series_and_Trackers is

  procedure Correct
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                t : in double_float; tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    p : Standard_Complex_Poly_Systems.Poly_Sys(hom'range)
      := Series_and_Homotopies.Eval(hom,t);
    jm : Standard_Complex_Jaco_Matrices.Jaco_Mat(hom'range,sol'range)
       := Standard_Complex_Jaco_Matrices.Create(p);
    f : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := Standard_Complex_Poly_SysFun.Create(p);
    jf : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat(hom'range,sol'range)
       := Standard_Complex_Jaco_Matrices.Create(jm);

    use Standard_Root_Refiners;

  begin
    nbrit := maxit;
    for k in 1..maxit loop
      Standard_Newton_Step(f,jf,sol,err,rco,res);
      if res <= tolres
       then nbrit := k; exit;
      end if;
    end loop;
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Poly_SysFun.Clear(f);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( file : in file_type;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                t : in double_float; tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float;
                verbose : in boolean := false ) is

    p : Standard_Complex_Poly_Systems.Poly_Sys(hom'range)
      := Series_and_Homotopies.Eval(hom,t);
    jm : Standard_Complex_Jaco_Matrices.Jaco_Mat(hom'range,sol'range)
       := Standard_Complex_Jaco_Matrices.Create(p);
    f : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := Standard_Complex_Poly_SysFun.Create(p);
    jf : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat(hom'range,sol'range)
       := Standard_Complex_Jaco_Matrices.Create(jm);

    use Standard_Root_Refiners;

  begin
    nbrit := maxit;
    for k in 1..maxit loop
      Standard_Newton_Step(f,jf,sol,err,rco,res);
      if verbose then
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
      if res <= tolres
       then nbrit := k; exit;
      end if;
    end loop;
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Poly_SysFun.Clear(f);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in double_double; tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double ) is

    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(hom'range)
      := Series_and_Homotopies.Eval(hom,t);
    jm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(hom'range,sol'range)
       := DoblDobl_Complex_Jaco_Matrices.Create(p);
    f : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := DoblDobl_Complex_Poly_SysFun.Create(p);
    jf : DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat(hom'range,sol'range)
       := DoblDobl_Complex_Jaco_Matrices.Create(jm);

    use DoblDobl_Root_Refiners;

  begin
    nbrit := maxit;
    for k in 1..maxit loop
      DoblDobl_Newton_Step(f,jf,sol,err,rco,res);
      if res <= tolres
       then nbrit := k; exit;
      end if;
    end loop;
    DoblDobl_Complex_Poly_Systems.Clear(p);
    DoblDobl_Complex_Poly_SysFun.Clear(f);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( file : in file_type;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in double_double; tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in boolean := false ) is

    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(hom'range)
      := Series_and_Homotopies.Eval(hom,t);
    jm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(hom'range,sol'range)
       := DoblDobl_Complex_Jaco_Matrices.Create(p);
    f : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := DoblDobl_Complex_Poly_SysFun.Create(p);
    jf : DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat(hom'range,sol'range)
       := DoblDobl_Complex_Jaco_Matrices.Create(jm);

    use DoblDobl_Root_Refiners;

  begin
    nbrit := maxit;
    for k in 1..maxit loop
      DoblDobl_Newton_Step(f,jf,sol,err,rco,res);
      if verbose then
        put(file,"  err : "); put(file,err,3);
        put(file,"  rco : "); put(file,rco,3);
        put(file,"  res : "); put(file,res,3); new_line(file);
      end if;
      if res <= tolres
       then nbrit := k; exit;
      end if;
    end loop;
    DoblDobl_Complex_Poly_Systems.Clear(p);
    DoblDobl_Complex_Poly_SysFun.Clear(f);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in quad_double; tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double ) is

    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(hom'range)
      := Series_and_Homotopies.Eval(hom,t);
    jm : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(hom'range,sol'range)
       := QuadDobl_Complex_Jaco_Matrices.Create(p);
    f : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := QuadDobl_Complex_Poly_SysFun.Create(p);
    jf : QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat(hom'range,sol'range)
       := QuadDobl_Complex_Jaco_Matrices.Create(jm);

    use QuadDobl_Root_Refiners;

  begin
    nbrit := maxit;
    for k in 1..nbrit loop
      QuadDobl_Newton_Step(f,jf,sol,err,rco,res);
      if res <= tolres
       then nbrit := k; exit;
      end if;
    end loop;
    QuadDobl_Complex_Poly_Systems.Clear(p);
    QuadDobl_Complex_Poly_SysFun.Clear(f);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( file : in file_type; 
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in quad_double; tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in boolean := false ) is

    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(hom'range)
      := Series_and_Homotopies.Eval(hom,t);
    jm : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(hom'range,sol'range)
       := QuadDobl_Complex_Jaco_Matrices.Create(p);
    f : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := QuadDobl_Complex_Poly_SysFun.Create(p);
    jf : QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat(hom'range,sol'range)
       := QuadDobl_Complex_Jaco_Matrices.Create(jm);

    use QuadDobl_Root_Refiners;

  begin
    nbrit := maxit;
    for k in 1..maxit loop
      QuadDobl_Newton_Step(f,jf,sol,err,rco,res);
      if verbose then
        put(file,"  err : "); put(file,err,3);
        put(file,"  rco : "); put(file,rco,3);
        put(file,"  res : "); put(file,res,3); new_line(file);
      end if;
      if res <= tolres
       then nbrit := k; exit;
      end if;
    end loop;
    QuadDobl_Complex_Poly_Systems.Clear(p);
    QuadDobl_Complex_Poly_SysFun.Clear(f);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

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
              ( hom : Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : Standard_Complex_Vectors.Vector;
                step : double_float ) return double_float is

  -- DESCRIPTION :
  --   Given a homotopy, a predicted solution, with a step,
  --   returns the norm of the evaluated homotopy at the solution.

    res : double_float;
    phm : Standard_Complex_Poly_Systems.Poly_Sys(hom'range);
    val : Standard_Complex_Vectors.Vector(phm'range);

  begin
    phm := Series_and_Homotopies.Eval(hom,step);
    val := Standard_Complex_Poly_SysFun.Eval(phm,sol);
    res := Standard_Complex_Vector_Norms.Max_Norm(val);
    Standard_Complex_Poly_Systems.Clear(phm);
    return res;
  end Residual_Prediction;

  function Residual_Prediction
              ( hom : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : DoblDobl_Complex_Vectors.Vector;
                step : double_float ) return double_float is

  -- DESCRIPTION :
  --   Given a homotopy, a predicted solution, with a step,
  --   returns the norm of the evaluated homotopy at the solution.

    res : double_double;
    phm : DoblDobl_Complex_Poly_Systems.Poly_Sys(hom'range);
    val : DoblDobl_Complex_Vectors.Vector(phm'range);
    dd_step : constant double_double := create(step);

  begin
    phm := Series_and_Homotopies.Eval(hom,dd_step);
    val := DoblDobl_Complex_Poly_SysFun.Eval(phm,sol);
    res := DoblDobl_Complex_Vector_Norms.Max_Norm(val);
    DoblDobl_Complex_Poly_Systems.Clear(phm);
    return hi_part(res);
  end Residual_Prediction;

  function Residual_Prediction
              ( hom : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : QuadDobl_Complex_Vectors.Vector;
                step : double_float ) return double_float is

  -- DESCRIPTION :
  --   Given a homotopy, a predicted solution, with a step,
  --   returns the norm of the evaluated homotopy at the solution.

    res : quad_double;
    phm : QuadDobl_Complex_Poly_Systems.Poly_Sys(hom'range);
    val : QuadDobl_Complex_Vectors.Vector(phm'range);
    qd_step : constant quad_double := create(step);

  begin
    phm := Series_and_Homotopies.Eval(hom,qd_step);
    val := QuadDobl_Complex_Poly_SysFun.Eval(phm,sol);
    res := QuadDobl_Complex_Vector_Norms.Max_Norm(val);
    QuadDobl_Complex_Poly_Systems.Clear(phm);
    return hihi_part(res);
  end Residual_Prediction;

  procedure Track_One_Path
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs : out natural32;
                minsize,maxsize : out double_float ) is

    wrk : Standard_CSeries_Poly_Systems.Poly_Sys(hom'range);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := integer32(pars.corsteps);
    srv : Standard_Complex_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Complex_Series_Vectors.Vector(hom'range);
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range);
    poles : Standard_Complex_VecVecs.VecVec(pv'range);
    tolcff : constant double_float := pars.tolcff;
    tolres : constant double_float := pars.alpha;
    t,step,update : double_float := 0.0;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res,frp,predres : double_float;
    nbrit : natural32 := 0;

  begin
    minsize := 1.0; maxsize := 0.0;
    Standard_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction(maxdeg,nit,wrk,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(numdeg,dendeg,srv,pv,poles,frp);
      Standard_Complex_VecVecs.Clear(poles);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,tolres);
      step := pars.sbeta*step;
      Standard_Complex_Series_Vectors.Clear(eva);
      if frp > 0.0
       then step := Series_and_Predictors.Cap_Step_Size(step,frp,pars.pbeta);
      end if;
      Set_Step(t,step,pars.maxsize,onetarget);
      exit when (step < pars.minsize);
      loop
        wrk_sol := Series_and_Predictors.Predicted_Solution(pv,step);
        predres := Residual_Prediction(wrk,wrk_sol,step);
        exit when (predres <= pars.alpha);
        t := t - step; step := step/2.0; t := t + step;
        exit when (step < pars.minsize);
      end loop;
      Update_Step_Sizes(minsize,maxsize,step);
      exit when (step < pars.minsize);
      Correct(wrk,step,tolcff,pars.corsteps,nbrit,wrk_sol,err,rco,res);
      nbrcorrs := nbrcorrs + nbrit;
      Standard_Complex_Series_Vectors.Clear(srv);
      Standard_Pade_Approximants.Clear(pv);
      Standard_CSeries_Poly_Systems.Clear(wrk);
      if t = 1.0
       then nbrsteps := k; exit;
      end if;
      wrk := Series_and_Homotopies.Shift(hom,-t);
    end loop;
    wrk := Series_and_Homotopies.Shift(hom,-1.0);
    Correct(wrk,0.0,tolcff,pars.corsteps,nbrit,wrk_sol,err,rco,res);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := Standard_Complex_Numbers.Create(t);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs : out natural32;
                minsize,maxsize : out double_float ) is

    wrk : DoblDobl_CSeries_Poly_Systems.Poly_Sys(hom'range);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := integer32(pars.corsteps);
    srv : DoblDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : DoblDobl_Complex_Series_Vectors.Vector(hom'range);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range);
    tolcff : constant double_float := pars.tolcff;
    tolres : constant double_float := pars.alpha;
    t,step,update : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : double_double;
    zero : constant double_double := create(0.0);
    frp : double_double;
    predres : double_float;
    nbrit : natural32 := 0;

  begin
    minsize := 1.0; maxsize := 0.0;
    DoblDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction(maxdeg,nit,wrk,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(numdeg,dendeg,srv,pv,poles,frp);
      DoblDobl_Complex_VecVecs.Clear(poles);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,tolres);
      step := pars.sbeta*step;
      if frp > 0.0 then
        step := Series_and_Predictors.Cap_Step_Size
                  (step,hi_part(frp),pars.pbeta);
      end if;
      DoblDobl_Complex_Series_Vectors.Clear(eva);
      Set_Step(t,step,pars.maxsize,onetarget);
      loop
        dd_step := create(step);
        wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
        predres := Residual_Prediction(wrk,wrk_sol,step);
        exit when (predres <= pars.alpha);
        t := t - step; step := step/2.0; t := t + step;
        exit when (step < pars.minsize);
      end loop;
      Update_Step_Sizes(minsize,maxsize,step);
      exit when (step < pars.minsize);
      Correct(wrk,dd_step,tolcff,pars.corsteps,nbrit,wrk_sol,err,rco,res);
      nbrcorrs := nbrcorrs + nbrit;
      DoblDobl_Complex_Series_Vectors.Clear(srv);
      DoblDobl_Pade_Approximants.Clear(pv);
      DoblDobl_CSeries_Poly_Systems.Clear(wrk);
      dd_t := create(-t);
      wrk := Series_and_Homotopies.Shift(hom,dd_t);
      if t = 1.0
       then nbrsteps := k; exit;
      end if;
    end loop;
    dd_t := create(-1.0);
    wrk := Series_and_Homotopies.Shift(hom,dd_t);
    dd_step := create(0.0);
    Correct(wrk,dd_step,tolcff,pars.corsteps,nbrit,wrk_sol,err,rco,res);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := DoblDobl_Complex_Numbers.Create(Double_Double_Numbers.Create(-t));
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    DoblDobl_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs : out natural32;
                minsize,maxsize : out double_float ) is

    wrk : QuadDobl_CSeries_Poly_Systems.Poly_Sys(hom'range);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := integer32(pars.corsteps);
    srv : QuadDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : QuadDobl_Complex_Series_Vectors.Vector(hom'range);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range);
    tolcff : constant double_float := pars.tolcff;
    tolres : constant double_float := pars.alpha;
    t,step,update : double_float := 0.0;
    qd_t,qd_step : quad_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : QuadDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : quad_double;
    zero : constant quad_double := create(0.0);
    frp : quad_double;
    predres : double_float;
    nbrit : natural32 := 0;

  begin
    minsize := 1.0; maxsize := 0.0;
    QuadDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      Series_and_Predictors.Newton_Prediction(maxdeg,nit,wrk,wrk_sol,srv,eva);
      Series_and_Predictors.Pade_Approximants(numdeg,dendeg,srv,pv,poles,frp);
      QuadDobl_Complex_VecVecs.Clear(poles);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,tolres);
      step := pars.sbeta*step;
      if frp > 0.0 then
        step := Series_and_Predictors.Cap_Step_Size
                  (step,hihi_part(frp),pars.pbeta);
      end if;
      QuadDobl_Complex_Series_Vectors.Clear(eva);
      Set_Step(t,step,pars.maxsize,onetarget);
      loop
        qd_step := create(step);
        wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
        predres := Residual_Prediction(wrk,wrk_sol,step);
        exit when (predres <= pars.alpha);
        t := t - step; step := step/2.0; t := t + step;
        exit when (step < pars.minsize);
      end loop;
      Update_Step_Sizes(minsize,maxsize,step);
      exit when (step < pars.minsize);
      Correct(wrk,qd_step,tolcff,pars.corsteps,nbrit,wrk_sol,err,rco,res);
      nbrcorrs := nbrcorrs + nbrit;
      QuadDobl_Pade_Approximants.Clear(pv);
      QuadDobl_Complex_Series_Vectors.Clear(srv);
      QuadDobl_CSeries_Poly_Systems.Clear(wrk);
      qd_t := create(-t);
      wrk := Series_and_Homotopies.Shift(hom,qd_t);
      if t = 1.0
       then nbrsteps := k; exit;
      end if;
    end loop;
    qd_t := create(-1.0);
    wrk := Series_and_Homotopies.Shift(hom,qd_t);
    qd_step := create(0.0);
    Correct(wrk,qd_step,tolcff,pars.corsteps,nbrit,wrk_sol,err,rco,res);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := QuadDobl_Complex_Numbers.Create(Quad_Double_Numbers.Create(-t));
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    QuadDobl_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false ) is

    wrk : Standard_CSeries_Poly_Systems.Poly_Sys(hom'range);
    nit : constant integer32 := integer32(pars.corsteps);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    srv : Standard_Complex_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Complex_Series_Vectors.Vector(hom'range);
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range);
    poles : Standard_Complex_VecVecs.VecVec(pv'range);
    tolcff : constant double_float := pars.tolcff;
    tolres : constant double_float := pars.alpha;
    t,step,update : double_float := 0.0;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res,frp,predres : double_float;
    nbrit : natural32 := 0;

  begin
    minsize := 1.0; maxsize := 0.0;
    Standard_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      if verbose then
        put(file,"Step "); put(file,k,1); put(file," : ");
      end if;
      Series_and_Predictors.Newton_Prediction
        (file,maxdeg,nit,wrk,wrk_sol,srv,eva,false); -- verbose);
      Series_and_Predictors.Pade_Approximants(numdeg,dendeg,srv,pv,poles,frp);
      if verbose then
        put(file,"The smallest forward pole radius : ");
        put(file,frp,3); new_line(file);
      end if;
      Standard_Complex_VecVecs.Clear(poles);
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,tolres,verbose);
      step := pars.sbeta*step;
      if frp > 0.0
       then step := Series_and_Predictors.Cap_Step_Size(step,frp,pars.pbeta);
      end if;
      Standard_Complex_Series_Vectors.Clear(eva);
      Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);
      end if;
      loop
        wrk_sol := Series_and_Predictors.Predicted_Solution(pv,step);
        predres := Residual_Prediction(wrk,wrk_sol,step);
        if verbose
         then put(file,"  residual : "); put(file,predres,3); new_line(file);
        end if;
        exit when (predres <= pars.alpha);
        t := t - step; step := step/2.0; t := t + step;
        if verbose then
          put(file,"Step size : "); put(file,step,3);
          put(file," t = "); put(file,t,3);
        end if;
        exit when (step < pars.minsize);
      end loop;
      Update_Step_Sizes(minsize,maxsize,step);
      exit when (step < pars.minsize);
      Correct(file,wrk,step,tolcff,pars.corsteps,nbrit,
              wrk_sol,err,rco,res,verbose);
      nbrcorrs := nbrcorrs + nbrit;
      Standard_Pade_Approximants.Clear(pv);
      Standard_Complex_Series_Vectors.Clear(srv);
      Standard_CSeries_Poly_Systems.Clear(wrk);
      if t = 1.0
       then nbrsteps := k; exit;
      end if;
      wrk := Series_and_Homotopies.Shift(hom,-t);
    end loop;
    wrk := Series_and_Homotopies.Shift(hom,-1.0);
    Correct(file,wrk,0.0,tolcff,pars.corsteps,nbrit,
            wrk_sol,err,rco,res,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := Standard_Complex_Numbers.Create(t);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false ) is

    wrk : DoblDobl_CSeries_Poly_Systems.Poly_Sys(hom'range);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := integer32(pars.corsteps);
    srv : DoblDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : DoblDobl_Complex_Series_Vectors.Vector(hom'range);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range);
    tolcff : constant double_float := pars.tolcff;
    tolres : constant double_float := pars.alpha;
    t,step,update : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : double_double;
    zero : constant double_double := create(0.0);
    frp : double_double;
    predres : double_float;
    nbrit : natural32 := 0;

  begin
    minsize := 1.0; maxsize := 0.0;
    DoblDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      if verbose then
        put(file,"Step "); put(file,k,1); put(file," : ");
      end if;
      Series_and_Predictors.Newton_Prediction
        (file,maxdeg,nit,wrk,wrk_sol,srv,eva,false); -- verbose);
      Series_and_Predictors.Pade_Approximants(numdeg,dendeg,srv,pv,poles,frp);
      DoblDobl_Complex_VecVecs.Clear(poles);
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,tolres,verbose);
      step := pars.sbeta*step;
      if frp > 0.0 then
        step := Series_and_Predictors.Cap_Step_Size
                  (step,hi_part(frp),pars.pbeta);
      end if;
      DoblDobl_Complex_Series_Vectors.Clear(eva);
      Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);  new_line(file);
      end if;
      loop
        dd_step := create(step);
        wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
        predres := Residual_Prediction(wrk,wrk_sol,step);
        if verbose
         then put(file,"  residual : "); put(file,predres,3); new_line(file);
        end if;
        exit when (predres <= pars.alpha);
        t := t - step; step := step/2.0; t := t + step;
        if verbose then
          put(file,"Step size : "); put(file,step,3);
          put(file," t = "); put(file,t,3);
        end if;
        exit when (step < pars.minsize);
      end loop;
      Update_Step_Sizes(minsize,maxsize,step);
      exit when (step < pars.minsize);
      Correct(file,wrk,dd_step,tolcff,pars.corsteps,nbrit,
              wrk_sol,err,rco,res,verbose);
      nbrcorrs := nbrcorrs + nbrit;
      DoblDobl_Complex_Series_Vectors.Clear(srv);
      DoblDobl_Pade_Approximants.Clear(pv);
      DoblDobl_CSeries_Poly_Systems.Clear(wrk);
      dd_t := create(-t);
      wrk := Series_and_Homotopies.Shift(hom,dd_t);
      if t = 1.0
       then nbrsteps := k; exit;
      end if;
    end loop;
    dd_t := create(-1.0);
    wrk := Series_and_Homotopies.Shift(hom,dd_t);
    dd_step := create(0.0);
    Correct(file,wrk,dd_step,tolcff,pars.corsteps,nbrit,
            wrk_sol,err,rco,res,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := DoblDobl_Complex_Numbers.Create(Double_Double_Numbers.Create(-t));
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    DoblDobl_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs : out natural32;
                minsize,maxsize : out double_float;
                verbose : in boolean := false ) is

    wrk : QuadDobl_CSeries_Poly_Systems.Poly_Sys(hom'range);
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := integer32(pars.corsteps);
    srv : QuadDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : QuadDobl_Complex_Series_Vectors.Vector(hom'range);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range);
    tolcff : constant double_float := pars.tolcff;
    tolres : constant double_float := pars.alpha;
    t,step,update : double_float := 0.0;
    qd_t,qd_step : quad_double;
    max_steps : constant natural32 := pars.maxsteps;
    wrk_sol : QuadDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : quad_double;
    zero : constant quad_double := create(0.0);
    frp : quad_double;
    predres : double_float;
    nbrit : natural32 := 0;

  begin
    minsize := 1.0; maxsize := 0.0;
    QuadDobl_CSeries_Poly_Systems.Copy(hom,wrk);
    nbrcorrs := 0;
    nbrsteps := max_steps;
    for k in 1..max_steps loop
      if verbose then
        put(file,"Step "); put(file,k,1); put(file," : ");
      end if;
      Series_and_Predictors.Newton_Prediction
        (file,maxdeg,nit,wrk,wrk_sol,srv,eva,false); -- verbose);
      Series_and_Predictors.Pade_Approximants(numdeg,dendeg,srv,pv,poles,frp);
      QuadDobl_Complex_VecVecs.Clear(poles);
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,tolres,verbose);
      step := pars.sbeta*step;
      if frp > 0.0 then
        step := Series_and_Predictors.Cap_Step_Size
                  (step,hihi_part(frp),pars.pbeta);
      end if;
      QuadDobl_Complex_Series_Vectors.Clear(eva);
      Set_Step(t,step,pars.maxsize,onetarget);
      if verbose then
        put(file,"Step size : "); put(file,step,3);
        put(file," t = "); put(file,t,3);  new_line(file);
      end if;
      loop
        qd_step := create(step);
        wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
        predres := Residual_Prediction(wrk,wrk_sol,step);
        if verbose
         then put(file,"  residual : "); put(file,predres,3); new_line(file);
        end if;
        exit when (predres <= pars.alpha);
        t := t - step; step := step/2.0; t := t + step;
        if verbose then
          put(file,"Step size : "); put(file,step,3);
          put(file," t = "); put(file,t,3);
        end if;
        exit when (step < pars.minsize);
      end loop;
      Update_Step_Sizes(minsize,maxsize,step);
      exit when (step < pars.minsize);
      Correct(file,wrk,qd_step,tolcff,pars.corsteps,nbrit,
              wrk_sol,err,rco,res,verbose);
      nbrcorrs := nbrcorrs + nbrit;
      QuadDobl_Complex_Series_Vectors.Clear(srv);
      QuadDobl_Pade_Approximants.Clear(pv);
      QuadDobl_CSeries_Poly_Systems.Clear(wrk);
      qd_t := create(-t);
      wrk := Series_and_Homotopies.Shift(hom,qd_t);
      if t = 1.0
       then nbrsteps := k; exit;
      end if;
    end loop;
    qd_t := create(-1.0);
    wrk := Series_and_Homotopies.Shift(hom,qd_t);
    qd_step := create(0.0);
    Correct(file,wrk,qd_step,tolcff,pars.corsteps,nbrit,
            wrk_sol,err,rco,res,verbose);
    nbrcorrs := nbrcorrs + nbrit;
    sol.t := QuadDobl_Complex_Numbers.Create(Quad_Double_Numbers.Create(-t));
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    QuadDobl_CSeries_Poly_Systems.Clear(wrk);
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

  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false ) is

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs : natural32;
    minsize,maxsize : double_float;

  begin
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := pars.corsteps+1; maxnbrcorrs := 0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      Track_One_Path
        (file,hom,ls.all,pars,nbrsteps,nbrcorrs,minsize,maxsize,verbose);
      if verbose
       then Write_Path_Statistics(file,nbrsteps,nbrcorrs,minsize,maxsize);
      end if;
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      Standard_Complex_Solutions_io.put(file,ls.all);
      tmp := Tail_Of(tmp);
      Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
    end loop;
    tstop(timer);
    Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs);
    new_line(file);
    print_times(file,timer,"Tracking in double precision.");
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false ) is

    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs : natural32;
    minsize,maxsize : double_float;

  begin
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := pars.corsteps+1; maxnbrcorrs := 0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      Track_One_Path
        (file,hom,ls.all,pars,nbrsteps,nbrcorrs,minsize,maxsize,verbose);
      if verbose
       then Write_Path_Statistics(file,nbrsteps,nbrcorrs,minsize,maxsize);
      end if;
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      DoblDobl_Complex_Solutions_io.put(file,ls.all);
      tmp := Tail_Of(tmp);
      Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
    end loop;
    tstop(timer);
    Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs);
    new_line(file);
    print_times(file,timer,"Tracking in double double precision.");
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := false ) is

    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    timer : Timing_Widget;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs : natural32;
    minsize,maxsize : double_float;

  begin
    minnbrsteps := pars.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := pars.corsteps+1; maxnbrcorrs := 0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      Track_One_Path
        (file,hom,ls.all,pars,nbrsteps,nbrcorrs,minsize,maxsize,verbose);
      if verbose
       then Write_Path_Statistics(file,nbrsteps,nbrcorrs,minsize,maxsize);
      end if;
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      QuadDobl_Complex_Solutions_io.put(file,ls.all);
      tmp := Tail_Of(tmp);
      Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
    end loop;
    tstop(timer);
    Write_Total_Path_Statistics
      (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs);
    new_line(file);
    print_times(file,timer,"Tracking in quad double precision.");
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    nbrsteps,nbrcorrs : natural32;
    minsize,maxsize : double_float;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Track_One_Path(hom,ls.all,pars,nbrsteps,nbrcorrs,minsize,maxsize);
      tmp := Tail_Of(tmp);
    end loop;
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    nbrsteps,nbrcorrs : natural32;
    minsize,maxsize : double_float;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Track_One_Path(hom,ls.all,pars,nbrsteps,nbrcorrs,minsize,maxsize);
      tmp := Tail_Of(tmp);
    end loop;
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    nbrsteps,nbrcorrs : natural32;
    minsize,maxsize : double_float;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Track_One_Path(hom,ls.all,pars,nbrsteps,nbrcorrs,minsize,maxsize);
      tmp := Tail_Of(tmp);
    end loop;
  end Track_Many_Paths;

  procedure Write_Path_Statistics
              ( file : in file_type;
                nbrsteps,nbrcorrs : in natural32;
                minsize,maxsize : in double_float ) is
  begin
    put(file,"The total number of steps on the path     : ");
    put(file,nbrsteps,1); new_line(file);
    put(file,"Total number of correct steps on the path : ");
    put(file,nbrcorrs,1); new_line(file);
    put(file,"The smallest step size on the path        :");
    put(file,minsize,2); new_line(file);
    put(file,"The largest step size on the path         :");
    put(file,maxsize,2); new_line(file);
  end Write_Path_Statistics;

  procedure Write_Total_Path_Statistics
              ( file : in file_type;
                minnbrsteps,maxnbrsteps : in natural32;
                minnbrcorrs,maxnbrcorrs : in natural32 ) is
  begin
    new_line(file);
    put(file,"The smallest number of total steps : ");
    put(file,minnbrsteps,1); new_line(file);
    put(file,"The largest number of total steps  : ");
    put(file,maxnbrsteps,1); new_line(file);
    put(file,"The smallest number of corrector iterations : ");
    put(file,minnbrsteps,1); new_line(file);
    put(file,"The largest number of corrector iterations  : ");
    put(file,maxnbrsteps,1); new_line(file);
  end Write_Total_Path_Statistics;

end Series_and_Trackers;
