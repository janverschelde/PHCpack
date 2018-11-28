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
                t : in double_float; nit : in natural32;
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
    for k in 1..nit loop
      Standard_Newton_Step(f,jf,sol,err,rco,res);
      exit when (res < 1.0e-12);
    end loop;
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Poly_SysFun.Clear(f);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( file : in file_type;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                t : in double_float; nit : in natural32;
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
    for k in 1..nit loop
      Standard_Newton_Step(f,jf,sol,err,rco,res);
      if verbose then
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
      exit when (res < 1.0e-12);
    end loop;
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Poly_SysFun.Clear(f);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in double_double; nit : in natural32;
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
    for k in 1..nit loop
      DoblDobl_Newton_Step(f,jf,sol,err,rco,res);
      exit when (res < 1.0e-12);
    end loop;
    DoblDobl_Complex_Poly_Systems.Clear(p);
    DoblDobl_Complex_Poly_SysFun.Clear(f);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( file : in file_type;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in double_double; nit : in natural32;
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
    for k in 1..nit loop
      DoblDobl_Newton_Step(f,jf,sol,err,rco,res);
      if verbose then
        put(file,"  err : "); put(file,err,3);
        put(file,"  rco : "); put(file,rco,3);
        put(file,"  res : "); put(file,res,3); new_line(file);
      end if;
      exit when (res < 1.0e-12);
    end loop;
    DoblDobl_Complex_Poly_Systems.Clear(p);
    DoblDobl_Complex_Poly_SysFun.Clear(f);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in quad_double; nit : in natural32;
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
    for k in 1..nit loop
      QuadDobl_Newton_Step(f,jf,sol,err,rco,res);
      exit when (res < 1.0e-12);
    end loop;
    QuadDobl_Complex_Poly_Systems.Clear(p);
    QuadDobl_Complex_Poly_SysFun.Clear(f);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( file : in file_type; 
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                t : in quad_double; nit : in natural32;
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
    for k in 1..nit loop
      QuadDobl_Newton_Step(f,jf,sol,err,rco,res);
      if verbose then
        put(file,"  err : "); put(file,err,3);
        put(file,"  rco : "); put(file,rco,3);
        put(file,"  res : "); put(file,res,3); new_line(file);
      end if;
      exit when (res < 1.0e-12);
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
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

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

  begin
    Standard_CSeries_Poly_Systems.Copy(hom,wrk);
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
      exit when (step < pars.minsize);
      Correct(wrk,step,pars.corsteps,wrk_sol,err,rco,res);
      Standard_Complex_Series_Vectors.Clear(srv);
      Standard_Pade_Approximants.Clear(pv);
      Standard_CSeries_Poly_Systems.Clear(wrk);
      exit when (t = 1.0);
      wrk := Series_and_Homotopies.Shift(hom,-t);
    end loop;
    wrk := Series_and_Homotopies.Shift(hom,-1.0);
    Correct(wrk,0.0,pars.corsteps,wrk_sol,err,rco,res);
    sol.t := Standard_Complex_Numbers.Create(t);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

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

  begin
    DoblDobl_CSeries_Poly_Systems.Copy(hom,wrk);
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
      exit when (step < pars.minsize);
      dd_step := create(step);
      wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
      Correct(wrk,dd_step,pars.corsteps,wrk_sol,err,rco,res);
      DoblDobl_Complex_Series_Vectors.Clear(srv);
      DoblDobl_Pade_Approximants.Clear(pv);
      DoblDobl_CSeries_Poly_Systems.Clear(wrk);
      dd_t := create(-t);
      wrk := Series_and_Homotopies.Shift(hom,dd_t);
      exit when (t = 1.0);
    end loop;
    sol.t := DoblDobl_Complex_Numbers.Create(Double_Double_Numbers.Create(-t));
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    DoblDobl_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

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

  begin
    QuadDobl_CSeries_Poly_Systems.Copy(hom,wrk);
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
      exit when (step < pars.minsize);
      qd_step := create(step);
      wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
      Correct(wrk,qd_step,pars.corsteps,wrk_sol,err,rco,res);
      QuadDobl_Pade_Approximants.Clear(pv);
      QuadDobl_Complex_Series_Vectors.Clear(srv);
      QuadDobl_CSeries_Poly_Systems.Clear(wrk);
      qd_t := create(-t);
      wrk := Series_and_Homotopies.Shift(hom,qd_t);
      exit when (t = 1.0);
    end loop;
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

  begin
    Standard_CSeries_Poly_Systems.Copy(hom,wrk);
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
      exit when (step < pars.minsize);
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
      exit when (step < pars.minsize);
      Correct(file,wrk,step,pars.corsteps,wrk_sol,err,rco,res,verbose);
      Standard_Pade_Approximants.Clear(pv);
      Standard_Complex_Series_Vectors.Clear(srv);
      Standard_CSeries_Poly_Systems.Clear(wrk);
      exit when (t = 1.0);
      wrk := Series_and_Homotopies.Shift(hom,-t);
    end loop;
    wrk := Series_and_Homotopies.Shift(hom,-1.0);
    Correct(file,wrk,0.0,pars.corsteps,wrk_sol,err,rco,res,verbose);
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

  begin
    DoblDobl_CSeries_Poly_Systems.Copy(hom,wrk);
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
      exit when (step < pars.minsize);
      dd_step := create(step);
      wrk_sol := Series_and_Predictors.Predicted_Solution(pv,dd_step);
      Correct(file,wrk,dd_step,pars.corsteps,wrk_sol,err,rco,res,verbose);
      DoblDobl_Complex_Series_Vectors.Clear(srv);
      DoblDobl_Pade_Approximants.Clear(pv);
      DoblDobl_CSeries_Poly_Systems.Clear(wrk);
      dd_t := create(-t);
      wrk := Series_and_Homotopies.Shift(hom,dd_t);
      exit when (t = 1.0);
    end loop;
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

  begin
    QuadDobl_CSeries_Poly_Systems.Copy(hom,wrk);
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
      exit when (step < pars.minsize);
      qd_step := create(step);
      wrk_sol := Series_and_Predictors.Predicted_Solution(pv,qd_step);
      Correct(file,wrk,qd_step,pars.corsteps,wrk_sol,err,rco,res,verbose);
      QuadDobl_Complex_Series_Vectors.Clear(srv);
      QuadDobl_Pade_Approximants.Clear(pv);
      QuadDobl_CSeries_Poly_Systems.Clear(wrk);
      qd_t := create(-t);
      wrk := Series_and_Homotopies.Shift(hom,qd_t);
      exit when (t = 1.0);
    end loop;
    sol.t := QuadDobl_Complex_Numbers.Create(Quad_Double_Numbers.Create(-t));
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    QuadDobl_CSeries_Poly_Systems.Clear(wrk);
  end Track_One_Path;

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

  begin
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      Track_One_Path(file,hom,ls.all,pars,verbose);
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      Standard_Complex_Solutions_io.put(file,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
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

  begin
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      Track_One_Path(file,hom,ls.all,pars,verbose);
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      DoblDobl_Complex_Solutions_io.put(file,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
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

  begin
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      Track_One_Path(file,hom,ls.all,pars,verbose);
      put(file,"Solution "); put(file,i,1); put_line(file," :");
      QuadDobl_Complex_Solutions_io.put(file,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
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

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Track_One_Path(hom,ls.all,pars);
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

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Track_One_Path(hom,ls.all,pars);
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

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Track_One_Path(hom,ls.all,pars);
      tmp := Tail_Of(tmp);
    end loop;
  end Track_Many_Paths;

end Series_and_Trackers;
