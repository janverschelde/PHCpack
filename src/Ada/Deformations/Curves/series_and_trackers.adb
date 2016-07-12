with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Standard_Root_Refiners;
with DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;
with Standard_Dense_Series_Vectors;
with DoblDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_Vectors;
with Series_and_Homotopies;
with Series_and_Predictors;

package body Series_and_Trackers is

  procedure Correct
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
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
    end loop;
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Poly_SysFun.Clear(f);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( file : in file_type;
                hom : in Standard_Series_Poly_Systems.Poly_Sys;
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
    end loop;
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Poly_SysFun.Clear(f);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
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
    end loop;
    DoblDobl_Complex_Poly_Systems.Clear(p);
    DoblDobl_Complex_Poly_SysFun.Clear(f);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( file : in file_type;
                hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
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
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
    end loop;
    DoblDobl_Complex_Poly_Systems.Clear(p);
    DoblDobl_Complex_Poly_SysFun.Clear(f);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
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
    end loop;
    QuadDobl_Complex_Poly_Systems.Clear(p);
    QuadDobl_Complex_Poly_SysFun.Clear(f);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Correct
              ( file : in file_type; 
                hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
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
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
    end loop;
    QuadDobl_Complex_Poly_Systems.Clear(p);
    QuadDobl_Complex_Poly_SysFun.Clear(f);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Track_One_Path
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                verbose : in boolean := false ) is

    wrk : Standard_Series_Poly_Systems.Poly_Sys(hom'range);
    nit : integer32 := 4;
    srv : Standard_Dense_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Dense_Series_Vectors.Vector(hom'range);
    tolcff : constant double_float := 1.0E-12;
    tolres : constant double_float := 1.0E-5;
    t,step,update : double_float := 0.0;
    max_steps : constant natural32 := 500;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : double_float;

  begin
    Standard_Series_Poly_Systems.Copy(hom,wrk);
    for k in 1..max_steps loop
      put("Step "); put(k,1); put_line(" in the path tracker :");
      Series_and_Predictors.Newton_Prediction(nit,wrk,wrk_sol,srv,eva,verbose);
      new_line;
      put_line("Setting the step size based on the power series ...");
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,tolres,verbose);
      put("The computed step size : "); put(step,3);
      if step > 0.01
       then step := 0.01;
      end if;
      update := t + step;
      if update <= onetarget then
        t := update;
      else
        step := onetarget - t;
        t := onetarget; 
      end if;
      put(" t = "); put(t,3);  new_line;
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,step);
      put_line("Shifting the polynomial system ...");
      Standard_Series_Poly_Systems.Clear(wrk);
      wrk := Series_and_Homotopies.Shift(hom,t);
      put_line("Correcting the solution ...");
      Correct(standard_output,wrk,0.0,3,wrk_sol,err,rco,res,verbose);
      exit when (t = 1.0);
    end loop;
    put_line("The solution vector :"); put_line(wrk_sol);
    Standard_Series_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                verbose : in boolean := false ) is

    wrk : DoblDobl_Series_Poly_Systems.Poly_Sys(hom'range);
    nit : integer32 := 4;
    srv : DoblDobl_Dense_Series_Vectors.Vector(1..sol.n);
    eva : DoblDobl_Dense_Series_Vectors.Vector(hom'range);
    tolcff : constant double_float := 1.0E-12;
    tolres : constant double_float := 1.0E-5;
    t,step,update : double_float := 0.0;
    dd_t,dd_step : double_double;
    max_steps : constant natural32 := 500;
    wrk_sol : DoblDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : double_double;
    zero : constant double_double := create(0.0);

  begin
    DoblDobl_Series_Poly_Systems.Copy(hom,wrk);
    for k in 1..max_steps loop
      put("Step "); put(k,1); put_line(" in the path tracker :");
      Series_and_Predictors.Newton_Prediction(nit,wrk,wrk_sol,srv,eva,verbose);
      new_line;
      put_line("Setting the step size based on the power series ...");
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,tolres,verbose);
      put("The computed step size : "); put(step,3);
      if step > 0.01
       then step := 0.01;
      end if;
      update := t + step;
      if update <= onetarget then
        t := update;
      else
        step := onetarget - t;
        t := onetarget; 
      end if;
      put(" t = "); put(t,3);  new_line;
      dd_step := create(step);
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,dd_step);
      put_line("Shifting the polynomial system ...");
      DoblDobl_Series_Poly_Systems.Clear(wrk);
      dd_t := create(t);
      wrk := Series_and_Homotopies.Shift(hom,dd_t);
      put_line("Correcting the solution ...");
      Correct(standard_output,wrk,zero,3,wrk_sol,err,rco,res,verbose);
      exit when (t = 1.0);
    end loop;
    put_line("The solution vector :"); put_line(wrk_sol);
    DoblDobl_Series_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                verbose : in boolean := false ) is

    wrk : QuadDobl_Series_Poly_Systems.Poly_Sys(hom'range);
    nit : integer32 := 4;
    srv : QuadDobl_Dense_Series_Vectors.Vector(1..sol.n);
    eva : QuadDobl_Dense_Series_Vectors.Vector(hom'range);
    tolcff : constant double_float := 1.0E-12;
    tolres : constant double_float := 1.0E-5;
    t,step,update : double_float := 0.0;
    qd_t,qd_step : quad_double;
    max_steps : constant natural32 := 500;
    wrk_sol : QuadDobl_Complex_Vectors.Vector(1..sol.n) := sol.v;
    onetarget : constant double_float := 1.0;
    err,rco,res : quad_double;
    zero : constant quad_double := create(0.0);

  begin
    QuadDobl_Series_Poly_Systems.Copy(hom,wrk);
    for k in 1..max_steps loop
      put("Step "); put(k,1); put_line(" in the path tracker :");
      Series_and_Predictors.Newton_Prediction(nit,wrk,wrk_sol,srv,eva,verbose);
      new_line;
      put_line("Setting the step size based on the power series ...");
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,tolres,verbose);
      put("The computed step size : "); put(step,3);
      if step > 0.01
       then step := 0.01;
      end if;
      update := t + step;
      if update <= onetarget then
        t := update;
      else
        step := onetarget - t;
        t := onetarget; 
      end if;
      put(" t = "); put(t,3);  new_line;
      qd_step := create(step);
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,qd_step);
      put_line("Shifting the polynomial system ...");
      QuadDobl_Series_Poly_Systems.Clear(wrk);
      qd_t := create(t);
      wrk := Series_and_Homotopies.Shift(hom,qd_t);
      put_line("Correcting the solution ...");
      Correct(standard_output,wrk,zero,3,wrk_sol,err,rco,res,verbose);
      exit when (t = 1.0);
    end loop;
    put_line("The solution vector :"); put_line(wrk_sol);
    QuadDobl_Series_Poly_Systems.Clear(wrk);
  end Track_One_Path;

end Series_and_Trackers;
