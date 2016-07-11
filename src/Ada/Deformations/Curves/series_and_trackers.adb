with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Root_Refiners;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Dense_Series_Vectors;
with Series_and_Homotopies;
with Series_and_Predictors;

package body Series_and_Trackers is

  procedure Correct
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                t : in double_float; nit : in natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                verbose : in boolean := false ) is

    p : Standard_Complex_Poly_Systems.Poly_Sys(hom'range)
      := Series_and_Homotopies.Eval(hom,t);
    jm : Standard_Complex_Jaco_Matrices.Jaco_Mat(hom'range,sol'range)
       := Standard_Complex_Jaco_Matrices.Create(p);
    f : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := Standard_Complex_Poly_SysFun.Create(p);
    jf : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat(hom'range,sol'range)
       := Standard_Complex_Jaco_Matrices.Create(jm);
    err,rco,res : double_float;

    use Standard_Root_Refiners;

  begin
    for k in 1..nit loop
      Standard_Newton_Step(f,jf,sol,err,rco,res);
      if verbose then
        put("  err :"); put(err,3);
        put("  rco :"); put(rco,3);
        put("  res :"); put(res,3); new_line;
      end if;
    end loop;
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Poly_SysFun.Clear(f);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(jf);
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
    t,step : double_float := 0.0;
    max_steps : constant natural32 := 500;
    wrk_sol : Standard_Complex_Vectors.Vector(1..sol.n) := sol.v;

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
      if t + step <= 1.0 then
        t := t + step;
      else
        step := 1.0 - t;
        t := 1.0; 
      end if;
      put(" t = "); put(t,3);  new_line;
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,step);
      put_line("Shifting the polynomial system ...");
      Standard_Series_Poly_Systems.Clear(wrk);
      wrk := Series_and_Homotopies.Shift(hom,t);
      put_line("Correcting the solution ...");
      Correct(hom,t,3,wrk_sol,true);
      exit when (t = 1.0);
    end loop;
    put_line("The solution vector :"); put_line(wrk_sol);
    Standard_Series_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( h : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Solutions.Solution;
                verbose : in boolean := false ) is
  begin
    null;
  end Track_One_Path;

  procedure Track_One_Path
              ( h : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                s : in out Quaddobl_Complex_Solutions.Solution;
                verbose : in boolean := false ) is
  begin
    null;
  end Track_One_Path;

end Series_and_Trackers;
