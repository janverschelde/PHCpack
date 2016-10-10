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
        put(file,"  err : "); put(file,err,3);
        put(file,"  rco : "); put(file,rco,3);
        put(file,"  res : "); put(file,res,3); new_line(file);
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
        put(file,"  err : "); put(file,err,3);
        put(file,"  rco : "); put(file,rco,3);
        put(file,"  res : "); put(file,res,3); new_line(file);
      end if;
    end loop;
    QuadDobl_Complex_Poly_Systems.Clear(p);
    QuadDobl_Complex_Poly_SysFun.Clear(f);
    QuadDobl_Complex_Jaco_Matrices.Clear(jm);
    QuadDobl_Complex_Jaco_Matrices.Clear(jf);
  end Correct;

  procedure Track_One_Path
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution ) is

    wrk : Standard_Series_Poly_Systems.Poly_Sys(hom'range);
    nit : constant integer32 := 3;
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
      Series_and_Predictors.Newton_Prediction(nit,wrk,wrk_sol,srv,eva);
      new_line;
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,tolres);
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
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,step);
      Standard_Series_Poly_Systems.Clear(wrk);
      wrk := Series_and_Homotopies.Shift(hom,t);
      Correct(wrk,0.0,3,wrk_sol,err,rco,res);
      exit when (t = 1.0);
    end loop;
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_Series_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution ) is

    wrk : DoblDobl_Series_Poly_Systems.Poly_Sys(hom'range);
    nit : constant integer32 := 4;
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
      Series_and_Predictors.Newton_Prediction(nit,wrk,wrk_sol,srv,eva);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,tolres);
      if step > 0.1
       then step := 0.1;
      end if;
      update := t + step;
      if update <= onetarget then
        t := update;
      else
        step := onetarget - t;
        t := onetarget; 
      end if;
      dd_step := create(step);
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,dd_step);
      DoblDobl_Series_Poly_Systems.Clear(wrk);
      dd_t := create(t);
      wrk := Series_and_Homotopies.Shift(hom,dd_t);
      Correct(wrk,zero,3,wrk_sol,err,rco,res);
      exit when (t = 1.0);
    end loop;
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    DoblDobl_Series_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution ) is

    wrk : QuadDobl_Series_Poly_Systems.Poly_Sys(hom'range);
    nit : constant integer32 := 4;
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
      Series_and_Predictors.Newton_Prediction(nit,wrk,wrk_sol,srv,eva);
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,tolres);
      if step > 0.1
       then step := 0.1;
      end if;
      update := t + step;
      if update <= onetarget then
        t := update;
      else
        step := onetarget - t;
        t := onetarget; 
      end if;
      qd_step := create(step);
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,qd_step);
      QuadDobl_Series_Poly_Systems.Clear(wrk);
      qd_t := create(t);
      wrk := Series_and_Homotopies.Shift(hom,qd_t);
      Correct(wrk,zero,3,wrk_sol,err,rco,res);
      exit when (t = 1.0);
    end loop;
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    QuadDobl_Series_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                verbose : in boolean := false ) is

    wrk : Standard_Series_Poly_Systems.Poly_Sys(hom'range);
    nit : constant integer32 := 4;
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
      put(file,"Step "); put(file,k,1);
      put_line(file," in the path tracker :");
      Series_and_Predictors.Newton_Prediction
        (file,nit,wrk,wrk_sol,srv,eva,verbose);
      new_line(file);
      put_line(file,"Setting the step size based on the power series ...");
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,tolres,verbose);
      put(file,"The computed step size : "); put(file,step,3);
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
      put(file," t = "); put(file,t,3);  new_line(file);
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,step);
      put_line(file,"Shifting the polynomial system ...");
      Standard_Series_Poly_Systems.Clear(wrk);
      wrk := Series_and_Homotopies.Shift(hom,t);
      put_line(file,"Correcting the solution ...");
      Correct(file,wrk,0.0,3,wrk_sol,err,rco,res,verbose);
      exit when (t = 1.0);
    end loop;
    put_line(file,"The solution vector :"); put_line(file,wrk_sol);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    Standard_Series_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                sol : in out DoblDobl_Complex_Solutions.Solution;
                verbose : in boolean := false ) is

    wrk : DoblDobl_Series_Poly_Systems.Poly_Sys(hom'range);
    nit : constant integer32 := 4;
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
      put(file,"Step "); put(file,k,1);
      put_line(file," in the path tracker :");
      Series_and_Predictors.Newton_Prediction
        (file,nit,wrk,wrk_sol,srv,eva,verbose);
      new_line(file);
      put_line(file,"Setting the step size based on the power series ...");
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,tolres,verbose);
      put(file,"The computed step size : "); put(file,step,3);
      if step > 0.1
       then step := 0.1;
      end if;
      update := t + step;
      if update <= onetarget then
        t := update;
      else
        step := onetarget - t;
        t := onetarget; 
      end if;
      put(file," t = "); put(file,t,3);  new_line(file);
      dd_step := create(step);
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,dd_step);
      put_line(file,"Shifting the polynomial system ...");
      DoblDobl_Series_Poly_Systems.Clear(wrk);
      dd_t := create(t);
      wrk := Series_and_Homotopies.Shift(hom,dd_t);
      put_line(file,"Correcting the solution ...");
      Correct(file,wrk,zero,3,wrk_sol,err,rco,res,verbose);
      exit when (t = 1.0);
    end loop;
    put_line(file,"The solution vector :"); put_line(file,wrk_sol);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    DoblDobl_Series_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_One_Path
              ( file : in file_type;
                hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                sol : in out Quaddobl_Complex_Solutions.Solution;
                verbose : in boolean := false ) is

    wrk : QuadDobl_Series_Poly_Systems.Poly_Sys(hom'range);
    nit : constant integer32 := 4;
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
      put(file,"Step "); put(file,k,1);
      put_line(file," in the path tracker :");
      Series_and_Predictors.Newton_Prediction
        (file,nit,wrk,wrk_sol,srv,eva,verbose);
      new_line(file);
      put_line(file,"Setting the step size based on the power series ...");
      step := Series_and_Predictors.Set_Step_Size
                (file,eva,tolcff,tolres,verbose);
      put(file,"The computed step size : "); put(file,step,3);
      if step > 0.1
       then step := 0.1;
      end if;
      update := t + step;
      if update <= onetarget then
        t := update;
      else
        step := onetarget - t;
        t := onetarget; 
      end if;
      put(file," t = "); put(file,t,3);  new_line(file);
      qd_step := create(step);
      wrk_sol := Series_and_Predictors.Predicted_Solution(srv,qd_step);
      put_line(file,"Shifting the polynomial system ...");
      QuadDobl_Series_Poly_Systems.Clear(wrk);
      qd_t := create(t);
      wrk := Series_and_Homotopies.Shift(hom,qd_t);
      put_line(file,"Correcting the solution ...");
      Correct(file,wrk,zero,3,wrk_sol,err,rco,res,verbose);
      exit when (t = 1.0);
    end loop;
    put_line(file,"The solution vector :"); put_line(file,wrk_sol);
    sol.v := wrk_sol;
    sol.err := err; sol.rco := rco; sol.res := res;
    QuadDobl_Series_Poly_Systems.Clear(wrk);
  end Track_One_Path;

  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                verbose : in boolean := false ) is

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      Series_and_Trackers.Track_One_Path(file,hom,ls.all,verbose);
      tmp := Tail_Of(tmp);
    end loop;
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in boolean := false ) is

    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      Series_and_Trackers.Track_One_Path(file,hom,ls.all,verbose);
      tmp := Tail_Of(tmp);
    end loop;
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( file : in file_type;
                hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in boolean := false ) is

    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      put(file,"Tracking path "); put(file,i,1); put_line(file," ...");
      Series_and_Trackers.Track_One_Path(file,hom,ls.all,verbose);
      tmp := Tail_Of(tmp);
    end loop;
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Series_and_Trackers.Track_One_Path(hom,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Series_and_Trackers.Track_One_Path(hom,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Track_Many_Paths;

  procedure Track_Many_Paths
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;

  begin
    for i in 1..len loop
      ls := Head_Of(tmp);
      Series_and_Trackers.Track_One_Path(hom,ls.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Track_Many_Paths;

end Series_and_Trackers;
