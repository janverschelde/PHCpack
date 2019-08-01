with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;
with Standard_Integer_Vectors;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_Singular_Values;   use Standard_Complex_Singular_Values;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Standard_Homotopy;
with Standard_Coefficient_Homotopy;
with DoblDobl_Homotopy;
with DoblDobl_Coefficient_Homotopy;
with QuadDobl_Homotopy;
with QuadDobl_Coefficient_Homotopy;
with Homotopy_Mixed_Residuals;

package body Homotopy_Newton_Steps is

  procedure Standard_LU_Newton_Step
              ( nq : in integer32;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    y : Standard_Complex_Vectors.Vector(1..nq)
      := Standard_Homotopy.Eval(x,t);
    A : Standard_Complex_Matrices.Matrix(1..nq,x'range)
      := Standard_Homotopy.Diff(x,t);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant double_float := Norm1(A);

  begin
    Standard_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    Standard_Complex_Vectors.Add(x,y);
    err := Standard_Complex_Vector_Norms.Max_Norm(y);
    y := Standard_Homotopy.Eval(x,t);
    res := Standard_Complex_Vector_Norms.Max_Norm(y);
  end Standard_LU_Newton_Step;

  procedure DoblDobl_LU_Newton_Step
              ( nq : in integer32;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    use DoblDobl_Complex_Numbers_cv;

    dd_t : constant DoblDobl_Complex_Numbers.Complex_Number
         := Standard_to_DoblDobl_Complex(t);
    y : DoblDobl_Complex_Vectors.Vector(1..nq)
      := DoblDobl_Homotopy.Eval(x,dd_t);
    A : DoblDobl_Complex_Matrices.Matrix(1..nq,x'range)
      := DoblDobl_Homotopy.Diff(x,dd_t);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant double_double := Norm1(A);
    dd_rco : double_double;

  begin
    DoblDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,dd_rco);
    rco := hi_part(dd_rco);
    lusolve(A,A'last(1),ipvt,y);
    DoblDobl_Complex_Vectors.Add(x,y);
    err := hi_part(DoblDobl_Complex_Vector_Norms.Max_Norm(y));
    y := DoblDobl_Homotopy.Eval(x,dd_t);
    res := hi_part(DoblDobl_Complex_Vector_Norms.Max_Norm(y));
  end DoblDobl_LU_Newton_Step;

  procedure QuadDobl_LU_Newton_Step
              ( nq : in integer32;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    use QuadDobl_Complex_Numbers_cv;

    qd_t : constant QuadDobl_Complex_Numbers.Complex_Number
         := Standard_to_QuadDobl_Complex(t);
    y : QuadDobl_Complex_Vectors.Vector(1..nq)
      := QuadDobl_Homotopy.Eval(x,qd_t);
    A : QuadDobl_Complex_Matrices.Matrix(1..nq,x'range)
      := QuadDobl_Homotopy.Diff(x,qd_t);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant quad_double := Norm1(A);
    qd_rco : quad_double;

  begin
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,qd_rco);
    rco := hihi_part(qd_rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := hihi_part(QuadDobl_Complex_Vector_Norms.Max_Norm(y));
    y := QuadDobl_Homotopy.Eval(x,qd_t);
    res := hihi_part(QuadDobl_Complex_Vector_Norms.Max_Norm(y));
  end QuadDobl_LU_Newton_Step;

  procedure Standard_LU_Newton_Step
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    nq : constant integer32 := abh'last;
    y : Standard_Complex_Vectors.Vector(1..nq);
    A : Standard_Complex_Matrices.Matrix(1..nq,x'range);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : double_float;

  begin
    if Standard_Coefficient_Homotopy.Number_of_Equations = -1 then
      y := Standard_Homotopy.Eval(x,t);
      A := Standard_Homotopy.Diff(x,t);
    else
      y := Standard_Coefficient_Homotopy.Eval(x,t);
      A := Standard_Coefficient_Homotopy.Diff(x,t);
    end if;
    Anorm := Norm1(A);
    Standard_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    Standard_Complex_Vectors.Add(x,y);
    err := Standard_Complex_Vector_Norms.Max_Norm(y);
    res := Homotopy_mixed_Residuals.Residual(abh,x,t);
  end Standard_LU_Newton_Step;

  procedure Standard_SVD_Newton_Step
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    nq : constant integer32 := abh'last;
    y : Standard_Complex_Vectors.Vector(1..nq);
    A : Standard_Complex_Matrices.Matrix(1..nq,x'range);
    u : Standard_Complex_Matrices.Matrix(y'range,y'range);
    v : Standard_Complex_Matrices.Matrix(x'range,x'range);
    p : constant integer32 := x'length;
    m : constant integer32 := Min0(nq+1,p);
    e : Standard_Complex_Vectors.Vector(1..p);
    s : Standard_Complex_Vectors.Vector(1..m);
    info : integer32;
    dx : Standard_Complex_Vectors.Vector(x'range);

  begin
    if Standard_Coefficient_Homotopy.Number_of_Equations = -1 then
      y := Standard_Homotopy.Eval(x,t);
      A := Standard_Homotopy.Diff(x,t);
    else
      y := Standard_Coefficient_Homotopy.Eval(x,t);
      A := Standard_Coefficient_Homotopy.Diff(x,t);
    end if;
    SVD(A,nq,p,s,e,u,v,11,info);
    rco := Inverse_Condition_Number(s);
    Standard_Complex_Vectors.Min(y);
    dx := Solve(u,v,s,y);
    err := Standard_Complex_Vector_Norms.Max_Norm(dx);
    Standard_Complex_Vectors.Add(x,dx);
    res := Homotopy_mixed_Residuals.Residual(abh,x,t);
  end Standard_SVD_Newton_Step;

  procedure DoblDobl_LU_Newton_Step
              ( abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    use DoblDobl_Complex_Numbers_cv;

    nq : constant integer32 := abh'last;
    dd_t : constant DoblDobl_Complex_Numbers.Complex_Number
         := Standard_to_DoblDobl_Complex(t);
    y : DoblDobl_Complex_Vectors.Vector(1..nq);
    A : DoblDobl_Complex_Matrices.Matrix(1..nq,x'range);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm,dd_rco,dd_res : double_double;

  begin
    if DoblDobl_Coefficient_Homotopy.Number_of_Equations = -1 then
      y := DoblDobl_Homotopy.Eval(x,dd_t);
      A := DoblDobl_Homotopy.Diff(x,dd_t);
    else
      y := DoblDobl_Coefficient_Homotopy.Eval(x,dd_t);
      A := DoblDobl_Coefficient_Homotopy.Diff(x,dd_t);
    end if;
    Anorm := Norm1(A);
    DoblDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,dd_rco);
    rco := hi_part(dd_rco);
    lusolve(A,A'last(1),ipvt,y);
    DoblDobl_Complex_Vectors.Add(x,y);
    err := hi_part(DoblDobl_Complex_Vector_Norms.Max_Norm(y));
    dd_res := Homotopy_mixed_Residuals.Residual(abh,x,dd_t);
    res := hi_part(dd_res);
  end DoblDobl_LU_Newton_Step;

  procedure QuadDobl_LU_Newton_Step
              ( abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    use QuadDobl_Complex_Numbers_cv;

    nq : constant integer32 := abh'last;
    qd_t : constant QuadDobl_Complex_Numbers.Complex_Number
         := Standard_to_QuadDobl_Complex(t);
    y : QuadDobl_Complex_Vectors.Vector(1..nq);
    A : QuadDobl_Complex_Matrices.Matrix(1..nq,x'range);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm,qd_rco,qd_res : quad_double;

  begin
    if QuadDobl_Coefficient_Homotopy.Number_of_Equations = -1 then
      y := QuadDobl_Homotopy.Eval(x,qd_t);
      A := QuadDobl_Homotopy.Diff(x,qd_t);
    else
      y := QuadDobl_Coefficient_Homotopy.Eval(x,qd_t);
      A := QuadDobl_Coefficient_Homotopy.Diff(x,qd_t);
    end if;
    Anorm := Norm1(A);
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,qd_rco);
    rco := hihi_part(qd_rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := hihi_part(QuadDobl_Complex_Vector_Norms.Max_Norm(y));
    qd_res := Homotopy_mixed_Residuals.Residual(abh,x,qd_t);
    res := hihi_part(qd_res);
  end QuadDobl_LU_Newton_Step;

  procedure Correct
              ( nq : in integer32; t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 ) is

    prev_err,prev_res : double_float := 1.0;
   -- nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      Standard_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
     -- nrm := Standard_Complex_Vector_Norms.Max_Norm(sol);
     -- if res <= tolres then -- convergence
     -- left := err/(nrm+1.0);
      if err <= tolres and res <= tolres then -- convergence
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if ((res > prev_res10) and (err > prev_err10))
         then nbrit := k; exit;
        end if;
      end if;
      prev_err := err; -- previous forward error
      prev_res := res; -- previous backward error
    end loop;
  end Correct;

  procedure Correct
              ( nq : in integer32; t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 ) is

    prev_err,prev_res : double_float := 1.0;
   -- nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
   -- solnrm : double_double;
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      DoblDobl_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
     -- solnrm := DoblDobl_Complex_Vector_Norms.Max_Norm(sol);
     -- nrm := hi_part(solnrm);
     -- if res <= tolres then -- convergence
     -- left := err/(nrm+1.0);
      if err <= tolres and res <= tolres then
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if ((res > prev_res10) and (err > prev_err10))
         then nbrit := k; exit;
        end if;
      end if;
      prev_err := err; -- previous forward error
      prev_res := res; -- previous backward error
    end loop;
  end Correct;

  procedure Correct
              ( nq : in integer32; t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 ) is

    prev_err,prev_res : double_float := 1.0;
   -- nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
   -- solnrm : quad_double;
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      QuadDobl_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
     -- solnrm := QuadDobl_Complex_Vector_Norms.Max_Norm(sol);
     -- nrm := hihi_part(solnrm);
     -- if res <= tolres then -- convergence
     -- left := err/(nrm+1.0);
      if err <= tolres and res <= tolres then
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if ((res > prev_res10) and (err > prev_err10))
         then nbrit := k; exit;
        end if;
      end if;
      prev_err := err; -- previous forward error
      prev_res := res; -- previous backward error
    end loop;
  end Correct;

  procedure Correct
              ( file : in file_type;
                nq : in integer32; t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false ) is

    prev_err,prev_res : double_float := 1.0;
   -- nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      Standard_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
     -- nrm := Standard_Complex_Vector_Norms.Max_Norm(sol);
      if verbose then
       -- put(file,"  nrm :"); put(file,nrm,3);
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
     -- if res <= tolres then -- convergence
     -- left := err/(nrm+1.0);
      if err <= tolres and res <= tolres then
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if ((res > prev_res10) and (err > prev_err10))
         then nbrit := k; exit;
        end if;
      end if;
      prev_err := err; -- previous forward error
      prev_res := res; -- previous backward error
    end loop;
  end Correct;

  procedure Correct
              ( file : in file_type;
                nq : in integer32; t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false ) is

    prev_err,prev_res : double_float := 1.0;
   -- nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
   -- solnrm : double_double;
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      DoblDobl_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
     -- solnrm := DoblDobl_Complex_Vector_Norms.Max_Norm(sol);
     -- nrm := hi_part(solnrm);
      if verbose then
       -- put(file,"  nrm :"); put(file,nrm,3);
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
     -- if res <= tolres then -- convergence
     -- left := err/(nrm+1.0);
      if err <= tolres and err <= tolres then
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if ((res > prev_res10) and (err > prev_err10))
         then nbrit := k; exit;
        end if;
      end if;
      prev_err := err; -- previous forward error
      prev_res := res; -- previous backward error
    end loop;
  end Correct;

  procedure Correct
              ( file : in file_type;
                nq : in integer32; t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false ) is

    prev_err,prev_res : double_float := 1.0;
   -- nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
   -- solnrm : quad_double;
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      QuadDobl_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
     -- solnrm := QuadDobl_Complex_Vector_Norms.Max_Norm(sol);
     -- nrm := hihi_part(solnrm);
      if verbose then
       -- put(file,"  nrm :"); put(file,nrm,3);
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
     -- if res <= tolres then -- convergence
     -- left := err/(nrm+1.0);
      if err <= tolres and res <= tolres then
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_err10 := 10.0*prev_err;
        prev_res10 := 10.0*prev_res;
        if ((res > prev_res10) and (err > prev_err10))
         then nbrit := k; exit;
        end if;
      end if;
      prev_err := err; -- previous forward error
      prev_res := res; -- previous backward error
    end loop;
  end Correct;

-- CORRECTORS WITH MIXED RESIDUALS :

  procedure Correct
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 ) is

    prev_err,prev_res : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    prev_res10,prev_err10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      Standard_LU_Newton_Step(abh,cmplxt,sol,err,rco,res);
     -- Standard_SVD_Newton_Step(abh,cmplxt,sol,err,rco,res);
      if err <= tolres and res <= tolres then -- convergence
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if (res > prev_res10) and (err > prev_err10)
         then nbrit := k; exit;
        end if;
      end if;
      prev_res := res; -- previous residual
      prev_err := err;
    end loop;
  end Correct;

  procedure Correct
              ( abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 ) is

    prev_err,prev_res : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      DoblDobl_LU_Newton_Step(abh,cmplxt,sol,err,rco,res);
      if err <= tolres and res <= tolres then -- convergence
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if (res > prev_res10) and (err > prev_err10)
         then nbrit := k; exit;
        end if;
      end if;
      prev_res := res; -- previous residual
      prev_err := err;
    end loop;
  end Correct;

  procedure Correct
              ( abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 ) is

    prev_err,prev_res : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      QuadDobl_LU_Newton_Step(abh,cmplxt,sol,err,rco,res);
      if err <= tolres and res <= tolres then -- convergence
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if (res > prev_res10) and (err > prev_err10)
         then nbrit := k; exit;
        end if;
      end if;
      prev_res := res; -- previous residual
      prev_err := err;
    end loop;
  end Correct;

  procedure Correct
              ( file : in file_type;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false ) is

    prev_err,prev_res : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      Standard_LU_Newton_Step(abh,cmplxt,sol,err,rco,res);
     -- Standard_SVD_Newton_Step(abh,cmplxt,sol,err,rco,res);
      if verbose then
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
      if err <= tolres and res <= tolres then -- convergence
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if (res > prev_res10) and (err > prev_err10)
         then nbrit := k; exit;
        end if;
      end if;
      prev_res := res; -- previous residual
      prev_err := err;
    end loop;
  end Correct;

  procedure Correct
              ( file : in file_type;
                abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false ) is

    prev_err,prev_res : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      DoblDobl_LU_Newton_Step(abh,cmplxt,sol,err,rco,res);
      if verbose then
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
      if err <= tolres and res <= tolres then -- convergence
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if (res > prev_res10) and (err > prev_err10)
         then nbrit := k; exit;
        end if;
      end if;
      prev_res := res; -- previous residual
      prev_err := err;
    end loop;
  end Correct;

  procedure Correct
              ( file : in file_type;
                abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false ) is

    prev_err,prev_res : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    prev_err10,prev_res10 : double_float;
    cntextra : natural32 := 0;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit+extra loop
      QuadDobl_LU_Newton_Step(abh,cmplxt,sol,err,rco,res);
      if verbose then
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
      if err <= tolres and res <= tolres then -- convergence
        if extra = 0 or (err = 0.0) or (res = 0.0) or (cntextra = extra)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1;
      elsif k > 1 then      -- check for divergence
        prev_res10 := 10.0*prev_res;
        prev_err10 := 10.0*prev_err;
        if (res > prev_res10) and (err > prev_err10)
         then nbrit := k; exit;
        end if;
      end if;
      prev_res := res; -- previous residual
      prev_err := err;
    end loop;
  end Correct;

end Homotopy_Newton_Steps;
