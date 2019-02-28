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
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;

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

  procedure Correct
              ( nq : in integer32; t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean ) is

    prev_err,prev_res,nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    prev_err100,prev_res100 : double_float;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit loop
      Standard_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
      nrm := Standard_Complex_Vector_Norms.Max_Norm(sol);
     -- if res <= tolres then -- convergence
      left := err/(nrm+1.0);
      if left <= tolres then
        nbrit := k; fail := false; exit;
      elsif k > 1 then      -- check for divergence
        prev_res100 := 100.0*prev_res;
        prev_err100 := 100.0*prev_err;
        if ((res > prev_res100) and (err > prev_err100))
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
                err,rco,res : out double_float; fail : out boolean ) is

    prev_err,prev_res,nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    solnrm : double_double;
    prev_err100,prev_res100 : double_float;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit loop
      DoblDobl_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
      solnrm := DoblDobl_Complex_Vector_Norms.Max_Norm(sol);
      nrm := hi_part(solnrm);
     -- if res <= tolres then -- convergence
      left := err/(nrm+1.0);
      if left <= tolres then
        nbrit := k; fail := false; exit;
      elsif k > 1 then      -- check for divergence
        prev_res100 := 100.0*prev_res;
        prev_err100 := 100.0*prev_err;
        if ((res > prev_res100) and (err > prev_err100))
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
                err,rco,res : out double_float; fail : out boolean ) is

    prev_err,prev_res,nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    solnrm : quad_double;
    prev_err100,prev_res100 : double_float;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit loop
      QuadDobl_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
      solnrm := QuadDobl_Complex_Vector_Norms.Max_Norm(sol);
      nrm := hihi_part(solnrm);
     -- if res <= tolres then -- convergence
      left := err/(nrm+1.0);
      if left <= tolres then
        nbrit := k; fail := false; exit;
      elsif k > 1 then      -- check for divergence
        prev_res100 := 100.0*prev_res;
        prev_err100 := 100.0*prev_err;
        if ((res > prev_res100) and (err > prev_err100))
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
                verbose : in boolean := false ) is

    prev_err,prev_res,nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    prev_err100,prev_res100 : double_float;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit loop
      Standard_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
      nrm := Standard_Complex_Vector_Norms.Max_Norm(sol);
      if verbose then
        put(file,"  nrm :"); put(file,nrm,3);
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
     -- if res <= tolres then -- convergence
      left := err/(nrm+1.0);
      if left <= tolres then
        nbrit := k; fail := false; exit;
      elsif k > 1 then      -- check for divergence
        prev_res100 := 100.0*prev_res;
        prev_err100 := 100.0*prev_err;
        if ((res > prev_res100) and (err > prev_err100))
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
                verbose : in boolean := false ) is

    prev_err,prev_res,nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    solnrm : double_double;
    prev_err100,prev_res100 : double_float;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit loop
      DoblDobl_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
      solnrm := DoblDobl_Complex_Vector_Norms.Max_Norm(sol);
      nrm := hi_part(solnrm);
      if verbose then
        put(file,"  nrm :"); put(file,nrm,3);
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
     -- if res <= tolres then -- convergence
      left := err/(nrm+1.0);
      if left <= tolres then
        nbrit := k; fail := false; exit;
      elsif k > 1 then      -- check for divergence
        prev_res100 := 100.0*prev_res;
        prev_err100 := 100.0*prev_err;
        if ((res > prev_res100) and (err > prev_err100))
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
                verbose : in boolean := false ) is

    prev_err,prev_res,nrm,left : double_float := 1.0;
    cmplxt : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(t);
    solnrm : quad_double;
    prev_err100,prev_res100 : double_float;

  begin
    fail := true;
    nbrit := maxit;
    for k in 1..maxit loop
      QuadDobl_LU_Newton_Step(nq,cmplxt,sol,err,rco,res);
      solnrm := QuadDobl_Complex_Vector_Norms.Max_Norm(sol);
      nrm := hihi_part(solnrm);
      if verbose then
        put(file,"  nrm :"); put(file,nrm,3);
        put(file,"  err :"); put(file,err,3);
        put(file,"  rco :"); put(file,rco,3);
        put(file,"  res :"); put(file,res,3); new_line(file);
      end if;
     -- if res <= tolres then -- convergence
      left := err/(nrm+1.0);
      if left <= tolres then
        nbrit := k; fail := false; exit;
      elsif k > 1 then      -- check for divergence
        prev_err100 := 100.0*prev_err;
        prev_res100 := 100.0*prev_res;
        if ((res > prev_res100) and (err > prev_err100))
         then nbrit := k; exit;
        end if;
      end if;
      prev_err := err; -- previous forward error
      prev_res := res; -- previous backward error
    end loop;
  end Correct;

end Homotopy_Newton_Steps;
