with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Vector_Splitters;
with Standard_Complex_Linear_Solvers;
with Standard_Mixed_Residuals;

package body Standard_Newton_Circuits is

  procedure LU_Newton_Step
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32; res,err : out double_float ) is
  begin
    Standard_Vector_Splitters.Complex_Parts(v,xr,xi);
    Standard_Coefficient_Circuits.EvalDiff(s,xr,xi);
    res := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
    Standard_Complex_Linear_Solvers.lufac(s.jm,s.dim,ipvt,info);
    if info = 0 then
      Standard_Complex_Vectors.Min(s.fx);
      Standard_Complex_Linear_Solvers.lusolve(s.jm,s.dim,ipvt,s.fx);
      err := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
      for k in v'range loop
        v(k) := v(k) + s.fx(k);
      end loop;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32; res,err : out double_float;
                verbose : in boolean := true ) is
  begin
    Standard_Vector_Splitters.Complex_Parts(v,xr,xi);
    Standard_Coefficient_Circuits.EvalDiff(s,xr,xi);
    res := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
    if verbose then
      put_line(file,"The approximation : "); put_line(file,v);
      put_line(file,"The function value : "); put_line(file,s.fx);
      put(file,"The residual :"); put(file,res,3); new_line(file);
    end if;
    Standard_Complex_Linear_Solvers.lufac(s.jm,s.dim,ipvt,info);
    if info /= 0 then
      if verbose then
        put(file,"info : "); put(file,info,1);
        put_line(file," singular Jacobian?");
      end if;
    else
      Standard_Complex_Vectors.Min(s.fx);
      Standard_Complex_Linear_Solvers.lusolve(s.jm,s.dim,ipvt,s.fx);
      err := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
      for k in v'range loop
        v(k) := v(k) + s.fx(k);
      end loop;
      if verbose then
        put_line(file,"The update : "); put_line(file,s.fx);
        put_line(file,"The updated approximation : "); put_line(file,v);
        put(file,"Forward error :"); put(file,err,3); new_line(file);
      end if;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                res,rco,err : out double_float ) is

    singular : boolean;

  begin
    Standard_Vector_Splitters.Complex_Parts(v,xr,xi);
    Standard_Coefficient_Circuits.EvalDiff(s,xr,xi);
    res := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
    Standard_Complex_Linear_Solvers.lufco(s.jm,s.dim,ipvt,rco);
    singular := (1.0 + rco = 1.0);
    if not singular then
      Standard_Complex_Vectors.Min(s.fx);
      Standard_Complex_Linear_Solvers.lusolve(s.jm,s.dim,ipvt,s.fx);
      err := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
      for k in v'range loop
        v(k) := v(k) + s.fx(k);
      end loop;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                res,rco,err : out double_float;
                verbose : in boolean := true ) is

    singular : boolean;

  begin
    Standard_Vector_Splitters.Complex_Parts(v,xr,xi);
    Standard_Coefficient_Circuits.EvalDiff(s,xr,xi);
    res := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
    if verbose then
      put_line(file,"The approximation : "); put_line(file,v);
      put_line(file,"The function value : "); put_line(file,s.fx);
      put(file,"The residual :"); put(file,res,3); new_line(file);
    end if;
    Standard_Complex_Linear_Solvers.lufco(s.jm,s.dim,ipvt,rco);
    singular := (1.0 + rco = 1.0);
    if verbose then
      put(file,"rco :"); put(file,rco,3);
      if singular
       then put_line(file," singular Jacobian?");
       else new_line(file);
      end if;
    end if;
    if not singular then
      Standard_Complex_Vectors.Min(s.fx);
      Standard_Complex_Linear_Solvers.lusolve(s.jm,s.dim,ipvt,s.fx);
      err := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
      for k in v'range loop
        v(k) := v(k) + s.fx(k);
      end loop;
      if verbose then
        put_line(file,"The update : "); put_line(file,s.fx);
        put_line(file,"The updated approximation : "); put_line(file,v);
        put(file,"Forward error :"); put(file,err,3); new_line(file);
      end if;
    end if;
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH MIXED RESIDUAL CALCULATION :

  procedure LU_Newton_Step
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                abscfs : in Standard_Coefficient_Circuits.Link_to_System;
                v,radv : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32; res,err,mixres : out double_float ) is
  begin
    LU_Newton_Step(s,v,xr,xi,ipvt,info,res,err);
    Standard_Vector_Splitters.Complex_Parts(v,xr,xi);
    Standard_Coefficient_Circuits.Eval(s,xr,xi);
    res := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
    radv := Standard_Mixed_Residuals.AbsVal(v);
    Standard_Vector_Splitters.Complex_Parts(radv,xr,xi);
    Standard_Coefficient_Circuits.Eval(abscfs,xr,xi);
    mixres := Standard_Mixed_Residuals.Mixed_Residual(s.fx,abscfs.fx);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                abscfs : in Standard_Coefficient_Circuits.Link_to_System;
                v,radv : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32; res,err,mixres : out double_float;
                verbose : in boolean := true ) is
  begin
    LU_Newton_Step(file,s,v,xr,xi,ipvt,info,res,err,verbose);
    Standard_Vector_Splitters.Complex_Parts(v,xr,xi);
    Standard_Coefficient_Circuits.Eval(s,xr,xi);
    res := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
    radv := Standard_Mixed_Residuals.AbsVal(v);
    Standard_Vector_Splitters.Complex_Parts(radv,xr,xi);
    Standard_Coefficient_Circuits.Eval(abscfs,xr,xi);
    mixres := Standard_Mixed_Residuals.Mixed_Residual(s.fx,abscfs.fx);
    if verbose then
      put(file,"The mixed residual :"); put(file,mixres,3); new_line(file);
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                abscfs : in Standard_Coefficient_Circuits.Link_to_System;
                v,radv : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                res,rco,err,mixres : out double_float ) is
  begin
    LU_Newton_Step(s,v,xr,xi,ipvt,res,rco,err);
    Standard_Vector_Splitters.Complex_Parts(v,xr,xi);
    Standard_Coefficient_Circuits.Eval(s,xr,xi);
    res := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
    radv := Standard_Mixed_Residuals.AbsVal(v);
    Standard_Vector_Splitters.Complex_Parts(radv,xr,xi);
    Standard_Coefficient_Circuits.Eval(abscfs,xr,xi);
    mixres := Standard_Mixed_Residuals.Mixed_Residual(s.fx,abscfs.fx);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                abscfs : in Standard_Coefficient_Circuits.Link_to_System;
                v,radv : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                res,rco,err,mixres : out double_float;
                verbose : in boolean := true ) is
  begin
    LU_Newton_Step(file,s,v,xr,xi,ipvt,res,rco,err,verbose);
    Standard_Vector_Splitters.Complex_Parts(v,xr,xi);
    Standard_Coefficient_Circuits.Eval(s,xr,xi);
    res := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
    radv := Standard_Mixed_Residuals.AbsVal(v);
    Standard_Vector_Splitters.Complex_Parts(radv,xr,xi);
    Standard_Coefficient_Circuits.Eval(abscfs,xr,xi);
    mixres := Standard_Mixed_Residuals.Mixed_Residual(s.fx,abscfs.fx);
    if verbose then
      put(file,"The mixed residual :"); put(file,mixres,3); new_line(file);
    end if;
  end LU_Newton_Step;

-- MANY NEWTON STEPS :

  procedure LU_Newton_Steps
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                maxit : in natural32; tolres,tolerr : in double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32; initres,res,err : out double_float;
                numit : out natural32; fail : out boolean;
                extra : in natural32 := 0 ) is

    cntextra : natural32 := 0;

  begin
    for k in 1..maxit+extra loop
      LU_Newton_Step(s,v,xr,xi,ipvt,info,res,err);
      if k = 1
       then initres := res;
      end if;
      if res <= tolres and err <= tolerr then -- convergence
        if (cntextra = extra) or (res = 0.0) or (err = 0.0)
         then numit := k; fail := false; return;
        end if;
        cntextra := cntextra + 1; -- do an extra step
      end if;
    end loop;
    fail := true; numit := maxit;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                maxit : in natural32; tolres,tolerr : in double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32; initres,res,err : out double_float;
                numit : out natural32; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := true ) is

    cntextra : natural32 := 0;

  begin
    for k in 1..maxit+extra loop
      LU_Newton_Step(file,s,v,xr,xi,ipvt,info,res,err,verbose);
      if k = 1
       then initres := res;
      end if;
      if res <= tolres and err <= tolerr then -- convergence
        if (cntextra = extra) or (res = 0.0) or (err = 0.0)
         then numit := k; fail := false; return;
        end if;
        cntextra := cntextra + 1; -- do an extra step
      end if;
    end loop;
    fail := true; numit := maxit;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                maxit : in natural32; tolres,tolerr : in double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                initres,res,rco,err : out double_float;
                numit : out natural32; fail : out boolean;
                extra : in natural32 := 0 ) is

    cntextra : natural32 := 0;

  begin
    for k in 1..maxit+extra loop
      LU_Newton_Step(s,v,xr,xi,ipvt,res,rco,err);
      if k = 1
       then initres := res;
      end if;
      if res <= tolres and err <= tolerr then -- convergence
        if (cntextra = extra) or (res = 0.0) or (err = 0.0)
         then numit := k; fail := false; return;
        end if;
        cntextra := cntextra + 1; -- do an extra step
      end if;
    end loop;
    fail := true; numit := maxit;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                maxit : in natural32; tolres,tolerr : in double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                initres,res,rco,err : out double_float;
                numit : out natural32; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := true ) is

    cntextra : natural32 := 0;

  begin
    for k in 1..maxit+extra loop
      LU_Newton_Step(file,s,v,xr,xi,ipvt,res,rco,err,verbose);
      if k = 1
       then initres := res;
      end if;
      if res <= tolres and err <= tolerr then -- convergence
        if (cntextra = extra) or (res = 0.0) or (err = 0.0)
         then numit := k; fail := false; return;
        end if;
        cntextra := cntextra + 1; -- do an extra step
      end if;
    end loop;
    fail := true; numit := maxit;
  end LU_Newton_Steps;

-- MANY NEWTON STEPS WITH MIXED RESIDUAL CALCULATION :

  procedure LU_Newton_Steps
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                abscfs : in Standard_Coefficient_Circuits.Link_to_System;
                v,radv : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                maxit : in natural32; tolres,tolerr : in double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32;
                initres,res,err,mixres : out double_float;
                numit : out natural32; fail : out boolean;
                extra : in natural32 := 0 ) is

    cntextra : natural32 := 0;

  begin
    for k in 1..maxit+extra loop
      LU_Newton_Step(s,abscfs,v,radv,xr,xi,ipvt,info,res,err,mixres);
      if k = 1
       then initres := res;
      end if;
      if mixres <= tolres and err <= tolerr then -- convergence
        if (cntextra = extra) or (res = 0.0) or (err = 0.0)
         then numit := k; fail := false; return;
        end if;
        cntextra := cntextra + 1; -- do an extra step
      end if;
    end loop;
    fail := true; numit := maxit;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                abscfs : in Standard_Coefficient_Circuits.Link_to_System;
                v,radv : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                maxit : in natural32; tolres,tolerr : in double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32;
                initres,res,err,mixres : out double_float;
                numit : out natural32; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := true ) is

    cntextra : natural32 := 0;

  begin
    for k in 1..maxit+extra loop
      LU_Newton_Step
        (file,s,abscfs,v,radv,xr,xi,ipvt,info,res,err,mixres,verbose);
      if k = 1
       then initres := res;
      end if;
      if mixres <= tolres and err <= tolerr then -- convergence
        if (cntextra = extra) or (res = 0.0) or (err = 0.0)
         then numit := k; fail := false; return;
        end if;
        cntextra := cntextra + 1; -- do an extra step
      end if;
    end loop;
    fail := true; numit := maxit;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                abscfs : in Standard_Coefficient_Circuits.Link_to_System;
                v,radv : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                maxit : in natural32; tolres,tolerr : in double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                initres,res,rco,err,mixres : out double_float;
                numit : out natural32; fail : out boolean;
                extra : in natural32 := 0 ) is

    cntextra : natural32 := 0;

  begin
    for k in 1..maxit+extra loop
      LU_Newton_Step(s,abscfs,v,radv,xr,xi,ipvt,res,rco,err,mixres);
      if k = 1
       then initres := res;
      end if;
      if mixres <= tolres and err <= tolerr then -- convergence
        if (cntextra = extra) or (res = 0.0) or (err = 0.0)
         then numit := k; fail := false; return;
        end if;
        cntextra := cntextra + 1; -- do an extra step
      end if;
    end loop;
    fail := true; numit := maxit;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                abscfs : in Standard_Coefficient_Circuits.Link_to_System;
                v,radv : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                maxit : in natural32; tolres,tolerr : in double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                initres,res,rco,err,mixres : out double_float;
                numit : out natural32; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := true ) is

    cntextra : natural32 := 0;

  begin
    for k in 1..maxit+extra loop
      LU_Newton_Step
        (file,s,abscfs,v,radv,xr,xi,ipvt,res,rco,err,mixres,verbose);
      if k = 1
       then initres := res;
      end if;
      if mixres <= tolres and err <= tolerr then -- convergence
        if (cntextra = extra) or (res = 0.0) or (err = 0.0)
         then numit := k; fail := false; return;
        end if;
        cntextra := cntextra + 1; -- do an extra step
      end if;
    end loop;
    fail := true; numit := maxit;
  end LU_Newton_Steps;

end Standard_Newton_Circuits;
