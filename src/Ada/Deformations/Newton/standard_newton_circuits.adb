with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Vector_Splitters;
with Standard_Complex_Linear_Solvers;

package body Standard_Newton_Circuits is

  procedure LU_Newton_Step
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
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
      put_line("The function value : "); put_line(s.fx);
      put("The residual :"); put(res,3); new_line;
    end if;
    Standard_Complex_Linear_Solvers.lufac(s.jm,s.dim,ipvt,info);
    if info /= 0 then
      if verbose
       then put("info : "); put(info,1); put_line(" singular Jacobian?");
      end if;
    else
      Standard_Complex_Vectors.Min(s.fx);
      Standard_Complex_Linear_Solvers.lusolve(s.jm,s.dim,ipvt,s.fx);
      err := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
      if verbose then
        put_line("The update : "); put_line(s.fx);
        put("Forward error :"); put(err,3); new_line;
      end if;
      for k in v'range loop
        v(k) := v(k) + s.fx(k);
      end loop;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
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
      put_line("The function value : "); put_line(s.fx);
      put("The residual :"); put(res,3); new_line;
    end if;
    Standard_Complex_Linear_Solvers.lufco(s.jm,s.dim,ipvt,rco);
    singular := (1.0 + rco = 1.0);
    if verbose then
      put("rco :"); put(rco,3);
      if singular
       then put_line(" singular Jacobian?");
       else new_line;
      end if;
    end if;
    if not singular then
      Standard_Complex_Vectors.Min(s.fx);
      Standard_Complex_Linear_Solvers.lusolve(s.jm,s.dim,ipvt,s.fx);
      err := Standard_Complex_Vector_Norms.Max_Norm(s.fx);
      if verbose then
        put_line("The update : "); put_line(s.fx);
        put("Forward error :"); put(err,3); new_line;
      end if;
      for k in v'range loop
        v(k) := v(k) + s.fx(k);
      end loop;
    end if;
  end LU_Newton_Step;

end Standard_Newton_Circuits;
