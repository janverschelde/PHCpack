with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers_Polar;    use QuadDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Singular_Values;  use QuadDobl_Complex_Singular_Values;
with Process_io;
with QuadDobl_Plane_Representations;
with QuadDobl_Moving_Planes;
with QuadDobl_Intrinsic_Newton;         use QuadDobl_Intrinsic_Newton;
with QuadDobl_Rescaling_Coordinates;    use QuadDobl_Rescaling_Coordinates;

package body QuadDobl_Intrinsic_Trackers is

-- AUXILIARIES :

  function At_End ( s : Solu_Info; t1 : Complex_Number;
                    max : natural32 ) return boolean is

  -- DESCRIPTION :
  --   We are at the end of the path if s.sol.t equals t1,
  --   or if s.niter equals or exceeds max.

  begin
    return (Equal(s.sol.t,t1) or (s.niter > max));
  end At_End;

  procedure Step_Control
              ( fail : in boolean; pp : in Pred_Pars;
                step : in out quad_double; nbsuccess : in out natural32 ) is

  -- DESCRIPTION :
  --   Applies step control and management of data for step back.

  -- ON ENTRY :
  --   fail     if true, corrector step failed, otherwise if false;
  --   pp       settings for predictor values;
  --   step     current step size;
  --   nbsucces is the number of consecutive corrector successes.

  -- ON RETURN :
  --   step     updated step size;
  --   nbsucces is updated accordingly to fail.

  begin
    if fail then
      nbsuccess := 0;
      step := step*pp.redfac;
    else
      nbsuccess := nbsuccess + 1;
      if nbsuccess > pp.success_steps then
        step := step*pp.expfac;
        if step > pp.maxstep
         then step := create(pp.maxstep);
        end if;
      end if;
    end if;
  end Step_Control;

  procedure Step_Control
              ( file : in file_type; fail : in boolean; pp : in Pred_Pars;
                step : in out quad_double; nbsuccess : in out natural32 ) is

  -- DESCRIPTION :
  --   Applies step control and management of data for step back.

  -- ON ENTRY :
  --   fail     if true, corrector step failed, otherwise if false;
  --   pp       settings for predictor values;
  --   step     current step size;
  --   nbsucces is the number of consecutive corrector successes.

  -- ON RETURN :
  --   step     updated step size;
  --   nbsucces is updated accordingly to fail.

  begin
    if fail then
      nbsuccess := 0;
      step := step*pp.redfac;
    else
      nbsuccess := nbsuccess + 1;
      if nbsuccess > pp.success_steps then
        step := step*pp.expfac;
        if step > pp.maxstep
         then step := create(pp.maxstep);
        end if;
      end if;
    end if;
    if Process_io.Contains_Output_Code(Process_io.p) then
      put(file,"#consecutive successes : ");
      put(file,nbsuccess,1); new_line(file);
    end if;
  end Step_Control;

  procedure Step_Control
              ( fail : in boolean; pp : in Pred_Pars;
                step : in out quad_double; nbsuccess : in out natural32;
                t,pt0,pt1,pt2,pt3 : in out Complex_Number;
                x,px0,px1,px2,px3 : in out Vector ) is

  -- DESCRIPTION :
  --   Applies step control and management of data for step back.

  -- ON ENTRY :
  --   fail     if true, corrector step failed, otherwise if false;
  --   pp       settings for predictor values;
  --   step     current step size;
  --   nbsucces is the number of consecutive corrector successes;
  --   t        current value of continuation parameter;
  --   pt0      previous value of t where corrector had success;
  --   pt1      value previous to pt0 where corrector had success;
  --   pt2      value previous to pt1 where corrector had success;
  --   pt3      value previous to pt2 where corrector had success;
  --   x        current solution vector corresponding to t;
  --   px0      current solution vector corresponding to pt0;
  --   px1      current solution vector corresponding to pt1;
  --   px2      current solution vector corresponding to pt2;
  --   px3      current solution vector corresponding to pt3.

  -- ON RETURN :
  --   step     updated step size;
  --   nbsucces is updated accordingly to fail;
  --   t        equals pt0 when fail, otherwise is unchanged;
  --   pt0      equals t when not fail, otherwise is unchanged;
  --   pt1      equals pt0 when not fail, otherwise is unchanged;
  --   pt2      equals pt1 when not fail, otherwise is unchanged;
  --   pt3      equals pt2 when not fail, otherwise is unchanged;
  --   px0      equals x when not fail, otherwise is unchanged;
  --   px1      equals px0 when not fail, otherwise is unchanged;
  --   px2      equals px1 when not fail, otherwise is unchanged;
  --   px3      equals px2 when not fail, otherwise is unchanged.

  begin
    if fail then
      nbsuccess := 0;
      step := step*pp.redfac;
    else
      nbsuccess := nbsuccess + 1;
      if nbsuccess > pp.success_steps then
        step := step*pp.expfac;
        if step > pp.maxstep
         then step := create(pp.maxstep);
        end if;
      end if;
    end if;
    if step <= pp.minstep
     then return;
    end if;
    if fail then
      t := pt0; x := px0;
     --t := pt1; x := px1;      -- 01/13 why this ?
     --pt0 := pt1; px0 := px1;
    else 
      pt3 := pt2; px3 := px2;
      pt2 := pt1; px2 := px1;
      pt1 := pt0; px1 := px0;
      pt0 := t;   px0 := x;
    end if;
  end Step_Control;

  procedure Step_Control
              ( file : in file_type; fail : in boolean; pp : in Pred_Pars;
                step : in out quad_double; nbsuccess : in out natural32;
                t,pt0,pt1,pt2,pt3 : in out Complex_Number;
                x,px0,px1,px2,px3 : in out Vector ) is

  -- DESCRIPTION :
  --   Applies step control and management of data for step back.

  -- ON ENTRY :
  --   fail     if true, corrector step failed, otherwise if false;
  --   pp       settings for predictor values;
  --   step     current step size;
  --   nbsucces is the number of consecutive corrector successes;
  --   t        current value of continuation parameter;
  --   pt0      previous value of t where corrector had success;
  --   pt1      value previous to pt0 where corrector had success;
  --   pt2      value previous to pt1 where corrector had success;
  --   pt3      value previous to pt2 where corrector had success;
  --   x        current solution vector corresponding to t;
  --   px0      current solution vector corresponding to pt0;
  --   px1      current solution vector corresponding to pt1;
  --   px2      current solution vector corresponding to pt2;
  --   px3      current solution vector corresponding to pt3.

  -- ON RETURN :
  --   step     updated step size;
  --   nbsucces is updated accordingly to fail;
  --   t        equals pt0 when fail, otherwise is unchanged;
  --   pt0      equals t when not fail, otherwise is unchanged;
  --   pt1      equals pt0 when not fail, otherwise is unchanged;
  --   pt2      equals pt1 when not fail, otherwise is unchanged;
  --   pt3      equals pt2 when not fail, otherwise is unchanged;
  --   px0      equals x when not fail, otherwise is unchanged;
  --   px1      equals px0 when not fail, otherwise is unchanged;
  --   px2      equals px1 when not fail, otherwise is unchanged;
  --   px3      equals px2 when not fail, otherwise is unchanged.

  begin
    if fail then
      nbsuccess := 0;
      step := step*pp.redfac;
    else
      nbsuccess := nbsuccess + 1;
      if nbsuccess > pp.success_steps then
        step := step*pp.expfac;
        if step > pp.maxstep
         then step := create(pp.maxstep);
        end if;
      end if;
    end if;
    if Process_io.Contains_Output_Code(Process_io.p) then
      put(file,"#consecutive successes : ");
      put(file,nbsuccess,1); new_line(file);
    end if;
    if step <= pp.minstep
     then return;
    end if;
    if fail then
      t := pt0; x := px0;
     -- put_line(file,"Back up solution to vector :");
     -- put_line(file,x);
     --t := pt1; x := px1;      -- 01/13 why this ?
     --pt0 := pt1; px0 := px1;
    else
      pt3 := pt2; px3 := px2;
      pt2 := pt1; px2 := px1;
      pt1 := pt0; px1 := px0;
      pt0 := t;   px0 := x;
    end if;
   -- put_line(file,"Current vector of t-values : ");
   -- put(file,"t3 = "); put(file,pt3); new_line(file);
   -- put(file,"t2 = "); put(file,pt2); new_line(file);
   -- put(file,"t1 = "); put(file,pt1); new_line(file);
   -- put(file,"t0 = "); put(file,pt0); new_line(file);
  end Step_Control;

  function Quadratic_Predictor
               ( t,t0,t1,t2 : Complex_Number; x0,x1,x2 : Vector )
               return Vector is

  -- DESCRIPTION :
  --   Uses extrapolation through (t0,x0), (t1,x1), and (t2,x2) to
  --   predict the solution vector at t.

    v : Vector(x0'range);
    f01,f02,f12,dt01,dt02,dt12,pt0,pt01 : Complex_Number;

  begin
    dt01 := t1 - t0;
    dt02 := t2 - t0;
    dt12 := t2 - t1;
    pt0 := t - t0;
    pt01 := pt0*(t-t1);
    for i in v'range loop
      f01 := (x1(i) - x0(i))/dt01;
      f02 := (x2(i) - x0(i))/dt02;
      f12 := (f02 - f01)/dt12;
      v(i) := x0(i) + f01*pt0 + f12*pt01;
    end loop;
    return v;
  end Quadratic_Predictor;

  function Cubic_Predictor
               ( t,t0,t1,t2,t3 : Complex_Number; x0,x1,x2,x3 : Vector )
               return Vector is

  -- DESCRIPTION :
  --   Uses extrapolation through (t0,x0), (t1,x1), (t2,x2), and (t3,x3)
  --   to predict the solution vector at t.

    v : Vector(x0'range);
    f1,f2,f3,pt0,pt1,pt2,dt01,dt02,dt03,dt12,dt13,dt23 : Complex_Number;

  begin
    pt0 := t-t0;
    pt1 := pt0*(t-t1);
    pt2 := pt1*(t-t2);
    dt01 := t1-t0; dt02 := t2-t0; dt03 := t3-t0;
    dt12 := t2-t1; dt13 := t3-t1; dt23 := t3-t2;
    for i in v'range loop
      f1 := (x1(i) - x0(i))/dt01;
      f2 := (x2(i) - x0(i))/dt02;
      f3 := (x3(i) - x0(i))/dt03;
      f2 := (f2 - f1)/dt12;
      f3 := (f3 - f1)/dt13;
      f3 := (f3 - f2)/dt23;
      v(i) := x0(i) + f1*pt0 + f2*pt1 + f3*pt2;
    end loop;
    return v;
  end Cubic_Predictor;

  procedure Predictor
               ( t : in out Complex_Number; v : in out Vector;
                 step : in out quad_double;
                 pt0,pt1,pt2,pt3,t1 : in Complex_Number;
                 px1,px2,px3 : in Vector; fail : in boolean ) is

  -- DESCRIPTION :
  --   Predicts new value for continuation parameter t and,
  --   if not fail, also predicts the corresponding solution vector v.  

  -- ON ENTRY :
  --   t         current value of the continuation parameter t;
  --   v         current value of the solution vector v;
  --   step      current step size;
  --   pt0       backup value for the continuation parameter t;
  --   pt1       value for t previous to backup for t;
  --   pt2       value for t previous to pt1;
  --   pt3       value for t previous to pt2;
  --   t1        target value for continuation parameter (t1 = 1);
  --   px1       previous value for solution vector v;
  --   px2       value for solution vector v corresponding to pt2;
  --   px3       value for solution vector v corresponding to pt3;
  --   fail      true if last correction was not successful.

  -- ON RETURN :
  --   t         predicted value for the continuation parameter t;
  --   v         predicted value for the solution vector v;
  --   step      is actual step size taken.

    d0,d1,d2 : quad_double;
    one : constant quad_double := create(1.0);

  begin
    t := t + Create(step);
    if REAL_PART(t) > 1.0 then
      step := to_double(REAL_PART(t)) - 1.0;
      t := t1;
    end if;
    if not fail then
      d0 := AbsVal(pt0 - pt1);
      if one + d0 /= one then
        d1 := AbsVal(pt1 - pt2);
        if one + d1 = one then
          d0 := step/d0;
          v := v + Create(d0)*(v - px1);
        else
          d2 := AbsVal(pt2 - pt3);
          if one + d2 = one 
           then v := Quadratic_Predictor(t,pt0,pt1,pt2,v,px1,px2);
           else v := Cubic_Predictor(t,pt0,pt1,pt2,pt3,v,px1,px2,px3);
          end if;
        end if;
      end if;
    end if;
  end Predictor;

  procedure Predictor
               ( file : in file_type;
                 t : in out Complex_Number; v : in out Vector;
                 step : in out quad_double;
                 pt0,pt1,pt2,pt3,t1 : in Complex_Number;
                 px1,px2,px3 : in Vector; fail : in boolean ) is

  -- DESCRIPTION :
  --   Predicts new value for continuation parameter t and,
  --   if not fail, also predicts the corresponding solution vector v.  

  -- ON ENTRY :
  --   file      for intermediate diagnostics;
  --   t         current value of the continuation parameter t;
  --   v         current value of the solution vector v;
  --   step      current step size;
  --   pt0       backup value for the continuation parameter t;
  --   pt1       value for t previous to backup for t;
  --   pt2       value for t previous to pt1;
  --   pt3       value for t previous to pt2;
  --   t1        target value for continuation parameter (t1 = 1);
  --   px1       previous value for solution vector v;
  --   px2       value for solution vector v corresponding to pt2;
  --   px3       value for solution vector v corresponding to pt3;
  --   fail      true if last correction was not successful.

  -- ON RETURN :
  --   t         predicted value for the continuation parameter t;
  --   v         predicted value for the solution vector v;
  --   step      is actual step size taken.

    d0,d1,d2 : quad_double;
    one : constant quad_double := create(1.0);

  begin
    t := t + Create(step);
    if REAL_PART(t) > 1.0 then
      step := to_double(REAL_PART(t)) - 1.0;
      t := t1;
    end if;
    if fail then
      put_line(file,"No predictor for solution.");
    else
      d0 := AbsVal(pt0 - pt1);
      if one + d0 = one then
        put_line(file,"No predictor for solution.");
      else
        d1 := AbsVal(pt1 - pt2);
        if one + d1 = one then
          put_line(file,"Secant predictor for solution.");
          d0 := step/d0;
          v := v + Create(d0)*(v - px1);
        else
          d2 := AbsVal(pt2 - pt3);
          if one + d2 = one
           then put_line(file,"Quadratic predictor for solution.");
                v := Quadratic_Predictor(t,pt0,pt1,pt2,v,px1,px2);
           else put_line(file,"Cubic predictor for solution.");
                v := Cubic_Predictor(t,pt0,pt1,pt2,pt3,v,px1,px2,px3);
          end if;
        end if;
      end if;
    end if;
  end Predictor;

-- TARGET ROUTINES :

  procedure Silent_Affine_LU_Track 
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    pt0,pt1,pt2,pt3 : Complex_Number := s.sol.t;
    px0,px1,px2,px3 : Vector(s.sol.v'range) := s.sol.v;
    step : quad_double := create(pp.maxstep);
    plane : Matrix(jf'range(2),0..s.sol.v'last);
    fail : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    while not At_End(s,t1,cp.maxtot) or fail loop
      Predictor(s.sol.t,s.sol.v,step,pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
      plane := Path(s.sol.t);
      Affine_LU_Newton
        (f,jf,plane,s.sol.v,
         create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
         scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.niter := s.niter + nbit;
      s.nstep := s.nstep + 1;
      if fail then s.nfail := s.nfail + 1; end if;
      Step_Control(fail,pp,step,nbsuccess,s.sol.t,
                   pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      exit when (step <= pp.minstep);
    end loop;
    if (step > pp.minstep) then
      plane := Path(Create(one));
      Affine_LU_Newton
        (f,jf,plane,s.sol.v,
         create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
         scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,fail);
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.rcond := to_double(srcond);
    end if;
  end Silent_Affine_LU_Track;

  procedure Silent_Projective_LU_Track 
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info; k : in natural32;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    pt0,pt1,pt2,pt3 : Complex_Number := s.sol.t;
    px0,px1,px2,px3 : Vector(s.sol.v'range) := s.sol.v;
    step : quad_double := create(pp.maxstep);
    plane : Matrix(jf'range(2),0..s.sol.v'last);
    fail : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    while not At_End(s,t1,cp.maxtot) or fail loop
      Predictor(s.sol.t,s.sol.v,step,pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
      plane := Path(s.sol.t);
      Projective_LU_Newton
        (f,jf,plane,s.sol.v,k,
         create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
         scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.niter := s.niter + nbit;
      s.nstep := s.nstep + 1;
      if fail then s.nfail := s.nfail + 1; end if;
      Step_Control(fail,pp,step,nbsuccess,s.sol.t,
                   pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      exit when (step <= pp.minstep);
    end loop;
    if (step > pp.minstep) then
      plane := Path(Create(one));
      Projective_LU_Newton
        (f,jf,plane,s.sol.v,k,
         create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
         scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,fail);
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.rcond := to_double(srcond);
    end if;
  end Silent_Projective_LU_Track;

  procedure Reporting_Affine_LU_Track
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    pt0,pt1,pt2,pt3 : Complex_Number := s.sol.t;
    px0,px1,px2,px3 : Vector(s.sol.v'range) := s.sol.v;
    step : quad_double := create(pp.maxstep);
    plane : Matrix(jf'range(2),0..s.sol.v'last);
    fail : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    Process_io.sWrite(file,s.sol.all);
    while not At_End(s,t1,cp.maxtot) or fail loop
      fail := true; -- turns of predictor for solution
      if Process_io.Contains_Output_Code(Process_io.p) then
        Predictor(file,s.sol.t,s.sol.v,step,
                  pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
        Process_io.pWrite(file,s.nstep,step,s.sol.t);
      else
        Predictor(s.sol.t,s.sol.v,step,pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
      end if;
      plane := Path(s.sol.t);
      if Process_io.Contains_Output_Code(Process_io.c) then
        Affine_LU_Newton
          (file,f,jf,plane,s.sol.v,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      else
        Affine_LU_Newton
          (f,jf,plane,s.sol.v,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      end if;
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      Process_io.sWrite(file,s.sol.all);
      s.niter := s.niter + nbit;
      s.nstep := s.nstep + 1;
      if fail then s.nfail := s.nfail + 1; end if;
      if Process_io.Contains_Output_Code(Process_io.p) then
        Step_Control(file,fail,pp,step,nbsuccess,s.sol.t,
                     pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      else
        Step_Control(fail,pp,step,nbsuccess,s.sol.t,
                     pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      end if;
      exit when (step <= pp.minstep);
    end loop;
    if (step > pp.minstep) then
      plane := Path(Create(one));
      if Process_io.Contains_Output_Code(Process_io.c) then
        Affine_LU_Newton
          (file,f,jf,plane,s.sol.v,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,fail);
      else
        Affine_LU_Newton
          (f,jf,plane,s.sol.v,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,fail);
      end if;
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.rcond := to_double(srcond);
    end if;
  end Reporting_Affine_LU_Track;

  procedure Reporting_Projective_LU_Track
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info; k : in natural32;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    pt0,pt1,pt2,pt3 : Complex_Number := s.sol.t;
    px0,px1,px2,px3 : Vector(s.sol.v'range) := s.sol.v;
    step : quad_double := create(pp.maxstep);
    plane : Matrix(jf'range(2),0..s.sol.v'last);
    fail : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    Process_io.sWrite(file,s.sol.all);
    while not At_End(s,t1,cp.maxtot) or fail loop
      if Process_io.Contains_Output_Code(Process_io.p) then
        Predictor(file,s.sol.t,s.sol.v,step,
                  pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
        Process_io.pWrite(file,s.nstep,step,s.sol.t);
      else
        Predictor(s.sol.t,s.sol.v,step,pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
      end if;
      plane := Path(s.sol.t);
      if Process_io.Contains_Output_Code(Process_io.c) then
        Projective_LU_Newton
          (file,f,jf,plane,s.sol.v,k,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      else
        Projective_LU_Newton
          (f,jf,plane,s.sol.v,k,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      end if;
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      Process_io.sWrite(file,s.sol.all);
      s.niter := s.niter + nbit;
      s.nstep := s.nstep + 1;
      if fail then s.nfail := s.nfail + 1; end if;
      if Process_io.Contains_Output_Code(Process_io.p) then
        Step_Control(file,fail,pp,step,nbsuccess,s.sol.t,
                     pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      else
        Step_Control(fail,pp,step,nbsuccess,s.sol.t,
                     pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      end if;
      exit when (step <= pp.minstep);
    end loop;
    if (step > pp.minstep) then
      plane := Path(Create(one));
      if Process_io.Contains_Output_Code(Process_io.c) then
        Projective_LU_Newton
          (file,f,jf,plane,s.sol.v,k,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,fail);
      else
        Projective_LU_Newton
          (f,jf,plane,s.sol.v,k,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,fail);
      end if;
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.rcond := to_double(srcond);
    end if;
  end Reporting_Projective_LU_Track;

  procedure Silent_QR_Track
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    pt0,pt1,pt2,pt3 : Complex_Number := s.sol.t;
    px0,px1,px2,px3 : Vector(s.sol.v'range) := s.sol.v;
    step : quad_double := create(pp.maxstep);
    plane : Matrix(jf'range(2),0..s.sol.v'last);
    fail : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    scora,scorr,sresa,sresr : quad_double;

  begin
    while not At_End(s,t1,cp.maxtot) or fail loop
      Predictor(s.sol.t,s.sol.v,step,pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
      plane := Path(s.sol.t);
      Affine_QR_Newton(f,jf,plane,s.sol.v,
        create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
        scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.niter := s.niter + nbit;
      s.nstep := s.nstep + 1;
      if fail then s.nfail := s.nfail + 1; end if;
      Step_Control(fail,pp,step,nbsuccess,s.sol.t,
                   pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      exit when (step <= pp.minstep);
    end loop;
    if (step > pp.minstep) then
      plane := Path(Create(one));
      declare
        mm : constant integer32 := Min0(f'last+1,plane'last);
        sv : Vector(1..mm);
      begin
        Affine_SV_Newton
          (f,jf,plane,s.sol.v,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,sv,fail);
        s.cora := to_double(scora);
        s.corr := to_double(scorr);
        s.resa := to_double(sresa);
        s.resr := to_double(sresr);
        s.rcond := to_double(Radius(sv(s.sol.v'last)/sv(sv'first)));
      end;
    end if;
  end Silent_QR_Track;

  procedure Reporting_QR_Track
               ( file : in file_type; 
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    pt0,pt1,pt2,pt3 : Complex_Number := s.sol.t;
    px0,px1,px2,px3 : Vector(s.sol.v'range) := s.sol.v;
    step : quad_double := create(pp.maxstep);
    plane : Matrix(jf'range(2),0..s.sol.v'last);
    fail : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    scora,scorr,sresa,sresr : quad_double;

  begin
    Process_io.sWrite(file,s.sol.all);
    while not At_End(s,t1,cp.maxtot) or fail loop
      if Process_io.Contains_Output_Code(Process_io.p) then
        Predictor(file,s.sol.t,s.sol.v,step,
                  pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
        Process_io.pWrite(file,s.nstep,step,s.sol.t);
      else
        Predictor(s.sol.t,s.sol.v,step,pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
      end if;
      plane := Path(s.sol.t);
      if Process_io.Contains_Output_Code(Process_io.c) then
        Affine_QR_Newton
          (file,f,jf,plane,s.sol.v,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      else
        Affine_QR_Newton
          (f,jf,plane,s.sol.v,
           create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
           scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      end if;
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      Process_io.sWrite(file,s.sol.all);
      s.niter := s.niter + nbit;
      s.nstep := s.nstep + 1;
      if fail then s.nfail := s.nfail + 1; end if;
      if Process_io.Contains_Output_Code(Process_io.p) then
        Step_Control(file,fail,pp,step,nbsuccess,s.sol.t,
                     pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      else
        Step_Control(fail,pp,step,nbsuccess,s.sol.t,
                     pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      end if;
      exit when (step <= pp.minstep);
    end loop;
    if (step > pp.minstep) then
      plane := Path(Create(one));
      declare
        mm : constant integer32 := Min0(f'last+1,plane'last);
        sv : Vector(1..mm);
      begin
        if Process_io.Contains_Output_Code(Process_io.c) then
          Affine_SV_Newton
            (file,f,jf,plane,s.sol.v,create(cp.epsax),
             create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
             scora,scorr,sresa,sresr,nbit,cp.maxit,sv,fail);
        else
          Affine_SV_Newton
            (f,jf,plane,s.sol.v,create(cp.epsax),
             create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
             scora,scorr,sresa,sresr,nbit,cp.maxit,sv,fail);
        end if;
        s.cora := to_double(scora);
        s.corr := to_double(scorr);
        s.resa := to_double(sresa);
        s.resr := to_double(sresr);
        s.rcond := to_double(Radius(sv(s.sol.v'last)/sv(sv'first)));
      end;
    end if;
  end Reporting_QR_Track;

-- GENERIC TARGET ROUTINES :

  procedure G_Silent_LU_Track 
               ( ne,nv : in natural32; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    pt0,pt1,pt2,pt3 : Complex_Number := s.sol.t;
    px0,px1,px2,px3 : Vector(s.sol.v'range) := s.sol.v;
    step : quad_double := create(pp.maxstep);
    plane : Matrix(1..integer32(nv),0..s.sol.v'last);
    fail : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    scora,scorr,sresa,sresr,srcond : quad_double;

    procedure Newton is new Silent_Affine_LU_Newton(f,jf);
    procedure RCO_Newton is new Silent_Affine_LU_RCO_Newton(f,jf);

  begin
    while not At_End(s,t1,cp.maxtot) or fail loop
      Predictor(s.sol.t,s.sol.v,step,pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
      plane := Path(s.sol.t);
      Newton(ne,plane,s.sol.v,create(cp.epsax),create(cp.epsrx),
        create(cp.epsaf),create(cp.epsrf),
        scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.niter := s.niter + nbit;
      s.nstep := s.nstep + 1;
      if fail then s.nfail := s.nfail + 1; end if;
      Step_Control(fail,pp,step,nbsuccess,s.sol.t,
                   pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      exit when (step <= pp.minstep);
    end loop;
    if (step > pp.minstep) then
      plane := Path(Create(one));
      RCO_Newton(ne,plane,s.sol.v,create(cp.epsax),create(cp.epsrx),
        create(cp.epsaf),create(cp.epsrf),
        scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,fail);
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.rcond := to_double(srcond);
    end if;
  end G_Silent_LU_Track;

  procedure G_Reporting_LU_Track
               ( file : in file_type; 
                 ne,nv : in natural32; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    pt0,pt1,pt2,pt3 : Complex_Number := s.sol.t;
    px0,px1,px2,px3 : Vector(s.sol.v'range) := s.sol.v;
    step : quad_double := create(pp.maxstep);
    plane : Matrix(1..integer32(nv),0..s.sol.v'last);
    fail : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    scora,scorr,sresa,sresr,srcond : quad_double;

    procedure R_Newton is new Reporting_Affine_LU_Newton(f,jf);
    procedure S_Newton is new Silent_Affine_LU_Newton(f,jf);
    procedure R_RCO_Newton is new Reporting_Affine_LU_RCO_Newton(f,jf);
    procedure S_RCO_Newton is new Silent_Affine_LU_RCO_Newton(f,jf);

  begin
    while not At_End(s,t1,cp.maxtot) or fail loop
      if Process_io.Contains_Output_Code(Process_io.p) then
        Predictor(file,s.sol.t,s.sol.v,step,
                  pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
        Process_io.pWrite(file,s.nstep,step,s.sol.t);
      else
        Predictor(s.sol.t,s.sol.v,step,pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
      end if;
      plane := Path(s.sol.t);
      if Process_io.Contains_Output_Code(Process_io.c) then
        R_Newton(file,ne,plane,s.sol.v,
          create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
          scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      else
        S_Newton(ne,plane,s.sol.v,
          create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
          scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      end if;
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.niter := s.niter + nbit;
      s.nstep := s.nstep + 1;
      if fail then s.nfail := s.nfail + 1; end if;
      if Process_io.Contains_Output_Code(Process_io.p) then
        Step_Control(file,fail,pp,step,nbsuccess,s.sol.t,
                     pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      else
        Step_Control(fail,pp,step,nbsuccess,s.sol.t,
                     pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      end if;
      exit when (step <= pp.minstep);
    end loop;
    if (step > pp.minstep) then
      plane := Path(Create(one));
      if Process_io.Contains_Output_Code(Process_io.c) then
        R_RCO_Newton(file,ne,plane,s.sol.v,create(cp.epsax),
          create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
          scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,fail);
      else
        S_RCO_Newton(ne,plane,s.sol.v,create(cp.epsax),
          create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
          scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,fail);
      end if;
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.rcond := to_double(srcond);
    end if;
  end G_Reporting_LU_Track;

  procedure G_Silent_QR_Track
               ( ne,nv : in natural32; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    pt0,pt1,pt2,pt3 : Complex_Number := s.sol.t;
    px0,px1,px2,px3 : Vector(s.sol.v'range) := s.sol.v;
    step : quad_double := create(pp.maxstep);
    plane : Matrix(1..integer32(nv),0..s.sol.v'last);
    fail : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    scora,scorr,sresa,sresr : quad_double;

    procedure QR_Newton is new Silent_Affine_QR_Newton(f,jf);
    procedure SV_Newton is new Silent_Affine_SV_Newton(f,jf);

  begin
    while not At_End(s,t1,cp.maxtot) or fail loop
      Predictor(s.sol.t,s.sol.v,step,pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
      plane := Path(s.sol.t);
      QR_Newton(ne,plane,s.sol.v,
        create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
        scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.niter := s.niter + nbit;
      s.nstep := s.nstep + 1;
      if fail then s.nfail := s.nfail + 1; end if;
      Step_Control(fail,pp,step,nbsuccess,s.sol.t,
                   pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      exit when (step <= pp.minstep);
    end loop;
    if (step > pp.minstep) then
      plane := Path(Create(one));
      declare
        mm : constant integer32 := Min0(integer32(ne)+1,plane'last);
        sv : Vector(1..mm);
      begin
        SV_Newton(ne,plane,s.sol.v,create(cp.epsax),
          create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
          scora,scorr,sresa,sresr,nbit,cp.maxit,sv,fail);
        s.cora := to_double(scora);
        s.corr := to_double(scorr);
        s.resa := to_double(sresa);
        s.resr := to_double(sresr);
        s.rcond := to_double(Radius(sv(s.sol.v'last)/sv(sv'first)));
      end;
    end if;
  end G_Silent_QR_Track;

  procedure G_Reporting_QR_Track
               ( file : in file_type; 
                 ne,nv : in natural32; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars ) is

    pt0,pt1,pt2,pt3 : Complex_Number := s.sol.t;
    px0,px1,px2,px3 : Vector(s.sol.v'range) := s.sol.v;
    step : quad_double := create(pp.maxstep);
    plane : Matrix(1..integer32(nv),0..s.sol.v'last);
    fail : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    one : constant quad_double := create(1.0);
    t1 : constant Complex_Number := Create(one);
    scora,scorr,sresa,sresr : quad_double;

    procedure R_QR_Newton is new Reporting_Affine_QR_Newton(f,jf);
    procedure S_QR_Newton is new Silent_Affine_QR_Newton(f,jf);
    procedure R_SV_Newton is new Reporting_Affine_SV_Newton(f,jf);
    procedure S_SV_Newton is new Silent_Affine_SV_Newton(f,jf);

  begin
    while not At_End(s,t1,cp.maxtot) or fail loop
      if Process_io.Contains_Output_Code(Process_io.p) then
        Predictor(file,s.sol.t,s.sol.v,step,
                  pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
        Process_io.pWrite(file,s.nstep,step,s.sol.t);
      else
        Predictor(s.sol.t,s.sol.v,step,pt0,pt1,pt2,pt3,t1,px1,px2,px3,fail);
      end if;
      plane := Path(s.sol.t);
      if Process_io.Contains_Output_Code(Process_io.c) then
        R_QR_Newton(file,ne,plane,s.sol.v,
          create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
          scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      else
        S_QR_Newton(ne,plane,s.sol.v,
          create(cp.epsax),create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
          scora,scorr,sresa,sresr,nbit,cp.maxit,fail);
      end if;
      s.cora := to_double(scora);
      s.corr := to_double(scorr);
      s.resa := to_double(sresa);
      s.resr := to_double(sresr);
      s.niter := s.niter + nbit;
      s.nstep := s.nstep + 1;
      if fail then s.nfail := s.nfail + 1; end if;
      if Process_io.Contains_Output_Code(Process_io.p) then
        Step_Control(file,fail,pp,step,nbsuccess,s.sol.t,
                     pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      else
        Step_Control(fail,pp,step,nbsuccess,s.sol.t,
                     pt0,pt1,pt2,pt3,s.sol.v,px0,px1,px2,px3);
      end if;
      exit when (step <= pp.minstep);
    end loop;
    if (step > pp.minstep) then
      plane := Path(Create(one));
      declare
        mm : constant integer32 := Min0(integer32(ne)+1,plane'last);
        sv : Vector(1..mm);
      begin
        if Process_io.Contains_Output_Code(Process_io.c) then
          R_SV_Newton
            (file,ne,plane,s.sol.v,create(cp.epsax),
             create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
             scora,scorr,sresa,sresr,nbit,cp.maxit,sv,fail);
        else
          S_SV_Newton(ne,plane,s.sol.v,create(cp.epsax),
            create(cp.epsrx),create(cp.epsaf),create(cp.epsrf),
            scora,scorr,sresa,sresr,nbit,cp.maxit,sv,fail);
        end if;
        s.cora := to_double(scora);
        s.corr := to_double(scorr);
        s.resa := to_double(sresa);
        s.resr := to_double(sresr);
        s.rcond := to_double(Radius(sv(s.sol.v'last)/sv(sv'first)));
      end;
    end if;
  end G_Reporting_QR_Track;

-- VERSIONS USING LOCAL COORDINATES :

  function Predict_Offset
               ( p : Matrix; b,x : Vector; d,h : quad_double;
                 step : Complex_Number; fixed : boolean ) return Vector is

  -- DESCRIPTION :
  --   Returns the predicted offset vector of the k-plane,
  --   using extrinsic coordinates of the solution in x.

  -- ON ENTRY :
  --   p         representation of a k-plane in n-space;
  --   b         target offset vector;
  --   x         extrinsic coordinates of the current solution;
  --   d         distance to the target offset vector;
  --   h         current step size;
  --   step      step size as a complex number;
  --   fixed     flag to ignore distance and step size.

    res : Vector(b'range);
    to_b : constant Vector(b'range) := b - x(b'range);
    c : Vector(b'range) := Complement_of_Projection(p,to_b);

  begin
    if fixed then
      for i in c'range loop
        res(i) := x(i) + c(i);
      end loop;
    else
      if d < h then
        for i in c'range loop
          res(i) := x(i) + c(i);
        end loop;
      else
        Normalize(c);
        for i in c'range loop
          res(i) := x(i) + step*c(i);
        end loop;
      end if;
    end if;
    return res;
  end Predict_Offset;

  function Evaluate_Prediction
              ( f : Eval_Poly_Sys; x : Vector ) return quad_double is

  -- DESCRIPTION :
  --   Returns the 2-norm of f evaluated at x.
  --   The prediction will turn right if the number on return is
  --   of the same magnitude of the step size.

    y : constant Vector := Eval(f,x);
    res : constant quad_double := Norm2(y);

  begin
    return res;
  end Evaluate_Prediction;

  procedure Predict_Offset
               ( p : in out Matrix; f : in Eval_Poly_Sys; b,x : in Vector;
                 d,rf,m : in quad_double; fixed : in boolean;
                 h : in out quad_double;
                 step : in out Complex_Number; nb : in natural32 ) is

  -- DESCRIPTION :
  --   Predicts the offset vector of the k-plane, using the extrinsic
  --   coordinates of the solution in x.

  -- ON ENTRY :
  --   p         representation of a k-plane in n-space;
  --   f         polynomial system to evaluate the prediction;
  --   b         target offset vector;
  --   x         extrinsic coordinates of the current solution;
  --   d         distance to the target offset vector;
  --   rf        reduction factor to cut step size with;
  --   m         maximal step size;
  --   fixed     to ignore distance and step size;
  --   h         current step size;
  --   step      step size as a complex number;
  --   nb        number of iterations of previous Newton correction.

  -- ON RETURN :
  --   p         offset vector is the predicted solution;
  --   h         step size is adjusted according to the evaluation,
  --             will only be cut if nb > 2;
  --   step      corresponding complex value of the step size.

    pb : Vector(b'range) := Predict_Offset(p,b,x,d,h,step,fixed);
    v : constant quad_double := Evaluate_Prediction(f,pb);
    r : constant quad_double := v/m;

  begin
    if r > 1.0 and d > 0.1 then -- and nb > 2 then
      h := h*rf;
      step := Create(h);
      if not fixed
       then pb := Predict_Offset(p,b,x,d,h,step,false);
      end if;
    end if;
    for i in p'range(1) loop
      p(i,0) := pb(i);
    end loop;
  end Predict_Offset;

  procedure Predict_Offset
               ( file : in file_type;
                 p : in out Matrix; f : in Eval_Poly_Sys; b,x : in Vector;
                 d,rf,m : in quad_double; fixed : in boolean;
                 h : in out quad_double;
                 step : in out Complex_Number; nb : in natural32 ) is

  -- DESCRIPTION :
  --   Predicts the offset vector of the k-plane, using the extrinsic
  --   coordinates of the solution in x.

  -- ON ENTRY :
  --   file      for report on the evaluation;
  --   p         representation of a k-plane in n-space;
  --   f         polynomial system to evaluate the prediction;
  --   b         target offset vector;
  --   x         extrinsic coordinates of the current solution;
  --   d         distance to the target offset vector;
  --   rf        reduction factor to cut step size with;
  --   m         maximum step size;
  --   fixed     flag to ignore distance and prediction;
  --   h         current step size;
  --   step      step size as a complex number;
  --   nb        number of iterations of previous Newton correction.

  -- ON RETURN :
  --   p         offset vector is the predicted solution;
  --   h         step size is adjusted according to the evaluation,
  --             will only be cut if nb > 2;
  --   step      corresponding complex value of the step size.

    pb : Vector(b'range) := Predict_Offset(p,b,x,d,h,step,fixed);
    v : constant quad_double := Evaluate_Prediction(f,pb);
    r : constant quad_double := v/m;

  begin
    if Process_io.Contains_Output_Code(Process_io.p) then
      put(file,"value of prediction : "); put(file,v,3);
      put(file,"  ratio with step : "); put(file,r,3);
    end if;
    if r > 1.0 and d > 0.1 then -- and nb > 2 then
      if Process_io.Contains_Output_Code(Process_io.p)
       then put_line(file,"  shorten step");
      end if;
      h := h*rf;
      step := Create(h);
      if not fixed
       then pb := Predict_Offset(p,b,x,d,h,step,false);
      end if;
    else
      if Process_io.Contains_Output_Code(Process_io.p)
       then put_line(file,"  step size ok");
      end if;
    end if;
    for i in p'range(1) loop
      p(i,0) := pb(i);
    end loop;
  end Predict_Offset;

  procedure Update_from_Offset
               ( x : in out Vector; p : in Matrix; z : in Vector ) is

  -- DESCRIPTION :
  --   Updates the extrinsic coordinates for the solution, using the
  --   offset vector of p and a linear combination of the vector of p,
  --   multiplied with the local intrinsic coordinates of z.

  begin
    for i in p'range(1) loop
      x(i) := p(i,0);
      for j in z'range loop
        x(i) := x(i) + z(j)*p(i,j);
      end loop;
    end loop;
  end Update_from_Offset;

  procedure Silent_Local_LU_Track
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 p : in Matrix; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars;
                 fail : out boolean ) is

    k : constant integer32 := p'last(2);
    step : quad_double := create(pp.maxstep);
    h : Complex_Number := Create(step);
    plane : Matrix(p'range(1),p'range(2)) := p;
    b : Vector(p'range(1));
    nfl : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    z : Vector(1..k);
    d,scora,scorr,sresa,sresr,srcond : quad_double;
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in b'range loop
      b(i) := plane(i,0);
    end loop;
    s.nstep := 0;
    loop
      d := Distance(p,s.sol.v);
      exit when (d <= 1.0E-12);
      Predict_Offset
        (plane,f,b,s.sol.v,d,create(pp.redfac),create(pp.maxstep),false,
         step,h,nbit);
      z := (1..k => Create(zero));
      Affine_LU_Newton(f,jf,plane,z,
        create(cp.epsax),one,create(cp.epsaf),one,
        scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,nfl);
        s.cora := to_double(scora);
        s.corr := to_double(scorr);
        s.resa := to_double(sresa);
        s.resr := to_double(sresr);
        s.rcond := to_double(srcond);
      if not nfl
       then Update_from_Offset(s.sol.v,plane,z);
      end if;
      s.niter := s.niter + nbit; s.nstep := s.nstep + 1;
      exit when (s.niter >= cp.maxtot);
      Step_Control(nfl,pp,step,nbsuccess);
      exit when (step <= pp.minstep);
      h := Create(step);
    end loop;
    fail := (d > 1.0E-12) or ((s.corr > cp.epsrx) and (s.resa > cp.epsaf));
  end Silent_Local_LU_Track;

  procedure Reporting_Local_LU_Track
               ( file : in file_type;
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 p : in Matrix; s : in out Solu_Info;
                 pp : in Pred_Pars; cp : in Corr_Pars;
                 fail : out boolean ) is

    k : constant integer32 := p'last(2);
    step : quad_double := create(pp.maxstep);
    h : Complex_Number := Create(step);
    plane : Matrix(p'range(1),p'range(2)) := p;
    b : Vector(p'range(1));
    nfl : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    z : Vector(1..k);
    d : quad_double;
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    for i in b'range loop
      b(i) := p(i,0);
    end loop;
    s.nstep := 0;
    loop
      d := Distance(p,s.sol.v);
      if Process_io.Contains_Output_Code(Process_io.p) then
        put(file,"step "); put(file,s.nstep+1,1);
        put(file,"  distance to target : "); put(file,d,3);
        put(file,"  step size : "); put(file,step,3); new_line(file);
      end if;
      exit when (d < 1.0E-12);
      Predict_Offset
        (file,plane,f,b,s.sol.v,d,create(pp.redfac),create(pp.maxstep),false,
         step,h,nbit);
      z := (1..k => Create(zero));
      Affine_LU_Newton(file,f,jf,plane,z,
        create(cp.epsax),one,create(cp.epsaf),one,
        scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,nfl);
        s.cora := to_double(scora);
        s.corr := to_double(scorr);
        s.resa := to_double(sresa);
        s.resr := to_double(sresr);
        s.rcond := to_double(srcond);
      if not nfl
       then Update_from_Offset(s.sol.v,plane,z);
      end if;
      s.niter := s.niter + nbit; s.nstep := s.nstep + 1;
      exit when (s.niter >= cp.maxtot);
      Step_Control(file,nfl,pp,step,nbsuccess);
      exit when (step <= pp.minstep);
      h := Create(step);
    end loop;
    fail := (d > 1.0E-12) or ((s.corr > cp.epsrx) and (s.resa > cp.epsaf));
  end Reporting_Local_LU_Track;

  procedure Silent_Recentered_LU_Track
               ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 start,target : in Matrix; reoriented : in boolean;
                 s : in out Solu_Info; pp : in Pred_Pars; cp : in Corr_Pars;
                 fail : out boolean ) is

    k : constant integer32 := target'last(2);
    step : quad_double := create(pp.maxstep);
    inner_step : quad_double;
    h : Complex_Number := Create(step);
    plane : Matrix(target'range(1),target'range(2));
    b : Vector(target'range(1));
    nfl : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    z : Vector(1..k);
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);
    previous_rt,rt,d : quad_double := zero;
    t : Complex_Number := Create(zero);
    scora,scorr,sresa,sresr,srcond : quad_double;

  begin
    s.nstep := 0; s.nfail := 0; rt := zero;
    loop
      rt := rt + step;
      if rt > 1.0
       then step := to_double(rt)-1.0; rt := one; h := Create(step);
      end if;
      t := Create(rt);
      plane := QuadDobl_Moving_Planes.Moving_Plane(start,target,t);
     -- Assuming only one vector moves, orthogonalization not needed
     -- but we must normalize the last moving vector!
      if reoriented then
        if rt /= one
         then QuadDobl_Moving_Planes.Normalize(plane,plane'last(2));
        end if;
      else
        if rt /= one
         then QuadDobl_Plane_Representations.Orthogonalize(plane);
        end if;
      end if;
      for i in b'range loop
        b(i) := plane(i,0);
      end loop;
      d := Distance(plane,s.sol.v);
      if d < 1.0E-13 then
        nfl := false;
      else
        inner_step := create(1.0);
        h := Create(inner_step);
        Predict_Offset
          (plane,f,b,s.sol.v,d,create(pp.redfac),create(pp.maxstep),true,
           inner_step,h,nbit);
        if inner_step < 1.0 then
          nfl := true;
        else
          z := (1..k => Create(zero));
          Affine_LU_Newton(f,jf,plane,z,
            create(cp.epsax),one,create(cp.epsaf),one,
            scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,nfl);
            s.cora := to_double(scora);
            s.corr := to_double(scorr);
            s.resa := to_double(sresa);
            s.resr := to_double(sresr);
            s.rcond := to_double(srcond);
          s.niter := s.niter + nbit; s.nstep := s.nstep + 1;
          if not nfl
           then Update_from_Offset(s.sol.v,plane,z);
          end if;
        end if;
      end if;
      if nfl
       then rt := previous_rt; s.nfail := s.nfail + 1;
       else previous_rt := rt;
      end if;
      exit when (s.niter >= cp.maxtot);
      exit when (not nfl) and (abs(rt - 1.0) < 1.0E-13);
      Step_Control(standard_output,nfl,pp,step,nbsuccess);
      exit when (step <= pp.minstep);
      h := Create(step);
    end loop;
    fail := (rt < create(1.0-1.0E-12))
         or ((s.corr > cp.epsrx) and (s.resa > cp.epsaf));
  end Silent_Recentered_LU_Track;

  procedure Reporting_Recentered_LU_Track
               ( file : in file_type;
                 f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                 start,target : in Matrix; reoriented : in boolean;
                 s : in out Solu_Info; pp : in Pred_Pars; cp : in Corr_Pars;
                 fail : out boolean ) is

    k : constant integer32 := target'last(2);
    step : quad_double := create(pp.maxstep);
    inner_step : quad_double;
    h : Complex_Number := Create(step);
    plane : Matrix(target'range(1),target'range(2));
    b : Vector(target'range(1));
    nfl : boolean := true;
    nbit,nbsuccess : natural32 := 0;
    z : Vector(1..k);
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);
    previous_rt,rt,d : quad_double := zero;
    t : Complex_Number := Create(zero);
    scora,scorr,sresa,sresr,srcond : quad_double;
   -- nrm : quad_double;
   -- fl : boolean;

  begin
    s.nstep := 0; s.nfail := 0; rt := zero;
    loop
      rt := rt + step;
      if rt > 1.0
       then step := to_double(rt) - 1.0; rt := one; h := Create(step);
      end if;
      t := Create(rt);
      plane := QuadDobl_Moving_Planes.Moving_Plane(start,target,t);
     -- Assuming only one vector moves, orthogonalization not needed
     -- but we must normalize the last moving vector!
      if reoriented then
        if rt /= one then
          QuadDobl_Moving_Planes.Normalize(plane,plane'last(2));
          -- QuadDobl_Moving_Planes.Normalize(plane,plane'last(2),nrm);
          -- put(file,"norm of moving vector : ");
          -- put(file,nrm,3); new_line(file);
        end if;
      else
        if rt /= one 
         then QuadDobl_Plane_Representations.Orthogonalize(plane);
        end if;
      end if;
     -- Check_Orthonormality(plane,1.0E-12,fl);
     -- if fl
     --  then put_line(file,"Orthogonality check of plane failed!");
     -- end if;
      for i in b'range loop
        b(i) := plane(i,0);
      end loop;
      if Process_io.Contains_Output_Code(Process_io.p) then
        put(file,"outer step "); put(file,s.nstep+1,1);
        put(file,"  t : "); put(file,rt,3); new_line(file);
      end if;
      d := Distance(plane,s.sol.v);
      if Process_io.Contains_Output_Code(Process_io.p) then
        put(file,"inner step "); put(file,s.nstep+1,1);
        put(file,"  distance : "); put(file,d,3);
        put(file,"  step size : "); put(file,step,3); new_line(file);
      end if;
      if d < 1.0E-13 then
        nfl := false;
        if Process_io.Contains_Output_Code(Process_io.p)
         then put_line(file,"zero distance so no predictor-corrector");
        end if;
      else
        inner_step := one; h := Create(inner_step);
        Predict_Offset -- step size must equal one !
          (file,plane,f,b,s.sol.v,d,create(pp.redfac),create(pp.maxstep),true,
           inner_step,h,nbit);
        if inner_step < 1.0 then
          nfl := true;
          if Process_io.Contains_Output_Code(Process_io.p)
           then put_line(file,"inner step prediction failed");
          end if;
        else
          z := (1..k => Create(zero));
          Affine_LU_Newton(file,f,jf,plane,z,
            create(cp.epsax),one,create(cp.epsaf),one,
            scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,nfl);
            s.cora := to_double(scora);
            s.corr := to_double(scorr);
            s.resa := to_double(sresa);
            s.resr := to_double(sresr);
            s.rcond := to_double(srcond);
          if not nfl then
            Update_from_Offset(s.sol.v,plane,z);
            if Process_io.Contains_Output_Code(Process_io.c) then
              d := Distance(plane,s.sol.v);
              put(file,"distance after Newton : ");
              put(file,d,3); new_line(file);
            end if;
          end if;
          s.niter := s.niter + nbit; s.nstep := s.nstep + 1;
        end if;
      end if;
      exit when (s.niter >= cp.maxtot);
      exit when (not nfl) and (abs(rt - 1.0) < 1.0E-13);
      if nfl
       then rt := previous_rt; s.nfail := s.nfail + 1;
       else previous_rt := rt;
      end if;
      Step_Control(file,nfl,pp,step,nbsuccess);
      exit when (step <= pp.minstep);
      h := Create(step);
    end loop;
    plane := target;
    for i in plane'range(1) loop -- ignore slack variables ...
      plane(i,0) := s.sol.v(i);
    end loop;
    z := (1..k => Create(zero));
    if Process_io.Contains_Output_Code(Process_io.c) then
      put_line(file,"running sanity Newton check at end of tracking...");
    end if;
    Affine_LU_Newton(file,f,jf,plane,z,
       create(cp.epsax),one,create(cp.epsaf),one,
       scora,scorr,sresa,sresr,nbit,cp.maxit,srcond,nfl);
       s.cora := to_double(scora);
       s.corr := to_double(scorr);
       s.resa := to_double(sresa);
       s.resr := to_double(sresr);
       s.rcond := to_double(srcond);
    fail := (rt < create(1.0-1.0E-12))
         or ((s.corr > cp.epsrx) and (s.resa > cp.epsaf));
    d := Distance(target,s.sol.v);
    put(file,"distance to target plane : "); put(file,d,3); new_line(file);
  end Reporting_Recentered_LU_Track;

end QuadDobl_Intrinsic_Trackers;
