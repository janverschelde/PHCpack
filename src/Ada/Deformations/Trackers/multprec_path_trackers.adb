with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Multprec_Floating_Vectors_io;       use Multprec_Floating_Vectors_io;
with Multprec_Floating_VecVecs;          use Multprec_Floating_VecVecs;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Predictors;                use Multprec_Predictors;
with Multprec_Correctors;                use Multprec_Correctors;
with Multprec_Dispatch_Predictors;       use Multprec_Dispatch_Predictors;
with Continuation_Parameters;
with Process_io;                         use Process_io;
with Multprec_Data_on_Path;              use Multprec_Data_on_Path;

package body Multprec_Path_Trackers is

-- LINEAR PATH FOLLOWING FOR ONE PATH :

  procedure Linear_Single_Normal_Silent_Continue
              ( s : in out Solu_Info; target : in Complex_Number;
                tol : in Floating_Number; proj : in boolean;
                p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t,prev_t2,prev_t1,prev_t0 : Complex_Number;
    old_sol,prev_sol,prev_sol2,prev_sol1,prev_sol0 : Vector(s.sol.v'range);
    old_v,prev_v,vv : Vector(s.sol.v'range);
    step : Floating_Number;
    nsuccess,trial : natural32 := 0;
    success : boolean := true;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);

    begin
      if proj
       then Projective_Corrector(s,c);
       else Affine_Corrector(s,c);
      end if;
    end Corrector;

  begin
    Copy(p.maxstep,step);
    Linear_Single_Initialize(s,p,old_t,prev_t,prev_v,old_sol,prev_sol);
    if p.predictor_type > 6 then
      Copy(prev_t,prev_t0); Copy(prev_sol,prev_sol0);
      Copy(prev_t,prev_t1); Copy(prev_sol,prev_sol1);
      Copy(prev_t,prev_t2); Copy(prev_sol,prev_sol2);
    end if;
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success) 
                                           and (s.niter <= c.maxtot) loop
      if p.predictor_type < 7 then
        Predictor(s,p,true,prev_sol,prev_v,vv,prev_t,target,step,tol,trial);
      elsif p.predictor_type = 7 then
        Single_Quadratic_Predictor
          (s,p,true,prev_sol,prev_sol0,prev_t,prev_t0,target,step,tol);
      elsif p.predictor_type = 8 then
        Single_Cubic_Predictor
          (s,p,true,prev_sol,prev_sol1,prev_sol0,prev_t,prev_t1,prev_t0,
           target,step,tol);
      else
        Single_Quartic_Predictor
          (s,p,true,prev_sol,prev_sol2,prev_sol1,prev_sol0,prev_t,prev_t2,
           prev_t1,prev_t0,target,step,tol);
      end if;
      Corrector(s,c);
      if p.predictor_type < 7 then
        Linear_Single_Management
          (s,p,c,old_t,prev_t,old_sol,prev_sol,old_v,prev_v,vv,
           step,nsuccess,trial,success);
      elsif p.predictor_type = 7 then
        Linear_Single_Quadratic_Management
          (s,p,c,old_t,prev_t,prev_t0,old_sol,prev_sol,prev_sol0,
           step,nsuccess,trial,success);
      elsif p.predictor_type = 8 then
        Linear_Single_Cubic_Management
          (s,p,c,old_t,prev_t,prev_t1,prev_t0,old_sol,prev_sol,prev_sol1,
           prev_sol0,step,nsuccess,trial,success);
      else
        Linear_Single_Quartic_Management
          (s,p,c,old_t,prev_t,prev_t2,prev_t1,prev_t0,old_sol,prev_sol,
           prev_sol2,prev_sol1,prev_sol0,step,nsuccess,trial,success);
      end if;
      if Stop(p,s.sol.t,target,step) then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Corrector(s,cp);
    end;
    Clear(step);
    Clear(old_t); Clear(prev_t);
    Clear(prev_t0); Clear(prev_t1); Clear(prev_t2);
    Clear(old_sol); Clear(prev_sol);
    Clear(prev_sol0); Clear(prev_sol1); Clear(prev_sol2);
    Clear(old_v); Clear(prev_v); Clear(vv);
  end Linear_Single_Normal_Silent_Continue;

  procedure Linear_Single_Normal_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in Complex_Number; tol : in Floating_Number;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t,prev_t2,prev_t1,prev_t0 : Complex_Number;
    old_sol,prev_sol,prev_sol2,prev_sol1,prev_sol0 : Vector(s.sol.v'range);
    old_v,prev_v,vv : Vector(s.sol.v'range);
    step : Floating_Number;
    nsuccess,trial : natural32 := 0;
    success : boolean := true;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);

    begin
      if proj
       then Projective_Corrector(file,s,c);
       else Affine_Corrector(file,s,c);
      end if;
    end Corrector;

  begin
    Copy(p.maxstep,step);
    Linear_Single_Initialize(s,p,old_t,prev_t,prev_v,old_sol,prev_sol);
    if p.predictor_type > 6 then
      Copy(prev_t,prev_t0); Copy(prev_sol,prev_sol0);
      Copy(prev_t,prev_t1); Copy(prev_sol,prev_sol1);
      Copy(prev_t,prev_t2); Copy(prev_sol,prev_sol2);
    end if;
    sWrite(file,s.sol.all);
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success) 
                                           and (s.niter <= c.maxtot) loop
      if p.predictor_type < 7 then
        Predictor(s,p,true,prev_sol,prev_v,vv,prev_t,target,step,tol,trial);
      elsif p.predictor_type = 7 then
        Single_Quadratic_Predictor
          (s,p,true,prev_sol,prev_sol0,prev_t,prev_t0,target,step,tol);
      elsif p.predictor_type = 8 then
        Single_Cubic_Predictor
          (s,p,true,prev_sol,prev_sol1,prev_sol0,prev_t,prev_t1,prev_t0,
           target,step,tol);
      else
        Single_Quartic_Predictor
          (s,p,true,prev_sol,prev_sol2,prev_sol1,prev_sol0,prev_t,prev_t2,
           prev_t1,prev_t0,target,step,tol);
      end if;
     -- if p.predictor_type = 6
     --  then put_line(file,"previous and current direction : "); 
     --       for i in prev_v'range loop
     --         put(file,prev_v(i),3,3,3);  put(file,"   ");
     --         put(file,vv(i),3,3,3);      new_line(file);
     --       end loop;
     -- end if;
      pWrite(file,step,s.sol.t,s.sol.all);
      Corrector(s,c);
      sWrite(file,s.sol.all);
      if p.predictor_type < 7 then
        Linear_Single_Management
          (s,p,c,old_t,prev_t,old_sol,prev_sol,old_v,prev_v,vv,
           step,nsuccess,trial,success);
      elsif p.predictor_type = 7 then
        Linear_Single_Quadratic_Management
          (s,p,c,old_t,prev_t,prev_t0,old_sol,prev_sol,prev_sol0,
           step,nsuccess,trial,success);
      elsif p.predictor_type = 8 then
        Linear_Single_Cubic_Management
          (s,p,c,old_t,prev_t,prev_t1,prev_t0,old_sol,prev_sol,prev_sol1,
           prev_sol0,step,nsuccess,trial,success);
      else
        Linear_Single_Quartic_Management
          (s,p,c,old_t,prev_t,prev_t2,prev_t1,prev_t0,old_sol,prev_sol,
           prev_sol2,prev_sol1,prev_sol0,step,nsuccess,trial,success);
      end if;
      if Stop(p,s.sol.t,target,step) then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Corrector(s,cp);
      sWrite(file,s.sol.all);
    end;
    Clear(step);
    Clear(old_t); Clear(prev_t);
    Clear(prev_t0); Clear(prev_t1); Clear(prev_t2);
    Clear(old_sol); Clear(prev_sol);
    Clear(prev_sol0); Clear(prev_sol1); Clear(prev_sol2);
    Clear(old_v); Clear(prev_v); Clear(vv);
  end Linear_Single_Normal_Reporting_Continue;

  procedure Linear_Single_Conditioned_Silent_Continue
              ( s : in out Solu_Info; target : in Complex_Number;
                tol : in Floating_Number; proj : in boolean;
                rtoric : in integer32; 
                v : in out Multprec_Floating_Vectors.Link_to_Vector;
                errorv : in out Floating_Number;
                p : in Pred_Pars; c : in Corr_Pars ) is

    old_t,prev_t,prev_t2,prev_t1,prev_t0 : Complex_Number;
    old_sol,prev_sol,prev_sol2,prev_sol1,prev_sol0 : Vector(s.sol.v'range);
    old_v,prev_v,vv : Vector(s.sol.v'range);
    step : Floating_Number;
    nsuccess,trial : natural32 := 0;
    success : boolean := true;

    dls,stp : Floating_Number := Create(integer(0));
    tolsing : Floating_Number
            := Create(Continuation_Parameters.tol_endg_inverse_condition);
    diff : Multprec_Floating_Vectors.Vector(s.sol.v'range)
         := (s.sol.v'range => Create(integer(0)));

    r : integer32 := 0;
    er : integer32 := 0;
    dt,ds,logs : Multprec_Floating_Vectors.Vector(0..rtoric)
               := (0..rtoric => Create(integer(0)));
    logx : VecVec(0..rtoric);
    wvl0,wvl1,wvl2 : VecVec(1..rtoric);
    errv : Multprec_Floating_Vectors.Vector(0..rtoric)
         := (0..rtoric => Create(integer(0)));

    m : integer32 := 1;
    thresm : natural32 := p.success_steps;
    estm : integer32 := m;
    fcnt,cntm : natural32 := 0;
    defer : natural32 := thresm;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);

    begin
      if proj
       then Projective_Corrector(s,c);
       else Affine_Corrector(s,c);
      end if;
    end Corrector;

  begin
    Copy(p.maxstep,step);
    Linear_Single_Initialize(s,p,old_t,prev_t,prev_v,old_sol,prev_sol);
    if p.predictor_type > 6 then
      Copy(prev_t,prev_t0); Copy(prev_sol,prev_sol0);
      Copy(prev_t,prev_t1); Copy(prev_sol,prev_sol1);
      Copy(prev_t,prev_t2); Copy(prev_sol,prev_sol2);
    end if;
    if rtoric > 0
     then s.sol.m := m; 
    end if;
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success)
                                           and (s.niter <= c.maxtot) loop
      if (rtoric > 0) then 
        if success and then s.rcond > tolsing
                   and then (errorv < 100.0) -- avoid divergence
         then Update_Direction
                (proj,fcnt,defer,r,s.sol.m,estm,cntm,thresm,er,s.sol.t,target,
                 s.sol.v,dt,ds,logs,logx,wvl0,wvl1,wvl2,v.all,errv,errorv);
         else er := -2;
        end if;
      end if;
      if p.predictor_type < 7 then
        Predictor(s,p,true,prev_sol,prev_v,vv,prev_t,target,step,tol,trial);
      elsif p.predictor_type = 7 then
        Single_Quadratic_Predictor
          (s,p,true,prev_sol,prev_sol0,prev_t,prev_t0,target,step,tol);
      elsif p.predictor_type = 8 then
        Single_Cubic_Predictor
          (s,p,true,prev_sol,prev_sol1,prev_sol0,prev_t,prev_t1,prev_t0,
           target,step,tol);
      else
        Single_Quartic_Predictor
          (s,p,true,prev_sol,prev_sol2,prev_sol1,prev_sol0,prev_t,prev_t2,
           prev_t1,prev_t0,target,step,tol);
      end if;
      Corrector(s,c);
      if p.predictor_type < 7 then
        Linear_Single_Management
          (s,p,c,old_t,prev_t,old_sol,prev_sol,old_v,prev_v,vv,
           step,nsuccess,trial,success);
      elsif p.predictor_type = 7 then
        Linear_Single_Quadratic_Management
          (s,p,c,old_t,prev_t,prev_t0,old_sol,prev_sol,prev_sol0,
           step,nsuccess,trial,success);
      elsif p.predictor_type = 8 then
        Linear_Single_Cubic_Management
          (s,p,c,old_t,prev_t,prev_t1,prev_t0,old_sol,prev_sol,prev_sol1,
           prev_sol0,step,nsuccess,trial,success);
      else
        Linear_Single_Quartic_Management
          (s,p,c,old_t,prev_t,prev_t2,prev_t1,prev_t0,old_sol,prev_sol,
           prev_sol2,prev_sol1,prev_sol0,step,nsuccess,trial,success);
      end if;
      if Stop(p,s.sol.t,target,step) then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Corrector(s,cp);
    end;
    Clear(step);
    Clear(old_t); Clear(prev_t);
    Clear(prev_t0); Clear(prev_t1); Clear(prev_t2);
    Clear(old_sol); Clear(prev_sol);
    Clear(prev_sol0); Clear(prev_sol1); Clear(prev_sol2);
    Clear(old_v); Clear(prev_v); Clear(vv);
  end Linear_Single_Conditioned_Silent_Continue;

  procedure Linear_Single_Conditioned_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in Complex_Number; tol : in Floating_Number;
                proj : in boolean; rtoric : in integer32;
                v : in out Multprec_Floating_Vectors.Link_to_Vector;
                errorv : in out Floating_Number;
                p : in Pred_Pars; c : in Corr_Pars ) is

    old_t,prev_t,prev_t2,prev_t1,prev_t0 : Complex_Number;
    step : Floating_Number;
    nsuccess,trial : natural32 := 0;
    old_sol,prev_sol,prev_sol2,prev_sol1,prev_sol0 : Vector(s.sol.v'range);
    old_v,prev_v,vv : Vector(s.sol.v'range);
    success : boolean := true;

    dls,stp : Floating_Number := Create(integer(0));
    tolsing : Floating_Number
            := Create(Continuation_Parameters.tol_endg_inverse_condition);
    diff : Multprec_Floating_Vectors.Vector(s.sol.v'range)
         := (s.sol.v'range => Create(integer(0)));

    r : integer32 := 0;
    er : integer32 := 0;
    dt,ds,logs : Multprec_Floating_Vectors.Vector(0..rtoric)
               := (0..rtoric => Create(integer(0)));
    logx : VecVec(0..rtoric);
    wvl0,wvl1,wvl2 : VecVec(1..rtoric);
    errv : Multprec_Floating_Vectors.Vector(0..rtoric)
         := (0..rtoric => Create(integer(0)));

    m : integer32 := 1;
    thresm : natural32 := p.success_steps;
    estm : integer32 := m;
    fcnt,cntm : natural32:= 0;
    defer : natural32 := thresm;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);

    begin
      if proj
       then Projective_Corrector(file,s,c);
       else Affine_Corrector(file,s,c);
      end if;
    end Corrector;

  begin
    Copy(p.maxstep,step);
    Linear_Single_Initialize(s,p,old_t,prev_t,prev_v,old_sol,prev_sol);
    sWrite(file,s.sol.all);          -- writing the start solution
    if p.predictor_type > 6 then
      Copy(prev_t,prev_t0); Copy(prev_sol,prev_sol0);
      Copy(prev_t,prev_t1); Copy(prev_sol,prev_sol1);
      Copy(prev_t,prev_t2); Copy(prev_sol,prev_sol2);
    end if;
    if rtoric > 0
     then s.sol.m := m;
    end if;
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success)
                                           and (s.niter <= c.maxtot) loop
      if (rtoric > 0) then
        if success and then s.rcond > tolsing
                   and then (errorv < 100.0) -- avoid divergence
         then Update_Direction(file,
                 proj,fcnt,defer,r,s.sol.m,estm,cntm,thresm,er,s.sol.t,target,
                 s.sol.v,dt,ds,logs,logx,wvl0,wvl1,wvl2,v.all,errv,errorv);
         else er := -2;
        end if;
      end if;
      if p.predictor_type < 7 then
        Predictor(s,p,true,prev_sol,prev_v,vv,prev_t,target,step,tol,trial);
      elsif p.predictor_type = 7 then
        Single_Quadratic_Predictor
          (s,p,true,prev_sol,prev_sol0,prev_t,prev_t0,target,step,tol);
      elsif p.predictor_type = 8 then
        Single_Cubic_Predictor
          (s,p,true,prev_sol,prev_sol1,prev_sol0,prev_t,prev_t1,prev_t0,
           target,step,tol);
      else
        Single_Quartic_Predictor
          (s,p,true,prev_sol,prev_sol2,prev_sol1,prev_sol0,prev_t,prev_t2,
           prev_t1,prev_t0,target,step,tol);
      end if;
      pWrite(file,step,s.sol.t,s.sol.all);
      Corrector(s,c);
      sWrite(file,s.sol.all);
      if p.predictor_type < 7 then
        Linear_Single_Management
          (s,p,c,old_t,prev_t,old_sol,prev_sol,old_v,prev_v,vv,
           step,nsuccess,trial,success);
      elsif p.predictor_type = 7 then
        Linear_Single_Quadratic_Management
          (s,p,c,old_t,prev_t,prev_t0,old_sol,prev_sol,prev_sol0,
           step,nsuccess,trial,success);
      elsif p.predictor_type = 8 then
        Linear_Single_Cubic_Management
          (s,p,c,old_t,prev_t,prev_t1,prev_t0,old_sol,prev_sol,prev_sol1,
           prev_sol0,step,nsuccess,trial,success);
      else
        Linear_Single_Quartic_Management
          (s,p,c,old_t,prev_t,prev_t2,prev_t1,prev_t0,old_sol,prev_sol,
           prev_sol2,prev_sol1,prev_sol0,step,nsuccess,trial,success);
      end if;
      if Stop(p,s.sol.t,target,step) then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Corrector(s,cp);
      sWrite(file,s.sol.all);
    end;
    Clear(step);
    Clear(old_t); Clear(prev_t);
    Clear(prev_t0); Clear(prev_t1); Clear(prev_t2);
    Clear(old_sol); Clear(prev_sol);
    Clear(prev_sol0); Clear(prev_sol1); Clear(prev_sol2);
    Clear(old_v); Clear(prev_v); Clear(vv);
  end Linear_Single_Conditioned_Reporting_Continue;

-- LINEAR PATH FOLLOWING FOR A NUMBER OF PATHS :

  procedure Linear_Multiple_Normal_Silent_Continue
              ( s : in out Solu_Info_Array;
                target : in Complex_Number; tol,dist_sols : in Floating_Number;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : Complex_Number;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : Floating_Number;
    cnt,nsuccess,trial : natural32 := 0;
    pivot : integer32 := s'first;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Multiple_Loose_Normal_Silent_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Multiple_Loose_Normal_Silent_Corrector(Norm,H,dH);

    procedure Correct ( s : in out Solu_Info_Array;
                        pivot : in out integer32;
                        dist_sols : in Floating_Number;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(s,pivot,dist_sols,c,fail);
       else Affine_Corrector(s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
    Copy(p.maxstep,step);
    Linear_Multiple_Initialize(s,p,t,old_t,prev_t,sa,old_sa,prev_sa);
    while not (At_End(t,target,p.dist_target,tol) and success) 
                            and (s(s'first).niter <= c.maxtot) loop
      Predictor(s,p,true,sa,prev_sa,t,prev_t,target,step,tol,dist_sols,trial);
      Correct(s,pivot,dist_sols,c,fail); Copy(s,sa);
      success := not fail;
      Linear_Multiple_Management(s,sa,old_sa,prev_sa,t,old_t,prev_t,p,step,
                                 pivot,nsuccess,trial,success);
      if step < p.minstep then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Correct(s,pivot,dist_sols,cp,fail);
    end;
    Clear(step); Clear(t); Clear(old_t); Clear(prev_t);
  end Linear_Multiple_Normal_Silent_Continue;

  procedure Linear_Multiple_Normal_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info_Array;
                target : in Complex_Number; tol,dist_sols : in Floating_Number;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : Complex_Number;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : Floating_Number;
    pivot : integer32 := s'first;
    cnt,nsuccess,trial : natural32 := 0;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Multiple_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Multiple_Loose_Normal_Reporting_Corrector(Norm,H,dH);

    procedure Correct ( file : in file_type; s : in out Solu_Info_Array;
                        pivot : in out integer32;
                        dist_sols : in Floating_Number;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(file,s,pivot,dist_sols,c,fail);
       else Affine_Corrector(file,s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
    Copy(p.maxstep,step);
    Linear_Multiple_Initialize(s,p,t,old_t,prev_t,sa,old_sa,prev_sa);
    for k in s'range loop                       -- write the start solutions
      sWrite(file,sa(k).all);
    end loop;
    while not (At_End(t,target,p.dist_target,tol) and success) 
                            and (s(s'first).niter <= c.maxtot) loop
      Predictor(s,p,true,sa,prev_sa,t,prev_t,target,step,tol,dist_sols,trial);
      pWrite(file,step,t);
      Correct(file,s,pivot,dist_sols,c,fail); Copy(s,sa);
      success := not fail;
      Linear_Multiple_Management(s,sa,old_sa,prev_sa,t,old_t,prev_t,p,step,
                                 pivot,nsuccess,trial,success);
      if step < p.minstep then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Correct(file,s,pivot,dist_sols,cp,fail);
    end;
    Clear(step); Clear(t); Clear(old_t); Clear(prev_t);
  end Linear_Multiple_Normal_Reporting_Continue;

  procedure Linear_Multiple_Conditioned_Silent_Continue
              ( s : in out Solu_Info_Array;
                target : in Complex_Number; tol,dist_sols : in Floating_Number;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : Complex_Number;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : Floating_Number;
    pivot : integer32 := s'first;
    cnt,nsuccess,trial : natural32 := 0;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Multiple_Loose_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Multiple_Loose_Conditioned_Silent_Corrector(Norm,H,dH);

    procedure Correct ( s : in out Solu_Info_Array;
                        pivot : in out integer32;
                        dist_sols : in Floating_Number;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(s,pivot,dist_sols,c,fail);
       else Affine_Corrector(s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
    Copy(p.maxstep,step);
    Linear_Multiple_Initialize(s,p,t,old_t,prev_t,sa,old_sa,prev_sa);
    while not (At_End(t,target,p.dist_target,tol) and success) 
                            and (s(s'first).niter <= c.maxtot) loop
      Predictor(s,p,true,sa,prev_sa,t,prev_t,target,step,tol,dist_sols,trial);
      Correct(s,pivot,dist_sols,c,fail); Copy(s,sa);
      success := not fail;
      Linear_Multiple_Management(s,sa,old_sa,prev_sa,t,old_t,prev_t,p,step,
                                 pivot,nsuccess,trial,success);
      if step < p.minstep then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Correct(s,pivot,dist_sols,cp,fail);
    end;
    Clear(step); Clear(t); Clear(old_t); Clear(prev_t);
  end Linear_Multiple_Conditioned_Silent_Continue;

  procedure Linear_Multiple_Conditioned_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info_Array;
                target : in Complex_Number; tol,dist_sols : in Floating_Number;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : Complex_Number;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : Floating_Number;
    pivot : integer32 := s'first;
    cnt,nsuccess,trial : natural32 := 0;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is 
      new Affine_Multiple_Loose_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is 
      new Projective_Multiple_Loose_Conditioned_Reporting_Corrector(Norm,H,dH);

    procedure Correct ( file : in file_type; s : in out Solu_Info_Array;
                        pivot : in out integer32;
                        dist_sols : in Floating_Number;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(file,s,pivot,dist_sols,c,fail);
       else Affine_Corrector(file,s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
    Copy(p.maxstep,step);
    Linear_Multiple_Initialize(s,p,t,old_t,prev_t,sa,old_sa,prev_sa);
    for k in s'range loop                        -- write start solutions
      sWrite(file,sa(k).all);
    end loop;
    while not (At_End(t,target,p.dist_target,tol) and success) 
                            and (s(s'first).niter <= c.maxtot) loop
      Predictor(s,p,true,sa,prev_sa,t,prev_t,target,step,tol,dist_sols,trial);
      pWrite(file,step,t);
      Correct(file,s,pivot,dist_sols,c,fail); Copy(s,sa);
      success := not fail;
      Linear_Multiple_Management(s,sa,old_sa,prev_sa,t,old_t,prev_t,p,step,
                                 pivot,nsuccess,trial,success);
      if step < p.minstep then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Correct(file,s,pivot,dist_sols,cp,fail);
    end;
    Clear(step); Clear(t); Clear(old_t); Clear(prev_t);
  end Linear_Multiple_Conditioned_Reporting_Continue;

-- CIRCULAR PATH FOLLOWING FOR ONE PATH :

  procedure Circular_Single_Normal_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in Complex_Number; tol,epslop : in Floating_Number;
                wc : out natural32; max_wc : in natural32;
                sum,all_sum : out Vector;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t : Complex_Number;
    t0_min_target : Complex_Number := s.sol.t - target;
    theta,old_theta : Floating_Number := Create(integer(0));
    two : Floating_Number := Create(integer(2));
    twopi : constant Floating_Number := two*Create(PI);
    step : Floating_Number := twopi*p.maxstep;
    old_solution,prev_solution,start_solution : Vector(s.sol.v'range);
    w_c,nsuccess : natural32 := 0;
    success : boolean := true;
    stop : boolean := false;
    w_sum,w_all_sum : Vector(s.sol.v'range);
    n_sum,n_all_sum : natural32 := 0;
    acc : Floating_Number;

    procedure T_C_Predictor is new Tangent_Circular_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);

  begin
    old_t := s.sol.t; old_solution := s.sol.v;        -- INITIALIZATION
    for i in s.sol.v'range loop
      w_sum(i) := s.sol.v(i)/two;
      Copy(w_sum(i),w_all_sum(i));
    end loop;
    Clear(two);
    start_solution := s.sol.v;
    if p.predictor_type = 0 
     then prev_t := s.sol.t; prev_solution := s.sol.v;
    end if;
    sWrite(file,s.sol.all);                -- write the start solution
    while (s.niter <= c.maxtot) loop 
      case p.predictor_type is
        when 0 => Secant_Circular_Predictor
                     (s.sol.v,prev_solution,s.sol.t,theta,
                      prev_t,t0_min_target,target,step,tol);
        when 2 => T_C_Predictor(s.sol.v,s.sol.t,theta, t0_min_target,target,
                                step,tol);
                  s.nsyst := s.nsyst + 1;
        when others => null; -- these choices make no sense !!!
      end case;
      pWrite(file,step,s.sol.t,s.sol.all);
      if proj
       then Projective_Corrector(file,s,c);
       else Affine_Corrector(file,s,c);
      end if;
      sWrite(file,s.sol.all);
      Circular_Management
           (s,p,c,old_t,prev_t,start_solution,old_solution,prev_solution,
            w_sum,w_all_sum,twopi,epslop,tol,theta,old_theta,
            step,nsuccess,n_sum,n_all_sum,w_c,max_wc,stop,success);
      exit when stop;
      if step < p.minstep then wc := 0; return; end if;
    end loop;
    wc := w_c;
    if n_all_sum /= 0 then
      acc := Create(n_all_sum);
      for i in w_all_sum'range loop
        all_sum(i) := w_all_sum(i)/acc;
      end loop;
      Clear(acc);
    end if;
    if n_sum /= 0 then
      acc := Create(n_sum);
      for i in w_sum'range loop
        sum(i) := w_sum(i)/acc;
      end loop;
      Clear(acc);
    elsif n_all_sum /= 0 then
      acc := Create(n_all_sum);
      for i in w_all_sum'range loop
        all_sum(i) := w_all_sum(i)/acc;
      end loop;
      Clear(acc);
    end if;
  end Circular_Single_Normal_Reporting_Continue;

  procedure Circular_Single_Conditioned_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in Complex_Number; tol,epslop : in Floating_Number;
                wc : out natural32; max_wc : in natural32;
                sum,all_sum : out Vector;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t : Complex_Number;
    theta,old_theta : Floating_Number := Create(integer(0));
    two : Floating_Number := Create(integer(2));
    twopi : Floating_Number := two*Create(PI);
    step : Floating_Number := twopi*p.maxstep;
    t0_min_target : Complex_Number := s.sol.t - target;
    old_solution,prev_solution,start_solution : Vector(s.sol.v'range);
    w_c,nsuccess : natural32 := 0;
    success : boolean := true;
    stop : boolean := false;
    w_sum,w_all_sum : Vector(s.sol.v'range);
    n_sum,n_all_sum : natural32 := 0;
    acc : Floating_Number;

    procedure T_C_Predictor is new Tangent_Circular_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);

  begin
    old_t := s.sol.t; old_solution := s.sol.v;            -- INITIALIZATION
    for i in s.sol.v'range loop
      w_sum(i) := s.sol.v(i)/two;
      Copy(w_sum(i),w_all_sum(i));
    end loop;
    Clear(two);
    start_solution := s.sol.v;
    if p.predictor_type = 0
     then prev_t := s.sol.t; prev_solution := old_solution;
    end if;
    sWrite(file,s.sol.all);              -- writing the start solution
    while s.niter <= c.maxtot loop
      case p.predictor_type is
        when 0 => Secant_Circular_Predictor
                    (s.sol.v,prev_solution,s.sol.t,theta,
                     prev_t,t0_min_target,target,step,tol);
        when 2 => T_C_Predictor(s.sol.v,s.sol.t,theta,t0_min_target,
                                target,step,tol);
                  s.nsyst := s.nsyst + 1;
        when others => null; -- these choices make no sense !!!
      end case;
      pWrite(file,step,s.sol.t,s.sol.all);
      if proj
       then Projective_Corrector(file,s,c);
       else Affine_Corrector(file,s,c);
      end if;
      sWrite(file,s.sol.all);
      Circular_Management
          (s,p,c,old_t,prev_t,start_solution,old_solution,prev_solution,
           w_sum,w_all_sum,twopi,epslop,tol,theta,old_theta,
           step,nsuccess,n_sum,n_all_sum,w_c,max_wc,stop,success);
      exit when stop;
      if step < p.minstep then wc := 0; return; end if;
    end loop;
    wc := w_c;
    if n_all_sum /= 0 then
      acc := Create(n_all_sum);
      for i in w_all_sum'range loop
        all_sum(i) := w_all_sum(i)/acc;
      end loop;
      Clear(acc);
    end if;
    if n_sum /= 0 then
      acc := Create(n_sum);
      for i in w_sum'range loop
        sum(i) := w_sum(i)/acc;
      end loop;
      Clear(acc);
    elsif n_all_sum /= 0 then
      acc := Create(n_all_sum);
      for i in w_all_sum'range loop
        sum(i) := w_all_sum(i)/acc;
      end loop;
      Clear(acc);
    end if;
  end Circular_Single_Conditioned_Reporting_Continue;

end Multprec_Path_Trackers;
