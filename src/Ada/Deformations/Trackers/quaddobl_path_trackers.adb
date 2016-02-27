with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Mathematical_Functions;
with Quad_Double_Vectors_io;             use Quad_Double_Vectors_io;
with Quad_Double_VecVecs;                use Quad_Double_VecVecs;
with QuadDobl_Complex_Equality_Tests;    use QuadDobl_Complex_Equality_Tests;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;
with QuadDobl_Predictors;                use QuadDobl_Predictors;
with QuadDobl_Correctors;                use QuadDobl_Correctors;
with QuadDobl_Orthogonal_Correctors;     use QuadDobl_Orthogonal_Correctors;
with QuadDobl_Dispatch_Predictors;       use QuadDobl_Dispatch_Predictors;
with Process_io;                         use Process_io;
with Directions_of_QuadDobl_Paths;       use Directions_of_QuadDobl_Paths;
with Standard_Data_on_Path;
with Quaddobl_Data_on_Path;              use Quaddobl_Data_on_Path;

package body QuadDobl_Path_Trackers is

-- LINEAR PATH FOLLOWING FOR ONE PATH :

  procedure Linear_Single_Normal_Silent_Continue
              ( s : in out Solu_Info; target : in Complex_Number;
                tol : in double_float; proj : in boolean;
                p : in Pred_Pars; c : in Corr_Pars; nbq : in integer32 := 0;
                f : access procedure ( s : in Solu_Info ) := null ) is
 
    old_t,prev_t,prev_t2,prev_t1,prev_t0 : Complex_Number;
    old_sol,prev_sol,prev_sol2,prev_sol1,prev_sol0 : Vector(s.sol.v'range);
    old_v,prev_v,vv : Vector(s.sol.v'range);
    step : double_float := p.maxstep;
    nsuccess,trial : natural32 := 0;
    success : boolean := true;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);
      procedure QRLS_Corrector is
        new Silent_QRLS_Corrector(Norm,H,dH);

    begin
      if proj then
        Projective_Corrector(s,c);
      elsif nbq = 0 then
        Affine_Corrector(s,c);
      else
        QRLS_Corrector(nbq,s,c);
      end if;
    end Corrector;

  begin
    Linear_Single_Initialize(s,p,old_t,prev_t,prev_v,old_sol,prev_sol);
    if p.predictor_type > 6 then
      prev_t0 := prev_t; prev_sol0 := prev_sol;
      prev_t1 := prev_t; prev_sol1 := prev_sol;
      prev_t2 := prev_t; prev_sol2 := prev_sol;
    end if;
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success) 
                                           and (s.niter <= c.maxtot) loop
      if p.predictor_type < 7 then
        Predictor(s,p,success,prev_sol,prev_v,vv,prev_t,target,step,tol,trial);
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
      if (f /= null) then f(s); end if;
      if p.predictor_type < 7 then
        Linear_Single_Management(s,p,c,old_t,prev_t,old_sol,prev_sol,
                                 old_v,prev_v,vv,step,nsuccess,trial,success);
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
      cp : constant Corr_Pars
         := Standard_Data_on_Path.End_Game_Corrector_Parameters
              (c,p.dist_target,tol);
    begin
      Corrector(s,cp);
      if (f /= null) then f(s); end if;
    end;
  end Linear_Single_Normal_Silent_Continue;

  procedure Linear_Single_Normal_Reporting_Continue
              ( file : in file_type;
                s : in out Solu_Info; target : in Complex_Number;
                tol : in double_float; proj : in boolean;
                p : in Pred_Pars; c : in Corr_Pars; nbq : in integer32 := 0;
                f : access procedure ( s : in Solu_Info ) := null ) is
 
    old_t,prev_t,prev_t2,prev_t1,prev_t0 : Complex_Number;
    old_sol,prev_sol,prev_sol2,prev_sol1,prev_sol0 : Vector(s.sol.v'range);
    old_v,prev_v,vv : Vector(s.sol.v'range);
    step : double_float := p.maxstep;
    nsuccess,trial : natural32 := 0;
    success : boolean := true;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);
      procedure QRLS_Corrector is
        new Reporting_QRLS_Corrector(Norm,H,dH);

    begin
      if proj then
        Projective_Corrector(file,s,c);
      elsif nbq = 0 then
        Affine_Corrector(file,s,c);
      else
        QRLS_Corrector(file,nbq,s,c);
      end if;
    end Corrector;

  begin
    Linear_Single_Initialize(s,p,old_t,prev_t,prev_v,old_sol,prev_sol);
    if p.predictor_type > 6 then
      prev_t0 := prev_t; prev_sol0 := prev_sol;
      prev_t1 := prev_t; prev_sol1 := prev_sol;
      prev_t2 := prev_t; prev_sol2 := prev_sol;
    end if;
    sWrite(file,s.sol.all);
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success) 
                                           and (s.niter <= c.maxtot) loop
      if p.predictor_type < 7 then
        Predictor(s,p,success,prev_sol,prev_v,vv,prev_t,target,step,tol,trial);
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
      if (f /= null) then f(s); end if;
      sWrite(file,s.sol.all);
      if p.predictor_type < 7 then
        Linear_Single_Management(s,p,c,old_t,prev_t,old_sol,prev_sol,
                                 old_v,prev_v,vv,step,nsuccess,trial,success);
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
      cp : constant Corr_Pars
         := Standard_Data_on_Path.End_Game_Corrector_Parameters
              (c,p.dist_target,tol);
    begin
      Corrector(s,cp);
      if (f /= null) then f(s); end if;
      sWrite(file,s.sol.all);
    end;
  --exception
  --  when others =>
  --    put_line("exception raised in Linear Single Normal Reporting continue");
  --    raise;
  end Linear_Single_Normal_Reporting_Continue;

  procedure Linear_Single_Conditioned_Silent_Continue
              ( s : in out Solu_Info; target : in Complex_Number;
                tol : in double_float; proj : in boolean;
                rtoric : in integer32; w : in out integer32;
                v : in out Quad_Double_Vectors.Link_to_Vector;
                errorv : in out quad_double;
                p : in Pred_Pars; c : in Corr_Pars; nbq : in integer32 := 0;
                f : access procedure ( s : in Solu_Info ) := null ) is

    old_t,prev_t,prev_t2,prev_t1,prev_t0 : Complex_Number;
    old_sol,prev_sol,prev_sol2,prev_sol1,prev_sol0 : Vector(s.sol.v'range);
    old_v,prev_v,vv : Vector(s.sol.v'range);
    step : double_float := p.maxstep;
    nsuccess,trial : natural32 := 0;
    success : boolean := true;

    dls,stp : double_float := 0.0;
    tolsing : constant double_float
            := Continuation_Parameters.tol_endg_inverse_condition;
    r : integer32 := 0;
    er : integer32 := 0;
    zero : constant quad_double := create(0.0);
    dt,ds,logs : Quad_Double_Vectors.Vector(0..rtoric)
               := (0..rtoric => zero);
    logx : VecVec(0..rtoric);
    wvl0,wvl1,wvl2 : VecVec(1..rtoric);
    errv : Quad_Double_Vectors.Vector(0..rtoric) := (0..rtoric => zero);
    m : constant integer32 := 1;
    thresm : constant natural32 := p.success_steps;
    estm : integer32 := m;
    fcnt,cntm : natural32 := 0;
    defer : natural32 := thresm;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);
      procedure SVD_Corrector is
        new Silent_SVD_Corrector(Norm,H,dH);

    begin
      if proj then
        Projective_Corrector(s,c);
      elsif nbq = 0 then
        Affine_Corrector(s,c);
      else
        SVD_Corrector(nbq,s,c);
      end if;
    end Corrector;

  begin
    Linear_Single_Initialize(s,p,old_t,prev_t,prev_v,old_sol,prev_sol);
    if p.predictor_type > 6 then
      prev_t0 := prev_t; prev_sol0 := prev_sol;
      prev_t1 := prev_t; prev_sol1 := prev_sol;
      prev_t2 := prev_t; prev_sol2 := prev_sol;
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
              w := s.sol.m;
         else er := -2;
        end if;
      end if;
      if p.predictor_type < 7 then
        Predictor(s,p,success,prev_sol,prev_v,vv,prev_t,target,step,tol,trial);
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
      if (f /= null) then f(s); end if;
      if p.predictor_type < 7 then
        Linear_Single_Management(s,p,c,old_t,prev_t,old_sol,prev_sol,
                                 old_v,prev_v,vv,step,nsuccess,trial,success);
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
      cp : constant Corr_Pars 
         := Standard_Data_on_Path.End_Game_Corrector_Parameters
              (c,p.dist_target,tol);
    begin
      Corrector(s,cp);
      if (f /= null) then f(s); end if;
    end;
  end Linear_Single_Conditioned_Silent_Continue;

  procedure Linear_Single_Conditioned_Reporting_Continue
              ( file : in file_type;
                s : in out Solu_Info; target : in Complex_Number;
                tol : in double_float; proj : in boolean; 
                rtoric : in integer32; w : in out integer32;
                v : in out Quad_Double_Vectors.Link_to_Vector;
                errorv : in out quad_double;
                p : in Pred_Pars; c : in Corr_Pars; nbq : in integer32 := 0;
                f : access procedure ( s : in Solu_Info ) := null ) is

    old_t,prev_t,prev_t2,prev_t1,prev_t0 : Complex_Number;
    step : double_float := p.maxstep;
    nsuccess,trial : natural32 := 0;
    old_sol,prev_sol,prev_sol2,prev_sol1,prev_sol0 : Vector(s.sol.v'range);
    old_v,prev_v,vv : Vector(s.sol.v'range);
    success : boolean := true;
    tolsing : constant double_float
            := Continuation_Parameters.tol_endg_inverse_condition;
    r : integer32 := 0;
    er : integer32 := 0;
    zero : constant quad_double := create(0.0);
    dt,ds,logs : Quad_Double_Vectors.Vector(0..rtoric)
               := (0..rtoric => zero);
    logx : VecVec(0..rtoric);
    wvl0,wvl1,wvl2 : VecVec(1..rtoric);
    errv : Quad_Double_Vectors.Vector(0..rtoric) := (0..rtoric => zero);
    m : constant integer32 := 1;
    thresm : constant natural32 := p.success_steps;
    estm : integer32 := m;
    fcnt,cntm : natural32 := 0;
    defer : natural32 := thresm;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
      procedure SVD_Corrector is
        new Reporting_SVD_Corrector(Norm,H,dH);

    begin
      if proj then
        Projective_Corrector(file,s,c);
      elsif nbq = 0 then
        Affine_Corrector(file,s,c);
      else
        SVD_Corrector(file,nbq,s,c);
      end if;
    exception
      when others => --put("exception raised in ");
    --    if proj then put_line(" Projective Corrector");
    --            else put_line(" Affine Corrector"); end if;
        return;
    end Corrector;

  begin
    Linear_Single_Initialize(s,p,old_t,prev_t,prev_v,old_sol,prev_sol);
    if p.predictor_type > 6 then
      prev_t0 := prev_t; prev_sol0 := prev_sol;
      prev_t1 := prev_t; prev_sol1 := prev_sol;
      prev_t2 := prev_t; prev_sol2 := prev_sol;
    end if;
    sWrite(file,s.sol.all);          -- writing the start solution
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
              w := s.sol.m;
         else er := -2;
        end if;
      end if;
      if p.predictor_type < 7 then
        Predictor(s,p,success,prev_sol,prev_v,vv,prev_t,target,step,tol,trial);
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
      if (f /= null) then f(s); end if;
      sWrite(file,s.sol.all);
      if p.predictor_type < 7 then
        Linear_Single_Management(s,p,c,old_t,prev_t,old_sol,prev_sol,
                                 old_v,prev_v,vv,step,nsuccess,trial,success);
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
      cp : constant Corr_Pars
         := Standard_Data_on_Path.End_Game_Corrector_Parameters
              (c,p.dist_target,tol);
    begin
      Corrector(s,cp);
      if (f /= null) then f(s); end if;
    end;
  end Linear_Single_Conditioned_Reporting_Continue;

-- LINEAR PATH FOLLOWING FOR A NUMBER OF PATHS :

  procedure Linear_Multiple_Normal_Silent_Continue
              ( s : in out Solu_Info_Array;
                target : in Complex_Number; tol,dist_sols : in double_float;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : Complex_Number;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : double_float := p.maxstep;
    nsuccess,trial : natural32 := 0;
    pivot : integer32 := s'first;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Multiple_Loose_Normal_Silent_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Multiple_Loose_Normal_Silent_Corrector(Norm,H,dH);

    procedure Correct ( s : in out Solu_Info_Array;
                        pivot : in out integer32; dist_sols : in double_float;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(s,pivot,dist_sols,c,fail);
       else Affine_Corrector(s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
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
      cp : constant Corr_Pars
         := Standard_Data_on_Path.End_Game_Corrector_Parameters
              (c,p.dist_target,tol);
    begin
      Correct(s,pivot,dist_sols,cp,fail);
    end;
  end Linear_Multiple_Normal_Silent_Continue;

  procedure Linear_Multiple_Normal_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info_Array;
                target : in Complex_Number; tol,dist_sols : in double_float;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : Complex_Number;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : double_float := p.maxstep;
    pivot : integer32 := s'first;
    nsuccess,trial : natural32 := 0;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Multiple_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Multiple_Loose_Normal_Reporting_Corrector(Norm,H,dH);

    procedure Correct ( file : in file_type; s : in out Solu_Info_Array;
                        pivot : in out integer32; dist_sols : in double_float;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(file,s,pivot,dist_sols,c,fail);
       else Affine_Corrector(file,s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
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
      cp : constant Corr_Pars
         := Standard_Data_on_Path.End_Game_Corrector_Parameters
              (c,p.dist_target,tol);
    begin
      Correct(file,s,pivot,dist_sols,cp,fail);
    end;
  end Linear_Multiple_Normal_Reporting_Continue;

  procedure Linear_Multiple_Conditioned_Silent_Continue
              ( s : in out Solu_Info_Array;
                target : in Complex_Number; tol,dist_sols : in double_float;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : Complex_Number;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : double_float := p.maxstep;
    pivot : integer32 := s'first;
    nsuccess,trial : natural32 := 0;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Multiple_Loose_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Multiple_Loose_Conditioned_Silent_Corrector(Norm,H,dH);

    procedure Correct ( s : in out Solu_Info_Array;
                        pivot : in out integer32; dist_sols : in double_float;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(s,pivot,dist_sols,c,fail);
       else Affine_Corrector(s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
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
      cp : constant Corr_Pars 
         := Standard_Data_on_Path.End_Game_Corrector_Parameters
              (c,p.dist_target,tol);
    begin
      Correct(s,pivot,dist_sols,cp,fail);
    end;
  end Linear_Multiple_Conditioned_Silent_Continue;

  procedure Linear_Multiple_Conditioned_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info_Array;
                target : in Complex_Number; tol,dist_sols : in double_float;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : Complex_Number;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : double_float := p.maxstep;
    pivot : integer32 := s'first;
    nsuccess,trial : natural32 := 0;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is 
      new Affine_Multiple_Loose_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is 
      new Projective_Multiple_Loose_Conditioned_Reporting_Corrector(Norm,H,dH);

    procedure Correct ( file : in file_type; s : in out Solu_Info_Array;
                        pivot : in out integer32; dist_sols : in double_float;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(file,s,pivot,dist_sols,c,fail);
       else Affine_Corrector(file,s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
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
      cp : constant Corr_Pars
         := Standard_Data_on_Path.End_Game_Corrector_Parameters
              (c,p.dist_target,tol);
    begin
      Correct(file,s,pivot,dist_sols,cp,fail);
    end;
  end Linear_Multiple_Conditioned_Reporting_Continue;

-- CIRCULAR PATH FOLLOWING FOR ONE PATH :

  procedure Circular_Single_Normal_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in Complex_Number; tol,epslop : in double_float;
                wc : out natural32; max_wc : in natural32;
                sum,all_sum : out Vector;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t : Complex_Number;
    t0_min_target : constant Complex_Number := s.sol.t - target;
    theta,old_theta : double_float := 0.0;
    twopi : constant double_float := 2.0*Standard_Mathematical_Functions.PI;
    step : double_float := twopi*p.maxstep;
    old_solution,prev_solution,start_solution : Vector(s.sol.v'range);
    w_c,nsuccess : natural32 := 0;
    success : boolean := true;
    stop : boolean := false;
    a_half : constant quad_double := create(0.5);
    w_sum,w_all_sum : Vector(s.sol.v'range) := Create(a_half)*s.sol.v;
    n_sum,n_all_sum : natural32 := 0;

    procedure T_C_Predictor is new Tangent_Circular_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);

  begin
    wc := 0;
    old_t := s.sol.t; old_solution := s.sol.v;        -- INITIALIZATION
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
      if step < p.minstep then return; end if;
    end loop;
    wc := w_c;
    if n_all_sum /= 0
     then all_sum := w_all_sum*Create(1.0/double_float(n_all_sum));
    end if;
    if n_sum /= 0
     then sum := w_sum*Create(1.0/double_float(n_sum));
     elsif n_all_sum /= 0
         then all_sum := w_all_sum*Create(1.0/double_float(n_all_sum));
    end if;
  end Circular_Single_Normal_Reporting_Continue;

  procedure Circular_Single_Conditioned_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in Complex_Number; tol,epslop : in double_float;
                wc : out natural32; max_wc : in natural32;
                sum,all_sum : out Vector;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t : Complex_Number;
    theta,old_theta : double_float := 0.0;
    twopi : constant double_float := 2.0*Standard_Mathematical_Functions.PI;
    step : double_float := twopi*p.maxstep;
    t0_min_target : constant Complex_Number := s.sol.t - target;
    old_solution,prev_solution,start_solution : Vector(s.sol.v'range);
    w_c,nsuccess : natural32 := 0;
    success : boolean := true;
    stop : boolean := false;
    a_half : constant quad_double := create(0.5);
    w_sum,w_all_sum : Vector(s.sol.v'range) := Create(a_half)*s.sol.v;
    n_sum,n_all_sum : natural32 := 0;

    procedure T_C_Predictor is new Tangent_Circular_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);

  begin
    wc := 0;
    old_t := s.sol.t; old_solution := s.sol.v;            -- INITIALIZATION
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
      if step < p.minstep then return; end if;
    end loop;
    wc := w_c;
    if n_all_sum /= 0
     then all_sum := w_all_sum*Create(1.0/double_float(n_all_sum));
    end if;
    if n_sum /= 0
     then sum := w_sum*Create(1.0/double_float(n_sum));
     elsif n_all_sum /= 0
         then sum := w_all_sum*Create(1.0/double_float(n_all_sum));
    end if;
  end Circular_Single_Conditioned_Reporting_Continue;

end QuadDobl_Path_Trackers;
