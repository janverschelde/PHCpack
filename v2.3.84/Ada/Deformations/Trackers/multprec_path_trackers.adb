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
--with Multprec_Directions_of_Paths;       use Multprec_Directions_of_Paths;

package body Multprec_Path_Trackers is

-- AUXILIAIRIES :

  function At_End ( t,target : Complex_Number;
                    distance,tol : Floating_Number ) return boolean is

  -- DESCRIPTION :
  --   Decides whether at end of continuation stage.

    res : boolean;
    diff : Complex_Number;
    absdiff : Floating_Number;

  begin
    if distance < tol then
      res := Equal(t,target,tol);
    else
      diff := t-target;
      absdiff := AbsVal(diff);
      if absdiff < distance
       then res := true;
       else res := false;
      end if;
      Clear(diff); Clear(absdiff);
    end if;
    return res;
  end At_End;

  function Stop ( p : Pred_Pars; t,target : Complex_Number;
                  step : Floating_Number ) return boolean is

  -- DESCRIPTION :
  --   Returns true if either the step size is smaller than p.minstep, or
  --   or alternatively, in case of geometric predictor, if the distance to
  --   the target has become smaller than p.minstep.

    res : boolean;
    diff : Complex_Number;
    absdiff : Floating_Number;

  begin
    if step < p.minstep
     then res := true;
     else if (p.predictor_type = 2) or (p.predictor_type = 5)
           then diff := t-target;
                absdiff := AbsVal(diff);
                res := (absdiff < p.minstep);
                Clear(diff); Clear(absdiff);
           else res := false;
          end if;
    end if;
    return res;
  end Stop;

  function Is_Multiple ( a,b,tol : Floating_Number ) return natural32 is

  -- DESCRIPTION :
  --   Returns a/b if a is a multiple of b, returns 0 in the other case.

    res : natural32;
    dfq : double_float;
    fq,absb,diff,absdiff : Floating_Number;
    iq : integer32;

  begin
    absb := AbsVal(b);
    if absb < tol then
      res := 0;
    else
      fq := a/b;
      dfq := Round(fq);
      iq := integer32(dfq);
      diff := Create(iq);
      Sub(diff,fq);
      absdiff := AbsVal(diff);
      if absdiff < tol
       then res := natural32(iq);
       else res := 0;
      end if;
      Clear(diff); Clear(absdiff);
    end if;
    Clear(absb);
    return res;
  end Is_Multiple;

  procedure Linear_Step_Control 
              ( success : in boolean; p : in Pred_Pars;
                step : in out Floating_Number;
                nsuccess : in out natural32; trial : in natural32 ) is

  -- DESCRIPTION :
  --   Control of step size for linear path following.
  --   With geometric prediction, the ratio (=step) will be enlarged
  --   when not success.  In order not to break the sequence, the ratio
  --   is not reduced when success.

    acc : Floating_Number;

  begin
    if (p.predictor_type = 2) or (p.predictor_type = 5) then
      if success then
        nsuccess := nsuccess + 1;
      else
        nsuccess := 0;
        acc := 1.0 - step;
        if p.expfac < 1.0
         then Mul(acc,p.expfac); Add(step,acc);
         else Div(acc,p.expfac); Add(step,acc);
        end if;
        Clear(acc);
        if step > 1.0
         then step := p.minstep/2.0;  -- that ends the game
        end if;
      end if;
    else
      if success then
        nsuccess := nsuccess + 1;
        if nsuccess > p.success_steps then
          Mul(step,p.expfac);
          if step > p.maxstep
           then Copy(p.maxstep,step);
          end if;
        end if;
      else
        nsuccess := 0;
        if trial mod 3 = 0
         then Mul(step,p.redfac);
        end if;
      end if;
    end if;
  end Linear_Step_Control;

  procedure Circular_Step_Control
             ( success : in boolean; p : in Pred_Pars;
               twopi : in Floating_Number;
               step : in out Floating_Number; nsuccess : in out natural32 ) is

  -- DESCRIPTION :
  --   Control of step size for circular path following, note that the
  --   step size should be some multiple of pi.

    maxstep : constant Floating_Number := p.maxstep*twopi;

  begin
    if success
     then nsuccess := nsuccess + 1;
          if nsuccess > p.success_steps
           then Mul(step,2.0);             -- expansion factor = 2
                if step > maxstep
                 then Copy(maxstep,step);
                end if;
          end if;
     else nsuccess := 0;
          Div(step,2.0);                   -- reduction factor = 1/2
    end if;
  end Circular_Step_Control;

  procedure Set_Corrector_Parameters
              ( c : in out Corr_Pars; eps : Floating_Number ) is

  -- DESCRIPTION :
  --   All eps* parameters in c are set to eps.

  begin
    Copy(eps,c.epsrx); Copy(eps,c.epsax);
    Copy(eps,c.epsrf); Copy(eps,c.epsaf);
  end Set_Corrector_Parameters;

  function End_Game_Corrector_Parameters
             ( current : Corr_Pars; distance,tol : Floating_Number ) 
             return Corr_Pars is

  -- DESCRIPTION :
  --   Returns corrector parameters for the end game of the first or the 
  --   second continuation stage, depending on the distance from the target.

    res : Corr_Pars;

  begin
    if distance < tol                 -- correct further to detect clustering
     then Copy(current.epsrx,res.epsrx);
          Copy(current.epsax,res.epsax);
          Copy(current.epsrf,res.epsrf);
          Copy(current.epsaf,res.epsaf);
          res.maxit := current.maxit;
          res.maxtot := current.maxtot;
          Set_Corrector_Parameters(res,tol);
     else res := Convert(Continuation_Parameters.Create_End_Game);
         -- or to move smoothly to end game
    end if;
    return res;
  end End_Game_Corrector_Parameters;

-- MANAGEMENT OF DATA DURING PATH FOLLOWING :

  procedure Linear_Single_Initialize 
                      ( s : in Solu_Info; p : in Pred_Pars;
                        old_t,prev_t : out Complex_Number;
                        prev_v,old_solution,prev_solution : out Vector ) is

  -- DESCRIPTION :
  --   Initialization for linear path following of one path.

  -- ON ENTRY :
  --   s                solution at beginning of path;
  --   p                predictor parameters.

  -- ON RETURN :
  --   old_t            back up value for continuation parameter;
  --   prev_t           previous value of continuation parameter;
  --   old_solution     back up value for solution vector;
  --   prev_solution    previous value of solution vector;

  begin
    Copy(s.sol.t,old_t);
    Copy(s.sol.v,old_solution);
    if p.predictor_type <= 2                    -- for all secant predictors
     then Copy(s.sol.t,prev_t);
          Copy(s.sol.v,prev_solution);
    end if;
    prev_v := (prev_v'range => Create(integer(0)));
  end Linear_Single_Initialize;

  procedure Linear_Single_Management
                 ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                   old_t,prev_t : in out Complex_Number;
                   old_solution,prev_solution,old_v,prev_v,vv : in out Vector;
                   step : in out Floating_Number;
                   nsuccess,trial : in out natural32;
                   success : in out boolean ) is

  -- DESCRIPTION :
  --   Management of data after correction during linear path following.

  -- PARAMETERS :
  --   s                current solution;
  --   p                predictor parameters;
  --   c                corrector parameters;
  --   old_t            back up value for continuation parameter;
  --   prev_t           previous value of continuation parameter;
  --   old_solution     back up value for solution vector;
  --   prev_solution    previous value of solution vector;
  --   old_v            back up value for tangent direction;
  --   prev_v           previous value for tangent direction;
  --   vv               current tangent direction;
  --   step             current step size;
  --   nsuccess         number of consecutive successes;
  --   trial            number of trials after failure;
  --   success          successful correction step.

  begin
    success := (((s.resa < c.epsaf) or else (s.cora < c.epsax))
        or else ((s.resr < c.epsrf) or else (s.corr < c.epsrx)));
    if ((p.predictor_type <= 2) and then success)         -- secant predictors
     then Copy(old_t,prev_t);
          Copy(old_solution,prev_solution);
    end if;
    if ((p.predictor_type = 6) and then success)          -- Hermite predictor
     then Copy(old_t,prev_t);
          Copy(old_solution,prev_solution);
          Copy(old_v,prev_v);
    end if;
    if ((p.predictor_type = 1) or (p.predictor_type = 3)) -- complex predictors
     then if success
           then trial := 0;
           else trial := trial + 1;
          end if;
    end if;
    s.nstep := s.nstep + 1;
    if not success
     then s.nfail := s.nfail + 1;
    end if;
    Linear_Step_Control(success,p,step,nsuccess,trial);
    if step < p.minstep
     then return;
    end if;
    if not success
     then Copy(old_t,s.sol.t);
          Copy(old_solution,s.sol.v);
     else Copy(s.sol.t,old_t);
          Copy(s.sol.v,old_solution);
          if p.predictor_type = 6
           then Copy(vv,old_v);
          end if;
    end if;
  end Linear_Single_Management;

  procedure Linear_Multiple_Initialize 
                 ( s : in Solu_Info_Array; p : in Pred_Pars;
                   t,old_t,prev_t : out Complex_Number;
                   sa,old_sa,prev_sa : in out Solution_Array ) is

  -- DECRIPTION :
  --   Initialization for linear path following for more than one path.

  begin
    Copy(s(s'first).sol.t,t);
    Copy(s(s'first).sol.t,old_t);
    Copy(s,sa); Copy(sa,old_sa);
    case p.predictor_type is
      when 0 | 1 | 2 => Copy(s(s'first).sol.t,prev_t);
                        Copy(sa,prev_sa);
      when others => null;
    end case;
  end Linear_Multiple_Initialize;

  procedure Linear_Multiple_Management
                  ( s : in out Solu_Info_array;
                    sa,old_sa,prev_sa : in out Solution_Array;
                    t,old_t,prev_t : in out Complex_Number; p : in Pred_Pars; 
                    step : in out Floating_Number; pivot : in integer32; 
                    nsuccess,trial : in out natural32;
                    success : in boolean ) is

  -- DESCRIPTION :
  --   Management of data after correction during linear path following.

  -- PARAMETERS :
  --   s            current solutions with information statistics;
  --   sa           current solutions;
  --   old_sa       back up value for solutions;
  --   prev_sa      previous solutions;
  --   t            current value of continuation parameter;
  --   old_t        back up value for continuation parameter;
  --   prev_t       previous value of continuation parameter;
  --   p            predictor parameters;
  --   step         current step size;
  --   pivot        solution where correction is started;
  --   nsuccess     number of consecutive successes;
  --   trial        number of trials after failure;
  --   success      successful correction step.

  begin
    if ((p.predictor_type <= 2) and then success)      -- for secant predictors
     then Copy(old_t,prev_t);
          Copy(old_sa,prev_sa);
    end if;
    if ((p.predictor_type = 1) or (p.predictor_type = 3)) -- complex predictors
     then if success
           then trial := 0;
           else trial := trial + 1;
          end if;
    end if;
    for k in s'range loop
      s(k).nstep := s(k).nstep + 1;
    end loop;
    if not success
     then s(pivot).nfail := s(pivot).nfail + 1;
    end if;
    Linear_Step_Control(success,p,step,nsuccess,trial);
    if step < p.minstep then return; end if;
    if success
     then Copy(sa,old_sa); Copy(t,old_t);
     else Copy(old_sa,sa); Copy(old_t,t);
    end if;
  end Linear_Multiple_Management;

  procedure Circular_Management
                   ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                     old_t,prev_t : in out Complex_Number;
                     start_solution : in Vector;
                     old_solution,prev_solution,w_sum,w_all_sum : in out Vector;
                     twopi,epslop,tol : in Floating_Number;
                     theta,old_theta,step : in out Floating_Number;
                     nsuccess,n_sum,n_all_sum,w_c : in out natural32;
                     max_wc : in natural32; stop,success : in out boolean ) is

  -- DESCRIPTION :
  --   Management of circular path following.

  -- PARAMETERS :
  --   s                current solution;
  --   p                predictor parameters;
  --   c                corrector parameters;
  --   old_t            back up value for continuation parameter;
  --   prev_t           previous value of continuation parameter;
  --   start_solution   solution vector at start of continuation;
  --   old_solution     back up value for solution vector;
  --   prev_solution    previous value of solution vector;
  --   w_sum            sum of cycle;
  --   w_all_sum        sum of all cycles;
  --   twopi            two times PI;
  --   epslop           tolerance to decide whether two vectors are equal;
  --   theta            current value of theta;
  --   old_theta        back up value for theta;
  --   step             current step size;
  --   nsuccess         number of consecutive successes;
  --   n_sum            number of cycles;
  --   n_all_sum        total number of cycles;
  --   w_c              winding number;
  --   max_wc           upper bound on winding number;
  --   stop             true when winding number has been computed;
  --   success          successful correction step.

    tmp : natural32;
    two,acc : Complex_Number;

  begin
    success := (((s.resa < c.epsaf) or else (s.cora < c.epsax))
        or else ((s.resr < c.epsrf) or else (s.corr < c.epsrx)));
    if p.predictor_type = 0 and then success
     then Copy(old_t,prev_t);
          Copy(old_solution,prev_solution);
    end if;
    s.nstep := s.nstep + 1;
    if not success then s.nfail := s.nfail + 1; end if;
    if success then
      Copy(theta,old_theta);
      Copy(s.sol.t,old_t);
      Copy(s.sol.v,old_solution);
      Add(w_all_sum,s.sol.v);
      n_all_sum := n_all_sum + 1;
      tmp := Is_Multiple(theta,twopi*p.maxstep,tol);
       if tmp /= 0
        then Add(w_sum,s.sol.v); n_sum := n_sum + 1;
       end if;
       w_c := Is_Multiple(theta,twopi,tol);
       if w_c /= 0 then
         if Equal(s.sol.v,start_solution,epslop) then
           two := Create(integer(2));
           for i in w_sum'range loop
             acc := s.sol.v(i)/two;
             Sub(w_sum(i),acc);
             Sub(w_all_sum(i),acc);
             Clear(acc);
           end loop;
           Clear(two);
           stop := true;
         else
           stop := (w_c >= max_wc);
         end if;
       end if;
     else 
       Copy(old_theta,theta);
       Copy(old_t,s.sol.t);
       Copy(old_solution,s.sol.v);
    end if;
    if not stop
     then Circular_Step_Control(success,p,twopi,step,nsuccess);
    end if;
  end Circular_Management;

-- UPDATE OF PATH DIRECTION :

  procedure Update_Direction
                ( proj : in boolean;
                  freqcnt,defer : in out natural32;
                  r,m,estm : in out integer32;
                  cntm : in out natural32; thresm : in natural32;
                  er : in out integer32;
                  t,target : in Complex_Number; x : in Vector; 
                  dt,s,logs : in out Multprec_Floating_Vectors.Vector;
                  logx,wvl0,wvl1,wvl2 : in out VecVec;
                  v,errv : in out Multprec_Floating_Vectors.Vector;
                  err : in out Floating_Number ) is

  -- DESCRIPTION :
  --   Computes an approximation of the direction of the path.

  -- ON ENTRY :
  --   file       to write intermediate data on, may be omitted;
  --   proj       whether solution vector lies in projective space;
  --   freqcnt    counts the frequency of calls;
  --   defer      only if freqcnt = defer, calculations are done;
  --   r          current order in extrapolation formula;
  --   m          current value for multiplicity;
  --   estm       current estimated for multiplicity;
  --   cntm       number of consecutive times estm has been guessed;
  --   thresm     threshold for augmenting m to estm;
  --   er         order of extrapolator on the errors;
  --   t          current value of continuation parameter;
  --   target     target value of continuation parameter;
  --   x          current solution vector;
  --   dt         vector with distances to target;
  --   s          s-values w.r.t. the current value m;
  --   logs       logarithms of the s-values;
  --   logx       logarithms of the solution vectors;
  --   wvl0       consecutive estimates for previous direction;
  --   wvl1       consecutive estimates for current direction;
  --   wvl2       used as work space for wvl0 and wvl1;
  --   v          current approximate direction of the path;
  --   errv       vector of errors used to estimate m;
  --   err        norm of errv.

  -- ON RETURN :
  --   All in-out variables are updated, provided freqcnt = defer.

  begin
    if freqcnt >= defer
     then
       freqcnt := 0;
      -- if proj
      --  then Projective_Update_Direction
      --         (r,m,estm,cntm,thresm,er,t,target,x,dt,s,logs,logx,v,errv,err);
      --  else Affine_Update_Direction
      --         (r,m,estm,cntm,thresm,er,t,target,x,
      --          dt,s,logs,logx,wvl0,wvl1,wvl2,v,errv,err);
      -- end if;
      -- defer := defer + 1;  -- that is asking for troubles !
     else
       freqcnt := freqcnt + 1;
    end if;
  end Update_Direction;

  procedure Update_Direction
                ( file : in file_type; proj : in boolean;
                  freqcnt,defer : in out natural32;
                  r,m,estm : in out integer32; cntm : in out natural32;
                  thresm : in natural32; er : in out integer32;
                  t,target : in Complex_Number; x : in Vector; 
                  dt,s,logs : in out Multprec_Floating_Vectors.Vector;
                  logx,wvl0,wvl1,wvl2 : in out VecVec;
                  v,errv : in out Multprec_Floating_Vectors.Vector;
                  err : in out Floating_Number ) is

  -- DESCRIPTION :
  --   Computes an approximation of the direction of the path and produces
  --   intermediate output to the file.

  begin
    if freqcnt >= defer
     then
       freqcnt := 0;
      -- if proj
      --  then Projective_Update_Direction
      --         (file,r,m,estm,cntm,thresm,er,t,target,x,
      --               dt,s,logs,logx,v,errv,err);
      --  else Affine_Update_Direction
      --         (file,r,m,estm,cntm,thresm,er,t,target,x,
      --               dt,s,logs,logx,wvl0,wvl1,wvl2,v,errv,err);
      -- end if;
      -- defer := defer + 1;  -- asking for troubles !
       put(file,"direction : "); put(file,v); new_line(file);
       put(file,"difference to old direction : "); put(file,err);
       new_line(file);
       put(file,"++ current m : "); put(file,m,1); put(file," ++ "); 
       put(file,cntm,1); put(file," times estimated m : "); put(file,estm,1);
       put(file," ++ threshold : "); put(file,thresm,1); put_line(file," ++");
       new_line(file);
     else
       freqcnt := freqcnt + 1;
    end if;
  end Update_Direction;

-- LINEAR PATH FOLLOWING FOR ONE PATH :

  procedure Linear_Single_Normal_Silent_Continue
              ( s : in out Solu_Info; target : in Complex_Number;
                tol : in Floating_Number; proj : in boolean;
                p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t : Complex_Number;
    old_solution,prev_solution,old_v,prev_v,vv : Vector(s.sol.v'range);
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
    Linear_Single_Initialize
      (s,p,old_t,prev_t,prev_v,old_solution,prev_solution);
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success) 
                                           and (s.niter <= c.maxtot) loop
      Predictor(s,p,true,prev_solution,prev_v,vv,prev_t,target,step,tol,trial);
      Corrector(s,c);
      Linear_Single_Management(s,p,c,old_t,prev_t,old_solution,prev_solution,
                               old_v,prev_v,vv,step,nsuccess,trial,success);
      if Stop(p,s.sol.t,target,step) then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Corrector(s,cp);
    end;
    Clear(step);
    Clear(old_t); Clear(prev_t);
    Clear(old_solution); Clear(prev_solution);
    Clear(old_v); Clear(prev_v); Clear(vv);
  end Linear_Single_Normal_Silent_Continue;

  procedure Linear_Single_Normal_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in Complex_Number; tol : in Floating_Number;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t : Complex_Number;
    old_solution,prev_solution,old_v,prev_v,vv : Vector(s.sol.v'range);
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
    Linear_Single_Initialize
      (s,p,old_t,prev_t,prev_v,old_solution,prev_solution);
    sWrite(file,s.sol.all);
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success) 
                                           and (s.niter <= c.maxtot) loop
      Predictor(s,p,true,prev_solution,prev_v,vv,prev_t,target,step,tol,trial);
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
      Linear_Single_Management(s,p,c,old_t,prev_t,old_solution,prev_solution,
                               old_v,prev_v,vv,step,nsuccess,trial,success);
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
    Clear(old_solution); Clear(prev_solution);
    Clear(old_v); Clear(prev_v); Clear(vv);
  end Linear_Single_Normal_Reporting_Continue;

  procedure Linear_Single_Conditioned_Silent_Continue
              ( s : in out Solu_Info; target : in Complex_Number;
                tol : in Floating_Number; proj : in boolean;
                rtoric : in integer32; 
                v : in out Multprec_Floating_Vectors.Link_to_Vector;
                errorv : in out Floating_Number;
                p : in Pred_Pars; c : in Corr_Pars ) is

    old_t,prev_t : Complex_Number;
    old_solution,prev_solution,old_v,prev_v,vv : Vector(s.sol.v'range);
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
    Linear_Single_Initialize
      (s,p,old_t,prev_t,prev_v,old_solution,prev_solution);
    if rtoric > 0
     then s.sol.m := m; 
    end if;
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success)
                                           and (s.niter <= c.maxtot) loop
      if (rtoric > 0) 
       then
         if success and then s.rcond > tolsing
                    and then (errorv < 100.0) -- avoid divergence
          then Update_Direction
                 (proj,fcnt,defer,r,s.sol.m,estm,cntm,thresm,er,s.sol.t,target,
                  s.sol.v,dt,ds,logs,logx,wvl0,wvl1,wvl2,v.all,errv,errorv);
          else er := -2;
         end if;
      end if;
      Predictor(s,p,true,prev_solution,prev_v,vv,prev_t,target,step,tol,trial);
      Corrector(s,c);
      Linear_Single_Management(s,p,c,old_t,prev_t,old_solution,prev_solution,
                               old_v,prev_v,vv,step,nsuccess,trial,success);
      if Stop(p,s.sol.t,target,step) then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Corrector(s,cp);
    end;
    Clear(step);
    Clear(old_t); Clear(prev_t);
    Clear(old_solution); Clear(prev_solution);
    Clear(old_v); Clear(prev_v); Clear(vv);
  end Linear_Single_Conditioned_Silent_Continue;

  procedure Linear_Single_Conditioned_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in Complex_Number; tol : in Floating_Number;
                proj : in boolean; rtoric : in integer32;
                v : in out Multprec_Floating_Vectors.Link_to_Vector;
                errorv : in out Floating_Number;
                p : in Pred_Pars; c : in Corr_Pars ) is

    old_t,prev_t : Complex_Number;
    step : Floating_Number;
    nsuccess,trial : natural32 := 0;
    old_solution,prev_solution,old_v,prev_v,vv : Vector(s.sol.v'range);
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
    Linear_Single_Initialize
      (s,p,old_t,prev_t,prev_v,old_solution,prev_solution);
    sWrite(file,s.sol.all);          -- writing the start solution
    if rtoric > 0
     then s.sol.m := m;
    end if;
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success)
                                           and (s.niter <= c.maxtot) loop
      if (rtoric > 0) 
       then
         if success and then s.rcond > tolsing
                    and then (errorv < 100.0) -- avoid divergence
          then Update_Direction(file,
                  proj,fcnt,defer,r,s.sol.m,estm,cntm,thresm,er,s.sol.t,target,
                  s.sol.v,dt,ds,logs,logx,wvl0,wvl1,wvl2,v.all,errv,errorv);
          else er := -2;
         end if;
      end if;
      Predictor(s,p,true,prev_solution,prev_v,vv,prev_t,target,step,tol,trial);
      pWrite(file,step,s.sol.t,s.sol.all);
      Corrector(s,c);
      sWrite(file,s.sol.all);
      Linear_Single_Management(s,p,c,old_t,prev_t,old_solution,prev_solution,
                               old_v,prev_v,vv,step,nsuccess,trial,success);
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
    Clear(old_solution); Clear(prev_solution);
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
    if n_all_sum /= 0
     then acc := Create(n_all_sum);
          for i in w_all_sum'range loop
            all_sum(i) := w_all_sum(i)/acc;
          end loop;
          Clear(acc);
    end if;
    if n_sum /= 0
     then acc := Create(n_sum);
          for i in w_sum'range loop
            sum(i) := w_sum(i)/acc;
          end loop;
          Clear(acc);
     elsif n_all_sum /= 0
         then acc := Create(n_all_sum);
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
    if n_all_sum /= 0
     then acc := Create(n_all_sum);
          for i in w_all_sum'range loop
            all_sum(i) := w_all_sum(i)/acc;
          end loop;
          Clear(acc);
    end if;
    if n_sum /= 0
     then acc := Create(n_sum);
          for i in w_sum'range loop
            sum(i) := w_sum(i)/acc;
          end loop;
          Clear(acc);
     elsif n_all_sum /= 0
         then acc := Create(n_all_sum);
              for i in w_all_sum'range loop
                sum(i) := w_all_sum(i)/acc;
              end loop;
              Clear(acc);
    end if;
  end Circular_Single_Conditioned_Reporting_Continue;

end Multprec_Path_Trackers;
