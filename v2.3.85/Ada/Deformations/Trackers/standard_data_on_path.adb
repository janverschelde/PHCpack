with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Directions_of_Standard_Paths;       use Directions_of_Standard_Paths;

package body Standard_Data_on_Path is

-- UTILITIES :

  function At_End ( t,target : Complex_Number; distance,tol : double_float )
                  return boolean is
  begin
    if distance < tol then
      return Equal(t,target,tol);
    elsif AbsVal(t-target) <= distance then
      return true;
    else 
      return false;
    end if;
  end At_End;

  function Stop ( p : Pred_Pars; t,target : Complex_Number;
                  step : double_float ) return boolean is
  begin
    if step <= p.minstep then
      return true;
    else
      if (p.predictor_type = 2) or (p.predictor_type = 5)
       then return (AbsVal(t-target) <= p.minstep);
       else return false;
      end if;
    end if;
  end Stop;

  function Is_Multiple ( a,b,tol : double_float ) return natural32 is

    fq : double_float;
    iq : integer32;

  begin
    if abs(b) < tol then
      return 0;
    else
      fq := a/b;
      iq := integer32(fq);
      if abs(fq - double_float(iq)) <= tol
       then return natural32(iq);
       else return 0;
      end if;
    end if;
  end Is_Multiple;

  procedure Linear_Step_Control 
              ( success : in boolean; p : in Pred_Pars;
                step : in out double_float;
                nsuccess : in out natural32; trial : in natural32 ) is
  begin
    if (p.predictor_type = 2) or (p.predictor_type = 5) then
      if success then
        nsuccess := nsuccess + 1;
      else
        nsuccess := 0;
        if p.expfac < 1.0
         then step := step + p.expfac*(1.0 - step);
         else step := step + (1.0 - step)/p.expfac;
        end if;
        if step >= 1.0
         then step := p.minstep/2.0;  -- that ends the game
        end if;
      end if;
    elsif success then
      nsuccess := nsuccess + 1;
      if nsuccess > p.success_steps then
        step := step*p.expfac;
        if step > p.maxstep
         then step := p.maxstep;
        end if;
      end if;
    else
      nsuccess := 0;
      if trial mod 3 = 0
       then step := step*p.redfac;
      end if;
    end if;
  end Linear_Step_Control;

  procedure Circular_Step_Control
              ( success : in boolean; p : in Pred_Pars;
                twopi : in double_float;
                step : in out double_float; nsuccess : in out natural32 ) is

    maxstep : constant double_float := p.maxstep*twopi;

  begin
    if success then
      nsuccess := nsuccess + 1;
      if nsuccess > p.success_steps then
        step := step*2.0;              -- expansion factor = 2
        if step > maxstep
         then step := maxstep;
        end if;
      end if;
    else
      nsuccess := 0;
      step := step*0.5;                -- reduction factor = 1/2
    end if;
  end Circular_Step_Control;

  procedure Set_Corrector_Parameters
              ( c : in out Corr_Pars; eps : double_float ) is
  begin
    c.epsrx := eps; c.epsax := eps; c.epsrf := eps; c.epsaf := eps;
  end Set_Corrector_Parameters;

  function End_Game_Corrector_Parameters
             ( current : Corr_Pars; distance,tol : double_float ) 
             return Corr_Pars is

    res : Corr_Pars := current;

  begin
    if distance < tol                 -- correct further to detect clustering
     then res := current;
          Set_Corrector_Parameters(res,tol);
     else res := Create_End_Game;     -- or to move to end game more smoothly
    end if;
    return res;
  end End_Game_Corrector_Parameters;

-- MANAGEMENT OF DATA DURING PATH FOLLOWING :

  procedure Linear_Single_Initialize 
              ( s : in Solu_Info; p : in Pred_Pars;
                old_t,prev_t : out Complex_Number;
                prev_v,old_solution,prev_solution
                  : out Standard_Complex_Vectors.Vector ) is
  begin
    old_t := s.sol.t; old_solution := s.sol.v;
    if p.predictor_type <= 2 or p.predictor_type >= 6 then
      prev_t := s.sol.t;              -- for secant predictors
      prev_solution := s.sol.v;
    end if;
    prev_v := (prev_v'range => Create(0.0));
  end Linear_Single_Initialize;

  procedure Linear_Single_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t : in out Complex_Number;
                old_solution,prev_solution,old_v,prev_v,vv
                  : in out Standard_Complex_Vectors.Vector;
                step : in out double_float;
                nsuccess,trial : in out natural32;
                success : in out boolean ) is
  begin
    success := (((s.resa <= c.epsaf) or else (s.cora <= c.epsax))
        or else ((s.resr <= c.epsrf) or else (s.corr <= c.epsrx)));
    if ((p.predictor_type <= 2) and then success) then  -- secant predictors
      prev_t := old_t;
      prev_solution := old_solution;
    end if;
    if ((p.predictor_type = 6) and then success) then   -- Hermite predictor
      prev_t := old_t;
      prev_solution := old_solution;
      prev_v := old_v;
    end if;
    if ((p.predictor_type = 1) or (p.predictor_type = 3)) then
      if success                                        -- complex predictors
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
    if not success then
      s.sol.t := old_t;
      s.sol.v := old_solution;
    else
      old_t := s.sol.t;
      old_solution := s.sol.v;
      if p.predictor_type = 6
       then old_v := vv;
      end if;
    end if;
  end Linear_Single_Management;

  procedure Linear_Single_Quadratic_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t,prev_t0 : in out Complex_Number;
                old_s,prev_s,prev_s0 : in out Standard_Complex_Vectors.Vector;
                step : in out double_float;
                nsuccess,trial : in out natural32;
                success : in out boolean ) is
  begin
    success := (((s.resa <= c.epsaf) or else (s.cora <= c.epsax))
        or else ((s.resr <= c.epsrf) or else (s.corr <= c.epsrx)));
    if success then 
      prev_t0 := prev_t;
      prev_s0 := prev_s;
      prev_t := old_t;
      prev_s := old_s;
    end if;
    s.nstep := s.nstep + 1;
    if not success
     then s.nfail := s.nfail + 1;
    end if;
    Linear_Step_Control(success,p,step,nsuccess,trial);
    if step < p.minstep
     then return;
    end if;
    if not success then
      s.sol.t := old_t;
      s.sol.v := old_s;
    else
      old_t := s.sol.t;
      old_s := s.sol.v;
    end if;
  end Linear_Single_Quadratic_Management;

  procedure Linear_Single_Cubic_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t,prev_t1,prev_t0 : in out Complex_Number;
                old_s,prev_s,prev_s1,prev_s0
                  : in out Standard_Complex_Vectors.Vector;
                step : in out double_float;
                nsuccess,trial : in out natural32;
                success : in out boolean ) is
  begin
    success := (((s.resa <= c.epsaf) or else (s.cora <= c.epsax))
        or else ((s.resr <= c.epsrf) or else (s.corr <= c.epsrx)));
    if success then 
      prev_t0 := prev_t1; prev_t1 := prev_t;
      prev_s0 := prev_s1; prev_s1 := prev_s;
      prev_t := old_t;
      prev_s := old_s;
    end if;
    s.nstep := s.nstep + 1;
    if not success
     then s.nfail := s.nfail + 1;
    end if;
    Linear_Step_Control(success,p,step,nsuccess,trial);
    if step < p.minstep
     then return;
    end if;
    if not success then
      s.sol.t := old_t;
      s.sol.v := old_s;
    else
      old_t := s.sol.t;
      old_s := s.sol.v;
    end if;
  end Linear_Single_Cubic_Management;

  procedure Linear_Single_Quartic_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prv_t,prv_t2,prv_t1,prv_t0 : in out Complex_Number;
                old_s,prv_s,prv_s2,prv_s1,prv_s0
                  : in out Standard_Complex_Vectors.Vector;
                step : in out double_float;
                nsuccess,trial : in out natural32;
                success : in out boolean ) is
  begin
    success := (((s.resa <= c.epsaf) or else (s.cora <= c.epsax))
        or else ((s.resr <= c.epsrf) or else (s.corr <= c.epsrx)));
    if success then 
      prv_t0 := prv_t1; prv_t1 := prv_t2; prv_t2 := prv_t;
      prv_s0 := prv_s1; prv_s1 := prv_s2; prv_s2 := prv_s;
      prv_t := old_t;
      prv_s := old_s;
    end if;
    s.nstep := s.nstep + 1;
    if not success
     then s.nfail := s.nfail + 1;
    end if;
    Linear_Step_Control(success,p,step,nsuccess,trial);
    if step < p.minstep
     then return;
    end if;
    if not success then
      s.sol.t := old_t;
      s.sol.v := old_s;
    else
      old_t := s.sol.t;
      old_s := s.sol.v;
    end if;
  end Linear_Single_Quartic_Management;

  procedure Linear_Multiple_Initialize 
              ( s : in Solu_Info_Array; p : in Pred_Pars;
                t,old_t,prev_t : out Complex_Number;
                sa,old_sa,prev_sa : in out Solution_Array ) is
  begin
    t := s(s'first).sol.t;
    old_t := s(s'first).sol.t;
    Copy(s,sa); Copy(sa,old_sa);
    case p.predictor_type is
      when 0 | 1 | 2 | 6 => prev_t := s(s'first).sol.t; Copy(sa,prev_sa);
      when others => null;
    end case;
  end Linear_Multiple_Initialize;

  procedure Linear_Multiple_Management
              ( s : in out Solu_Info_array;
                sa,old_sa,prev_sa : in out Solution_Array;
                t,old_t,prev_t : in out Complex_Number; p : in Pred_Pars; 
                step : in out double_float; pivot : in integer32; 
                nsuccess,trial : in out natural32;
                success : in boolean ) is
  begin
    if ((p.predictor_type <= 2) and then success)  -- for secant predictors
     then prev_t := old_t; Copy(old_sa,prev_sa);
    end if;
    if ((p.predictor_type = 1) or (p.predictor_type = 3)) then
      if success                                   -- complex predictors
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
     then Copy(sa,old_sa); old_t := t;
     else Copy(old_sa,sa); t := old_t;
    end if;
  end Linear_Multiple_Management;

  procedure Circular_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t : in out Complex_Number;
                start_solution : in Standard_Complex_Vectors.Vector;
                old_solution,prev_solution,w_sum,w_all_sum
                  : in out Standard_Complex_Vectors.Vector;
                twopi,epslop,tol : in double_float;
                theta,old_theta,step : in out double_float;
                nsuccess,n_sum,n_all_sum,w_c : in out natural32;
                max_wc : in natural32; stop,success : in out boolean ) is

    tmp : natural32;
    use Standard_Complex_Vectors;

  begin
    success := (((s.resa <= c.epsaf) or else (s.cora <= c.epsax))
        or else ((s.resr <= c.epsrf) or else (s.corr <= c.epsrx)));
    if p.predictor_type = 0 and then success then
      prev_t := old_t;
      prev_solution := old_solution;
    end if;
    s.nstep := s.nstep + 1;
    if not success then s.nfail := s.nfail + 1; end if;
    if success then
      old_theta := theta; old_t := s.sol.t;
      old_solution := s.sol.v;
      w_all_sum := w_all_sum + s.sol.v;
      n_all_sum := n_all_sum + 1;
      tmp := Is_Multiple(theta,twopi*p.maxstep,tol);
      if tmp /= 0
       then w_sum := w_sum + s.sol.v; n_sum := n_sum + 1;
      end if;
      w_c := Is_Multiple(theta,twopi,tol);
      if w_c /= 0 then
        if Equal(s.sol.v,start_solution,epslop) then
          w_sum := w_sum - s.sol.v * Create(0.5);
          w_all_sum := w_all_sum - s.sol.v * Create(0.5);
          stop := true;
        else
          stop := (w_c >= max_wc);
        end if;
      end if;
    else
      theta := old_theta;
      s.sol.t := old_t;
      s.sol.v := old_solution;
    end if;
    if not stop
     then Circular_Step_Control(success,p,twopi,step,nsuccess);
    end if;
  end Circular_Management;

-- UPDATE OF PATH DIRECTION :

  procedure Update_Direction
              ( proj : in boolean; freqcnt,defer : in out natural32;
                r,m,estm : in out integer32;
                cntm : in out natural32; thresm : in natural32;
                er : in out integer32;
                t,target : in Complex_Number;
                x : in Standard_Complex_Vectors.Vector; 
                dt,s,logs : in out Standard_Floating_Vectors.Vector;
                logx,wvl0,wvl1,wvl2 : in out Standard_Floating_VecVecs.VecVec;
                v,errv : in out Standard_Floating_Vectors.Vector;
                err : in out double_float ) is

  begin
    if freqcnt >= defer then
      freqcnt := 0;
      if proj then
        Projective_Update_Direction
          (r,m,estm,cntm,thresm,er,t,target,x,dt,s,logs,logx,v,errv,err);
      else
        Affine_Update_Direction
          (r,m,estm,cntm,thresm,er,t,target,x,
           dt,s,logs,logx,wvl0,wvl1,wvl2,v,errv,err);
      end if;
     -- defer := defer + 1;  -- that is asking for troubles !
    else
       freqcnt := freqcnt + 1;
    end if;
  end Update_Direction;

  procedure Update_Direction
              ( file : in file_type; proj : in boolean;
                freqcnt,defer : in out natural32;
                r,m,estm : in out integer32;
                cntm : in out natural32; thresm : in natural32;
                er : in out integer32;
                t,target : in Complex_Number;
                x : in Standard_Complex_Vectors.Vector; 
                dt,s,logs : in out Standard_Floating_Vectors.Vector;
                logx,wvl0,wvl1,wvl2 : in out Standard_Floating_VecVecs.VecVec;
                v,errv : in out Standard_Floating_Vectors.Vector;
                err : in out double_float ) is
  begin
    if freqcnt >= defer then
      freqcnt := 0;
      if proj then
        Projective_Update_Direction
          (file,r,m,estm,cntm,thresm,er,t,target,x,dt,s,logs,logx,v,errv,err);
      else
        Affine_Update_Direction
          (file,r,m,estm,cntm,thresm,er,t,target,x,
           dt,s,logs,logx,wvl0,wvl1,wvl2,v,errv,err);
      end if;
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

end Standard_Data_on_Path;
