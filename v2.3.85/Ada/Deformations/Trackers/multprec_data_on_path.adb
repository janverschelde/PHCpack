with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Floating_Vectors_io;       use Multprec_Floating_Vectors_io;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Continuation_Parameters;
--with Multprec_Directions_of_Paths;       use Multprec_Directions_of_Paths;

package body Multprec_Data_on_Path is

-- UTILITIES :

  function At_End ( t,target : Complex_Number;
                    distance,tol : Floating_Number ) return boolean is

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

    res : boolean;
    diff : Complex_Number;
    absdiff : Floating_Number;

  begin
    if step < p.minstep then
      res := true;
    else
      if (p.predictor_type = 2) or (p.predictor_type = 5) then
        diff := t-target;
        absdiff := AbsVal(diff);
        res := (absdiff < p.minstep);
        Clear(diff); Clear(absdiff);
      else
        res := false;
      end if;
    end if;
    return res;
  end Stop;

  function Is_Multiple ( a,b,tol : Floating_Number ) return natural32 is

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

    maxstep : constant Floating_Number := p.maxstep*twopi;

  begin
    if success then
      nsuccess := nsuccess + 1;
      if nsuccess > p.success_steps then
        Mul(step,2.0);             -- expansion factor = 2
        if step > maxstep
         then Copy(maxstep,step);
        end if;
      end if;
    else
      nsuccess := 0;
      Div(step,2.0);                   -- reduction factor = 1/2
    end if;
  end Circular_Step_Control;

  procedure Set_Corrector_Parameters
              ( c : in out Corr_Pars; eps : Floating_Number ) is
  begin
    Copy(eps,c.epsrx); Copy(eps,c.epsax);
    Copy(eps,c.epsrf); Copy(eps,c.epsaf);
  end Set_Corrector_Parameters;

  function End_Game_Corrector_Parameters
             ( current : Corr_Pars; distance,tol : Floating_Number ) 
             return Corr_Pars is

    res : Corr_Pars;

  begin
    if distance < tol then           -- correct further to detect clustering
      Copy(current.epsrx,res.epsrx);
      Copy(current.epsax,res.epsax);
      Copy(current.epsrf,res.epsrf);
      Copy(current.epsaf,res.epsaf);
      res.maxit := current.maxit;
      res.maxtot := current.maxtot;
      Set_Corrector_Parameters(res,tol);
    else 
      res := Convert(Continuation_Parameters.Create_End_Game);
      -- or to move smoothly to end game
    end if;
    return res;
  end End_Game_Corrector_Parameters;

-- MANAGEMENT OF DATA DURING PATH FOLLOWING :

  procedure Linear_Single_Initialize 
              ( s : in Solu_Info; p : in Pred_Pars;
                old_t,prev_t : out Complex_Number;
                prev_v,old_solution,prev_solution
                  : out Multprec_Complex_Vectors.Vector ) is
  begin
    Copy(s.sol.t,old_t);
    Multprec_Complex_Vectors.Copy(s.sol.v,old_solution);
    if p.predictor_type <= 2 or p.predictor_type >= 6 then
      Copy(s.sol.t,prev_t);
      Multprec_Complex_Vectors.Copy(s.sol.v,prev_solution);
    end if;
    prev_v := (prev_v'range => Create(integer(0)));
  end Linear_Single_Initialize;

  procedure Linear_Single_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t : in out Complex_Number;
                old_solution,prev_solution,old_v,prev_v,vv
                  : in out Multprec_Complex_Vectors.Vector;
                step : in out Floating_Number;
                nsuccess,trial : in out natural32;
                success : in out boolean ) is
  begin
    success := (((s.resa < c.epsaf) or else (s.cora < c.epsax))
        or else ((s.resr < c.epsrf) or else (s.corr < c.epsrx)));
    if ((p.predictor_type <= 2) and then success) then    -- secant predictors
      Copy(old_t,prev_t);
      Multprec_Complex_Vectors.Copy(old_solution,prev_solution);
    end if;
    if ((p.predictor_type = 6) and then success) then     -- Hermite predictor
      Copy(old_t,prev_t);
      Multprec_Complex_Vectors.Copy(old_solution,prev_solution);
      Multprec_Complex_Vectors.Copy(old_v,prev_v);
    end if;
    if ((p.predictor_type = 1) or (p.predictor_type = 3)) then
      if success                                         -- complex predictors
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
      Copy(old_t,s.sol.t);
      Multprec_Complex_Vectors.Copy(old_solution,s.sol.v);
    else
      Copy(s.sol.t,old_t);
      Multprec_Complex_Vectors.Copy(s.sol.v,old_solution);
      if p.predictor_type = 6
       then Multprec_Complex_Vectors.Copy(vv,old_v);
      end if;
    end if;
  end Linear_Single_Management;

  procedure Linear_Single_Quadratic_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t,prev_t0 : in out Complex_Number;
                old_s,prev_s,prev_s0 : in out Multprec_Complex_Vectors.Vector;
                step : in out Floating_Number;
                nsuccess,trial : in out natural32;
                success : in out boolean ) is
  begin
    success := (((s.resa < c.epsaf) or else (s.cora < c.epsax))
        or else ((s.resr < c.epsrf) or else (s.corr < c.epsrx)));
    if success then 
      Copy(prev_t,prev_t0);
      Multprec_Complex_Vectors.Copy(prev_s,prev_s0);
      Copy(old_t,prev_t);
      Multprec_Complex_Vectors.Copy(old_s,prev_s);
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
      Copy(old_t,s.sol.t);
      Multprec_Complex_Vectors.Copy(old_s,s.sol.v);
    else
      Copy(s.sol.t,old_t);
      Multprec_Complex_Vectors.Copy(s.sol.v,old_s);
    end if;
  end Linear_Single_Quadratic_Management;

  procedure Linear_Single_Cubic_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t,prev_t1,prev_t0 : in out Complex_Number;
                old_s,prev_s,prev_s1,prev_s0
                  : in out Multprec_Complex_Vectors.Vector;
                step : in out Floating_Number;
                nsuccess,trial : in out natural32;
                success : in out boolean ) is
  begin
    success := (((s.resa < c.epsaf) or else (s.cora < c.epsax))
        or else ((s.resr < c.epsrf) or else (s.corr < c.epsrx)));
    if success then 
      Copy(prev_t1,prev_t0); Copy(prev_t,prev_t1);
      Multprec_Complex_Vectors.Copy(prev_s1,prev_s0);
      Multprec_Complex_Vectors.Copy(prev_s,prev_s1);
      Copy(old_t,prev_t);
      Multprec_Complex_Vectors.Copy(old_s,prev_s);
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
      Copy(old_t,s.sol.t);
      Multprec_Complex_Vectors.Copy(old_s,s.sol.v);
    else
      Copy(s.sol.t,old_t);
      Multprec_Complex_Vectors.Copy(s.sol.v,old_s);
    end if;
  end Linear_Single_Cubic_Management;

  procedure Linear_Single_Quartic_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prv_t,prv_t2,prv_t1,prv_t0 : in out Complex_Number;
                old_s,prv_s,prv_s2,prv_s1,prv_s0
                  : in out Multprec_Complex_Vectors.Vector;
                step : in out Floating_Number;
                nsuccess,trial : in out natural32;
                success : in out boolean ) is
  begin
    success := (((s.resa < c.epsaf) or else (s.cora < c.epsax))
        or else ((s.resr < c.epsrf) or else (s.corr < c.epsrx)));
    if success then 
      Copy(prv_t1,prv_t0); Multprec_Complex_Vectors.Copy(prv_s1,prv_s0);
      Copy(prv_t2,prv_t1); Multprec_Complex_Vectors.Copy(prv_s2,prv_s1);
      Copy(prv_t,prv_t2);  Multprec_Complex_Vectors.Copy(prv_s,prv_s2);
      Copy(old_t,prv_t);
      Multprec_Complex_Vectors.Copy(old_s,prv_s);
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
      Copy(old_t,s.sol.t);
      Multprec_Complex_Vectors.Copy(old_s,s.sol.v);
    else
      Copy(s.sol.t,old_t);
      Multprec_Complex_Vectors.Copy(s.sol.v,old_s);
    end if;
  end Linear_Single_Quartic_Management;

  procedure Linear_Multiple_Initialize 
              ( s : in Solu_Info_Array; p : in Pred_Pars;
                t,old_t,prev_t : out Complex_Number;
                sa,old_sa,prev_sa : in out Solution_Array ) is
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
  begin
    if ((p.predictor_type <= 2) and then success) then -- for secant predictors
      Copy(old_t,prev_t);
      Copy(old_sa,prev_sa);
    end if;
    if ((p.predictor_type = 1) or (p.predictor_type = 3)) then
      if success -- complex predictors
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
                start_solution : in Multprec_Complex_Vectors.Vector;
                old_solution,prev_solution,w_sum,w_all_sum
                  : in out Multprec_Complex_Vectors.Vector;
                twopi,epslop,tol : in Floating_Number;
                theta,old_theta,step : in out Floating_Number;
                nsuccess,n_sum,n_all_sum,w_c : in out natural32;
                max_wc : in natural32; stop,success : in out boolean ) is

    tmp : natural32;
    two,acc : Complex_Number;

  begin
    success := (((s.resa < c.epsaf) or else (s.cora < c.epsax))
        or else ((s.resr < c.epsrf) or else (s.corr < c.epsrx)));
    if p.predictor_type = 0 and then success then
      Copy(old_t,prev_t);
      Multprec_Complex_Vectors.Copy(old_solution,prev_solution);
    end if;
    s.nstep := s.nstep + 1;
    if not success then s.nfail := s.nfail + 1; end if;
    if success then
      Copy(theta,old_theta);
      Copy(s.sol.t,old_t);
      Multprec_Complex_Vectors.Copy(s.sol.v,old_solution);
      Multprec_Complex_Vectors.Add(w_all_sum,s.sol.v);
      n_all_sum := n_all_sum + 1;
      tmp := Is_Multiple(theta,twopi*p.maxstep,tol);
       if tmp /= 0 then
         Multprec_Complex_Vectors.Add(w_sum,s.sol.v);
         n_sum := n_sum + 1;
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
       Multprec_Complex_Vectors.Copy(old_solution,s.sol.v);
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
                t,target : in Complex_Number;
                x : in Multprec_Complex_Vectors.Vector; 
                dt,s,logs : in out Multprec_Floating_Vectors.Vector;
                logx,wvl0,wvl1,wvl2 : in out Multprec_Floating_VecVecs.VecVec;
                v,errv : in out Multprec_Floating_Vectors.Vector;
                err : in out Floating_Number ) is
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
                t,target : in Complex_Number;
                x : in Multprec_Complex_Vectors.Vector; 
                dt,s,logs : in out Multprec_Floating_Vectors.Vector;
                logx,wvl0,wvl1,wvl2 : in out Multprec_Floating_VecVecs.VecVec;
                v,errv : in out Multprec_Floating_Vectors.Vector;
                err : in out Floating_Number ) is
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

end Multprec_Data_on_Path;
