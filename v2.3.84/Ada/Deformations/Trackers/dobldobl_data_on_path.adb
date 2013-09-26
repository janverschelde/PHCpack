with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Double_Double_Vectors_io;           use Double_Double_Vectors_io;
with DoblDobl_Complex_Equality_Tests;    use DoblDobl_Complex_Equality_Tests;
with Directions_of_DoblDobl_Paths;       use Directions_of_DoblDobl_Paths;
with Standard_Data_on_Path;

package body DoblDobl_Data_on_Path is

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

-- MANAGEMENT OF DATA DURING PATH FOLLOWING :

  procedure Linear_Single_Initialize 
              ( s : in Solu_Info; p : in Pred_Pars;
                old_t,prev_t : out Complex_Number;
                prev_v,old_solution,prev_solution
                  : out DoblDobl_Complex_Vectors.Vector ) is

    zero : constant double_double := create(0.0);

  begin
    old_t := s.sol.t; old_solution := s.sol.v;
    if p.predictor_type <= 2 or p.predictor_type >= 6 -- for secant predictors
     then prev_t := s.sol.t; prev_solution := s.sol.v;
    end if;
    prev_v := (prev_v'range => Create(zero));
  end Linear_Single_Initialize;

  procedure Linear_Single_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t : in out Complex_Number;
                old_solution,prev_solution,old_v,prev_v,vv
                  : in out DoblDobl_Complex_Vectors.Vector;
                step : in out double_float;
                nsuccess,trial : in out natural32;
                success : in out boolean ) is
  begin
    success := (((s.resa <= c.epsaf) or else (s.cora <= c.epsax))
        or else ((s.resr <= c.epsrf) or else (s.corr <= c.epsrx)));
    if ((p.predictor_type <= 2) and then success)       -- secant predictors
     then prev_t := old_t; prev_solution := old_solution;
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
    Standard_Data_on_Path.Linear_Step_Control(success,p,step,nsuccess,trial);
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
                old_s,prev_s,prev_s0 : in out DoblDobl_Complex_Vectors.Vector;
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
    Standard_Data_on_Path.Linear_Step_Control(success,p,step,nsuccess,trial);
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
                  : in out DoblDobl_Complex_Vectors.Vector;
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
    Standard_Data_on_Path.Linear_Step_Control(success,p,step,nsuccess,trial);
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
                  : in out DoblDobl_Complex_Vectors.Vector;
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
    Standard_Data_on_Path.Linear_Step_Control(success,p,step,nsuccess,trial);
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
    if ((p.predictor_type <= 2) and then success)      -- for secant predictors
     then prev_t := old_t; Copy(old_sa,prev_sa);
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
    Standard_Data_on_Path.Linear_Step_Control(success,p,step,nsuccess,trial);
    if step < p.minstep then return; end if;
    if success
     then Copy(sa,old_sa); old_t := t;
     else Copy(old_sa,sa); t := old_t;
    end if;
  end Linear_Multiple_Management;

  procedure Circular_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t : in out Complex_Number;
                start_solution : in DoblDobl_Complex_Vectors.Vector;
                old_solution,prev_solution,w_sum,w_all_sum
                  : in out DoblDobl_Complex_Vectors.Vector;
                twopi,epslop,tol : in double_float;
                theta,old_theta,step : in out double_float;
                nsuccess,n_sum,n_all_sum,w_c : in out natural32;
                max_wc : in natural32; stop,success : in out boolean ) is

    tmp : natural32;
    a_half : constant double_double := create(0.5);
    use DoblDobl_Complex_Vectors;

  begin
    success := (((s.resa <= c.epsaf) or else (s.cora <= c.epsax))
        or else ((s.resr <= c.epsrf) or else (s.corr <= c.epsrx)));
    if p.predictor_type = 0 and then success
     then prev_t := old_t;
          prev_solution := old_solution;
    end if;
    s.nstep := s.nstep + 1;
    if not success then s.nfail := s.nfail + 1; end if;
    if success then
      old_theta := theta; old_t := s.sol.t;
      old_solution := s.sol.v;
      w_all_sum := w_all_sum + s.sol.v;
      n_all_sum := n_all_sum + 1;
      tmp := Standard_Data_on_Path.Is_Multiple(theta,twopi*p.maxstep,tol);
      if tmp /= 0
       then w_sum := w_sum + s.sol.v; n_sum := n_sum + 1;
      end if;
      w_c := Standard_Data_on_Path.Is_Multiple(theta,twopi,tol);
      if w_c /= 0 then
        if Equal(s.sol.v,start_solution,epslop)
         then w_sum := w_sum - s.sol.v * Create(a_half);
              w_all_sum := w_all_sum - s.sol.v * Create(a_half);
              stop := true;
         else stop := (w_c >= max_wc);
        end if;
      end if;
    else
      theta := old_theta;
      s.sol.t := old_t;
      s.sol.v := old_solution;
    end if;
    if not stop then
      Standard_Data_on_Path.Circular_Step_Control
        (success,p,twopi,step,nsuccess);
    end if;
  end Circular_Management;

-- UPDATE OF PATH DIRECTION :

  procedure Update_Direction
              ( proj : in boolean;
                freqcnt,defer : in out natural32; 
                r,m,estm : in out integer32; cntm : in out natural32;
                thresm : in natural32; er : in out integer32;
                t,target : in Complex_Number;
                x : in DoblDobl_Complex_Vectors.Vector; 
                dt,s,logs : in out Double_Double_Vectors.Vector;
                logx,wvl0,wvl1,wvl2 : in out Double_Double_VecVecs.VecVec;
                v,errv : in out Double_Double_Vectors.Vector;
                err : in out double_double ) is
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
               r,m,estm : in out integer32; cntm : in out natural32;
               thresm : in natural32; er : in out integer32;
               t,target : in Complex_Number;
               x : in DoblDobl_Complex_Vectors.Vector; 
               dt,s,logs : in out Double_Double_Vectors.Vector;
               logx,wvl0,wvl1,wvl2 : in out Double_Double_VecVecs.VecVec;
               v,errv : in out Double_Double_Vectors.Vector;
               err : in out double_double ) is

  -- DESCRIPTION :
  --   Computes an approximation of the direction of the path and produces
  --   intermediate output to the file.

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

end DoblDobl_Data_on_Path;
