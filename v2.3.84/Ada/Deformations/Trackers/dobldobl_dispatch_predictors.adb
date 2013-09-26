with DoblDobl_Predictors;                use DoblDobl_Predictors;

package body DoblDobl_Dispatch_Predictors is

  procedure Single_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars; xt : in boolean;
                prev_x,prev_v : in Vector; v : in out Vector;
                prev_t,target : in Complex_Number;
                step,tol : in double_float; trial : in out natural32 ) is

    h : double_float;

    procedure TR_Predictor is new Tangent_Single_Real_Predictor(Norm,dH,dH);
    procedure TC_Predictor is new Tangent_Single_Complex_Predictor(Norm,dH,dH);
    procedure TG_Predictor is new Tangent_Geometric_Predictor(Norm,dH,dH);
    procedure HR_Predictor is new Hermite_Single_Real_Predictor(Norm,dH,dH);

  begin
    if not xt
     then
       case p.predictor_type is
         when 0 | 3 | 6 => Real_Predictor(s.sol.t,target,step,tol,p.power,h);
         when 1 | 4     => Complex_Predictor
                             (s.sol.t,target,step,tol,h,0.0,trial);
         when 2 | 5     => Geometric_Predictor(s.sol.t,target,step,tol);
         when others    => null;
       end case;      
     else
       case p.predictor_type is
         when 0 => Secant_Single_Real_Predictor
                     (s.sol.v,prev_x,s.sol.t,prev_t,target,step,tol,p.power);
         when 1 => Secant_Single_Complex_Predictor
                     (s.sol.v,prev_x,s.sol.t,prev_t,target,step,tol,
                      p.dist_target,trial);
         when 2 => Secant_Geometric_Predictor
                     (s.sol.v,prev_x,s.sol.t,prev_t,target,step,tol);
         when 3 => TR_Predictor(s.sol.v,s.sol.t,target,step,tol,p.power);
                   s.nsyst := s.nsyst + 1;
         when 4 => TC_Predictor
                     (s.sol.v,s.sol.t,target,step,tol,p.dist_target,trial);
                   s.nsyst := s.nsyst + 1;
         when 5 => TG_Predictor(s.sol.v,s.sol.t,target,step,tol);
                   s.nsyst := s.nsyst + 1;
         when 6 => HR_Predictor
                     (s.sol.v,prev_x,s.sol.t,prev_t,target,v,prev_v,step,tol);
                   s.nsyst := s.nsyst + 1;
         when others => null;
       end case;
    end if;
  end Single_Predictor;

  function Real_Equal ( a,b : Complex_Number; tol : in double_float )
                      return boolean is

    dd_a_re : constant double_double := REAL_PART(a);
    dd_b_re : constant double_double := REAL_PART(b);
    df_a_re : constant double_float := to_double(dd_a_re);
    df_b_re : constant double_float := to_double(dd_b_re);
    v : constant double_float := abs(df_a_re - df_b_re);

  begin
    if v > tol
     then return false;
     else return true;
    end if;
  end Real_Equal;

  procedure Single_Quadratic_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars; xt : in boolean;
                x1,x0 : in Vector;
                t1,t0,target : in Complex_Number;
                step,tol : in double_float ) is

  -- Assumed is that t0 <= t1 <= t2 = s.sol.t
  -- and the secant predictor is used for bootstrapping.

    h : double_float;
    real_tol : constant double_float := 1.0E-8;

  begin
    if not xt then
      Real_Predictor(s.sol.t,target,step,tol,p.power,h);
    else
      if Real_Equal(t0,t1,real_tol) then
        if Real_Equal(t1,s.sol.t,real_tol) then
          Real_Predictor(s.sol.t,target,step,tol,p.power,h);
        else
          Secant_Single_Real_Predictor
            (s.sol.v,x1,s.sol.t,t1,target,step,tol,p.power);
        end if;
      else -- t0 /= t1
        if Real_Equal(t1,s.sol.t,real_tol) then
          Secant_Single_Real_Predictor
            (s.sol.v,x0,s.sol.t,t0,target,step,tol,p.power);
        else -- t1 /= t2
          Quadratic_Single_Real_Predictor
            (s.sol.v,x1,x0,s.sol.t,t1,t0,target,step,tol,p.power);
        end if;
      end if;
    end if;
  end Single_Quadratic_Predictor;

  procedure Single_Cubic_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars; xt : in boolean;
                x2,x1,x0 : in Vector;
                t2,t1,t0,target : in Complex_Number;
                step,tol : in double_float ) is

  -- Assumed is that t0 <= t1 <= t2 <= t3 = s.sol.t
  -- and the quadratic predictor is used for bootstrapping.

    h : double_float;
    real_tol : constant double_float := 1.0E-8;

  begin
    if not xt then
      Real_Predictor(s.sol.t,target,step,tol,p.power,h);
    else
      if Real_Equal(t0,t1,real_tol) then -- bootstrap: skip t0
        Single_Quadratic_Predictor(s,p,xt,x2,x1,t2,t1,target,step,tol);
      else -- t0 /= t1
        if Real_Equal(t1,t2,real_tol) then -- bootstrap: skip t1
          Single_Quadratic_Predictor(s,p,xt,x2,x0,t2,t0,target,step,tol);
        else -- t1 /= t2
          if Real_Equal(t2,s.sol.t,real_tol) then -- skip t2
            Quadratic_Single_Real_Predictor
              (s.sol.v,x1,x0,s.sol.t,t1,t0,target,step,tol,p.power);
          else -- t2 /= t3 
            Cubic_Single_Real_Predictor
              (s.sol.v,x2,x1,x0,s.sol.t,t2,t1,t0,target,step,tol,p.power);
          end if;
        end if;
      end if;
    end if;
  end Single_Cubic_Predictor;

  procedure Single_Quartic_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars; xt : in boolean;
                x3,x2,x1,x0 : in Vector;
                t3,t2,t1,t0,target : in Complex_Number;
                step,tol : in double_float ) is

  -- Assumed is that t0 <= t1 <= t2 <= t3 <= t4 = s.sol.t
  -- and the cubic predictor is used for bootstrapping.

    h : double_float;
    real_tol : constant double_float := 1.0E-8;

  begin
    if not xt then
      Real_Predictor(s.sol.t,target,step,tol,p.power,h);
    else
      if Real_Equal(t0,t1,real_tol) then -- bootstrap: skip t0
        Single_Cubic_Predictor(s,p,xt,x3,x2,x1,t3,t2,t1,target,step,tol);
      else -- t0 /= t1
        if Real_Equal(t1,t2,real_tol) then -- bootstrap: skip t1
          Single_Cubic_Predictor(s,p,xt,x3,x2,x0,t3,t2,t0,target,step,tol);
        else -- t1 /= t2
          if Real_Equal(t2,t3,real_tol) then -- bootstrap: skip t2
            Single_Cubic_Predictor(s,p,xt,x3,x1,x0,t3,t1,t0,target,step,tol);
          else -- t2 /= t3
            if Real_Equal(t3,s.sol.t,real_tol) then -- skip t3
              Cubic_Single_Real_Predictor
                (s.sol.v,x2,x1,x0,s.sol.t,t2,t1,t0,target,step,tol,p.power);
            else -- t3 /= t4
              Quartic_Single_Real_Predictor
                (s.sol.v,x3,x2,x1,x0,s.sol.t,t3,t2,t1,t0,
                 target,step,tol,p.power);
            end if;
          end if;
        end if;
      end if;
    end if;
  end Single_Quartic_Predictor;

  procedure Multiple_Predictor
              ( s : in out Solu_Info_Array; p : in Pred_Pars; xt : in boolean;
                sa : in out Solution_Array; prev_sa : in Solution_Array;
                t : in out Complex_Number; prev_t,target : in Complex_Number;
                step,tol,dist : in double_float; trial : in natural32 ) is

    cnt : natural32 := 0;
    h : double_float;

    procedure TR_Predictor is new Tangent_Multiple_Real_Predictor(Norm,dH,dH);
    procedure TC_Predictor is
      new Tangent_Multiple_Complex_Predictor(Norm,dH,dH);

  begin
    if not xt
     then
       case p.predictor_type is
         when 0 | 3 | 6 => Real_Predictor(t,target,step,tol,p.power,h);
         when 1 | 4     => Complex_Predictor(t,target,step,tol,h,dist,trial);
         when 2 | 5     => Geometric_Predictor(t,target,step,tol);
         when others    => null;
       end case;
     else
       case p.predictor_type is
         when 0 => Secant_Multiple_Real_Predictor
                     (sa,prev_sa,t,prev_t,target,step,tol,dist,p.power);
         when 1 => Secant_Multiple_Complex_Predictor
                     (sa,prev_sa,t,prev_t,target,step,tol,dist,
                      p.dist_target,trial);
         when 2 => Geometric_Predictor(t,target,step,tol);
         when 3 => TR_Predictor(sa,t,target,step,tol,dist,cnt,p.power);
                   for k in s'range loop
                     s(k).nsyst := s(k).nsyst + 1;
                   end loop;
         when 4 => TC_Predictor
                     (sa,t,target,step,tol,dist,p.dist_target,trial,cnt);
                   for k in s'range loop
                     s(k).nsyst := s(k).nsyst + 1;
                   end loop;
         when 5 => Geometric_Predictor(t,target,step,tol);
         when 6 => Real_Predictor(t,target,step,tol,p.power,h);
         when others => null;
       end case;
    end if;
  end Multiple_Predictor;

end DoblDobl_Dispatch_Predictors;
