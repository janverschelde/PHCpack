with Multprec_Predictors;                use Multprec_Predictors;

package body Multprec_Dispatch_Predictors is

  procedure Single_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars; xt : in boolean;
                prev_x,prev_v : in Vector; v : in out Vector;
                prev_t,target : in Complex_Number;
                step,tol : in Floating_Number; trial : in out natural32 ) is

    h,zero : Floating_Number;

    procedure TR_Predictor is new Tangent_Single_Real_Predictor(Norm,dH,dH);
    procedure TC_Predictor is new Tangent_Single_Complex_Predictor(Norm,dH,dH);
    procedure TG_Predictor is new Tangent_Geometric_Predictor(Norm,dH,dH);
    procedure HR_Predictor is new Hermite_Single_Real_Predictor(Norm,dH,dH);

  begin
    if not xt
     then
       case p.predictor_type is
         when 0 | 3 | 6 => Real_Predictor(s.sol.t,target,step,tol,p.power,h);
         when 1 | 4     => zero := Create(integer(0));
                           Complex_Predictor
                             (s.sol.t,target,step,tol,h,zero,trial);
                           Clear(zero);
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

  procedure Multiple_Predictor
              ( s : in out Solu_Info_Array; p : in Pred_Pars; xt : in boolean;
                sa : in out Solution_Array; prev_sa : in Solution_Array;
                t : in out Complex_Number; prev_t,target : in Complex_Number;
                step,tol,dist : in Floating_Number; trial : in natural32 ) is

    cnt : natural32 := 0;
    h : Floating_Number;

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

end Multprec_Dispatch_Predictors;
