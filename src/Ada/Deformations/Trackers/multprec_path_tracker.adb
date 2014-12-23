with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Homotopy;
with Multprec_Dispatch_Predictors;       use Multprec_Dispatch_Predictors;
with Multprec_Correctors;                use Multprec_Correctors;
with Multprec_Continuation_Data;         use Multprec_Continuation_Data;
with Multprec_Data_on_Path;              use Multprec_Data_on_Path;
with Continuation_Parameters;

package body Multprec_Path_Tracker is

-- INTERNAL DATA :

  current : Link_to_Solution;
  point : Solu_Info;
  tol : Floating_Number := create(1.0E-12);
  old_t,prev_t,prev_t0,prev_t1,prev_t2 : Complex_Number;
  old_sol,prev_sol,prev_sol0,prev_sol1,prev_sol2,old_v,prev_v,vv
    : Multprec_Complex_Vectors.Link_to_Vector;
  step : Floating_Number;
  nsuccess,trial : natural32 := 0;
  success : boolean := true;

-- CONSTRUCTORS :

  procedure Clear_Solution_Data is

  -- DESCRIPTION :
  --   For the initialization of a second solution with the same homotopy,
  --   the intermediate solution data must be cleared first.

    use Multprec_Complex_Vectors;

  begin
    if old_sol /= null then Clear(old_sol); end if;
    if prev_sol /= null then Clear(prev_sol); end if;
    if prev_sol0 /= null then Clear(prev_sol0); end if;
    if prev_sol1 /= null then Clear(prev_sol1); end if;
    if prev_sol2 /= null then Clear(prev_sol2); end if;
    if old_v /= null then Clear(old_v); end if;
    if prev_v /= null then Clear(prev_v); end if;
    if vv /= null then Clear(vv); end if;
  end Clear_Solution_Data;

  procedure Init_Solution_Data is

  -- DESCRIPTION :
  --   Initializes the solution data: the backup solution,
  --   pointers to four previous solution vectors,
  --   and their corresponding t values.

  begin
    old_sol := new Multprec_Complex_Vectors.Vector(current.v'range);
    prev_sol := new Multprec_Complex_Vectors.Vector(current.v'range);
    prev_sol0 := new Multprec_Complex_Vectors.Vector(current.v'range);
    prev_sol1 := new Multprec_Complex_Vectors.Vector(current.v'range);
    prev_sol2 := new Multprec_Complex_Vectors.Vector(current.v'range);
    old_v := new Multprec_Complex_Vectors.Vector(current.v'range);
    prev_v := new Multprec_Complex_Vectors.Vector(current.v'range);
    vv := new Multprec_Complex_Vectors.Vector(current.v'range);
    step := Create(Continuation_Parameters.max_path_step_size);
    Copy(current.t,old_t);
    Multprec_Complex_Vectors.Copy(current.v,old_sol.all);
    Copy(old_t,prev_t);   
    Copy(prev_t,prev_t0); 
    Copy(prev_t,prev_t1);
    Copy(prev_t,prev_t2);
    Multprec_Complex_Vectors.Copy(old_sol.all,prev_sol.all);
    Multprec_Complex_Vectors.Copy(prev_sol.all,prev_sol0.all);
    Multprec_Complex_Vectors.Copy(prev_sol.all,prev_sol1.all);
    Multprec_Complex_Vectors.Copy(prev_sol.all,prev_sol2.all);
    for i in prev_v'range loop
      declare
        mp_zero : Floating_Number := create(0.0);
      begin
        prev_v(i) := Create(mp_zero);
        Clear(mp_zero);
      end;
    end loop;
  end Init_Solution_Data;

  procedure Init ( s : in Link_to_Solution ) is
  begin
    current := s;
    point := Shallow_Create(current);
    nsuccess := 0;
    trial := 0;
    success := true;
    Clear_Solution_Data;
    Init_Solution_Data;
  end Init;

  procedure Init ( p,q : in Link_to_Poly_Sys; fixed_gamma : in boolean;
                   deci : in natural32 ) is

    mp_re : Floating_Number;
    mp_im : Floating_Number;
    gamma : Complex_Number;

  begin
    if fixed_gamma then
      mp_re := create(0.57670012968461137);
      mp_im := create(0.8169559109411918);
    else
      declare
        st_gamma : constant Standard_Complex_Numbers.Complex_Number
                 := Standard_Random_Numbers.Random1;
        st_re : constant double_float
              := Standard_Complex_Numbers.REAL_PART(st_gamma);
        st_im : constant double_float
              := Standard_Complex_Numbers.IMAG_PART(st_gamma);
      begin
        mp_re := create(st_re);
        mp_im := create(st_im);
      end;
    end if;
    gamma := Create(mp_re,mp_im);
    Init(p,q,gamma,2,deci);
    Clear(mp_re); Clear(mp_im);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_Sys; fixed_gamma : in boolean;
                   s : in Link_to_Solution; deci : in natural32 ) is

    mp_re : Floating_Number;
    mp_im : Floating_Number;
    gamma : Complex_Number;

  begin
    if fixed_gamma then
      mp_re := create(0.57670012968461137);
      mp_im := create(0.8169559109411918);
    else
      declare
        st_gamma : constant Standard_Complex_Numbers.Complex_Number
                 := Standard_Random_Numbers.Random1;
        st_re : constant double_float
              := Standard_Complex_Numbers.REAL_PART(st_gamma);
        st_im : constant double_float
              := Standard_Complex_Numbers.IMAG_PART(st_gamma);
      begin
        mp_re := create(st_re);
        mp_im := create(st_im);
      end;
    end if;
    gamma := Create(mp_re,mp_im);
    Init(p,q,s,gamma,2,deci);
    Clear(mp_re); Clear(mp_im);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_Sys;
                   gamma : in Complex_Number; k,deci : in natural32 ) is
  begin
    Init(p,q,gamma,k,deci,0);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_Sys; s : in Link_to_Solution;
                   gamma : in Complex_Number; k,deci : in natural32 ) is
  begin
    Init(p,q,s,gamma,k,deci,0);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_Sys;
                   gamma : in Complex_Number; k,deci,cp : in natural32 ) is
  begin
    Multprec_Homotopy.Clear;
    Multprec_Homotopy.Create(p.all,q.all,k,gamma);
    Continuation_Parameters.Tune(cp,deci);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_Sys; s : in Link_to_Solution;
                   gamma : in Complex_Number; k,deci,cp : in natural32 ) is
  begin
    Init(p,q,gamma,k,deci,cp);
    Init(s);
  end Init;

  procedure Init ( h : in Link_to_Poly_Sys; txk : in integer32;
                   cp,deci : in natural32 ) is
  begin
    Multprec_Homotopy.Clear;
    Multprec_Homotopy.Create(h.all,txk);
    Continuation_Parameters.Tune(cp,deci);
  end Init;

  procedure Init ( h : in Link_to_Poly_Sys; txk : in integer32;
                   cp,deci : in natural32; s : in Link_to_Solution ) is
  begin
    Multprec_Homotopy.Clear;
    Multprec_Homotopy.Create(h.all,txk);
    Continuation_Parameters.Tune(cp,deci);
    Init(s);
  end Init;

-- PREDICTOR-CORRECTOR STAGE :

  procedure Predictor_Corrector_Stage
              ( target : in Complex_Number;
                p : in Multprec_Continuation_Data.Pred_Pars;
                c : in Multprec_Continuation_Data.Corr_Pars ) is

  -- DESCRIPTION :
  --   Runs one stage of a predictor-corrector method towards the target,
  --   for the given predictor and corrector parameters.

    procedure Predictor is 
      new Single_Predictor
            (Max_Norm,Multprec_Homotopy.diff,Multprec_Homotopy.diff);
    procedure Affine_Corrector is
      new Affine_Single_Severe_Normal_Silent_Corrector
            (Max_Norm,Multprec_Homotopy.Eval,Multprec_Homotopy.diff);

  begin
    if p.predictor_type < 7 then
      Predictor
        (point,p,success,prev_sol.all,prev_v.all,vv.all,
         prev_t,target,step,tol,trial);
    elsif p.predictor_type = 7 then
      Single_Quadratic_Predictor
        (point,p,true,prev_sol.all,prev_sol0.all,
         prev_t,prev_t0,target,step,tol);
    elsif p.predictor_type = 8 then
      Single_Cubic_Predictor
        (point,p,true,prev_sol.all,prev_sol1.all,prev_sol0.all,
         prev_t,prev_t1,prev_t0,target,step,tol);
    else
      Single_Quartic_Predictor
        (point,p,true,prev_sol.all,prev_sol2.all,prev_sol1.all,prev_sol0.all,
         prev_t,prev_t2,prev_t1,prev_t0,target,step,tol);
    end if;
    Affine_Corrector(point,c);
    if p.predictor_type < 7 then
      Linear_Single_Management
        (point,p,c,old_t,prev_t,old_sol.all,prev_sol.all,
         old_v.all,prev_v.all,vv.all,step,nsuccess,trial,success);
    elsif p.predictor_type = 7 then
      Linear_Single_Quadratic_Management
        (point,p,c,old_t,prev_t,prev_t0,old_sol.all,prev_sol.all,
         prev_sol0.all,step,nsuccess,trial,success);
    elsif p.predictor_type = 8 then
      Linear_Single_Cubic_Management
        (point,p,c,old_t,prev_t,prev_t1,prev_t0,old_sol.all,prev_sol.all,
         prev_sol1.all,prev_sol0.all,step,nsuccess,trial,success);
    else
      Linear_Single_Quartic_Management
        (point,p,c,old_t,prev_t,prev_t2,prev_t1,prev_t0,
         old_sol.all,prev_sol.all,prev_sol2.all,prev_sol1.all,prev_sol0.all,
         step,nsuccess,trial,success);
    end if;
  end Predictor_Corrector_Stage;

  procedure Track_along_Path ( target : in Complex_Number ) is

  -- DESCRIPTION :
  --   Runs one stage of a predictor-corrector method,
  --   along a solution path, at a safe distance of the target.

    p : constant Continuation_Parameters.Pred_Pars
      := Continuation_Parameters.Create_for_Path;
    patpp : constant Multprec_Continuation_Data.Pred_Pars
          := Multprec_Continuation_Data.Convert(p);
    c : constant Continuation_Parameters.Corr_Pars
      := Continuation_Parameters.Create_for_Path;
    patcp : constant Multprec_Continuation_Data.Corr_Pars
          := Multprec_Continuation_Data.Convert(c);

  begin
    Predictor_Corrector_Stage(target,patpp,patcp);
  end Track_along_Path;

  procedure End_Game ( target : in Complex_Number ) is

  -- DESCRIPTION :
  --   Runs one stage of a predictor-corrector method
  --   towards the end of a solution path.

    p : constant Continuation_Parameters.Pred_Pars
      := Continuation_Parameters.Create_End_Game;
    endpp : constant Multprec_Continuation_Data.Pred_Pars
          := Multprec_Continuation_Data.Convert(p);
    c : constant Continuation_Parameters.Corr_Pars
      := Continuation_Parameters.Create_End_Game;
    endcp : constant Multprec_Continuation_Data.Corr_Pars
          := Multprec_Continuation_Data.Convert(c);

  begin
    Predictor_Corrector_Stage(target,endpp,endcp);
  end End_Game;

-- SELECTORS :

  function get_current return Link_to_Solution is
  begin
    return current;
  end get_current;

  function get_current return Solu_Info is
  begin
    return point;
  end get_current;

  function get_next return Link_to_Solution is

    res : Link_to_Solution;
    one : Floating_Number := create(1.0);
    target_t : Complex_Number := Create(one);

  begin
    res := get_next(target_t);
    Clear(one); Clear(target_t);
    return res;
  end get_next;

  function get_next ( target_t : Complex_Number ) return Link_to_Solution is

    dist : constant double_float := Continuation_Parameters.start_end_game;
    diff : Floating_Number := AbsVal(current.t - target_t);
    abs_diff : constant double_float := Round(diff);

  begin
    if abs_diff <= dist
     then End_Game(target_t);
     else Track_along_Path(target_t);
    end if;
    current := Shallow_Create(point);
    Clear(diff);
    return current;
  end get_next;

-- DESTRUCTOR :

  procedure Clear is
  begin
    Multprec_Homotopy.Clear;
    current := null;
    Multprec_Complex_Vectors.Clear(old_sol);
    Multprec_Complex_Vectors.Clear(prev_sol);
    Multprec_Complex_Vectors.Clear(prev_sol0);
    Multprec_Complex_Vectors.Clear(prev_sol1);
    Multprec_Complex_Vectors.Clear(prev_sol2);
    Multprec_Complex_Vectors.Clear(old_v);
    Multprec_Complex_Vectors.Clear(prev_v);
    Multprec_Complex_Vectors.Clear(vv);
  end Clear;

begin
  current := null;
  old_sol := null;
  prev_sol := null;
  prev_sol0 := null;
  prev_sol1 := null;
  prev_sol2 := null;
  old_v := null;
  prev_v := null;
  vv := null;
end Multprec_Path_Tracker;
