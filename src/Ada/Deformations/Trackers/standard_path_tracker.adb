-- for testing purposes :
--with text_io;
--with Standard_Complex_Numbers_io;
-- use Standard_Complex_Numbers_io;
--with Standard_Complex_Vectors_io;
-- use Standard_Complex_Vectors_io;

with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Homotopy;
with Standard_Dispatch_Predictors;       use Standard_Dispatch_Predictors;
with Standard_Correctors;                use Standard_Correctors;
with Standard_Data_on_Path;              use Standard_Data_on_Path;
with Continuation_Parameters;

package body Standard_Path_Tracker is

-- INTERNAL DATA :

  current : Link_to_Solution;
  point : Solu_Info;
  tol : constant double_float := 1.0E-12;
  old_t,prev_t,prev_t0,prev_t1,prev_t2 : Complex_Number;
  old_sol,prev_sol,prev_sol0,prev_sol1,prev_sol2,old_v,prev_v,vv
    : Standard_Complex_Vectors.Link_to_Vector;
  step : double_float;
  nsuccess,trial : natural32 := 0;
  success : boolean := true;

-- CONSTRUCTORS :

  procedure Clear_Solution_Data is

  -- DESCRIPTION :
  --   For the initialization of a second solution with the same homotopy,
  --   the intermediate solution data must be cleared first.

    use Standard_Complex_Vectors;

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
    old_sol := new Standard_Complex_Vectors.Vector(current.v'range);
    prev_sol := new Standard_Complex_Vectors.Vector(current.v'range);
    prev_sol0 := new Standard_Complex_Vectors.Vector(current.v'range);
    prev_sol1 := new Standard_Complex_Vectors.Vector(current.v'range);
    prev_sol2 := new Standard_Complex_Vectors.Vector(current.v'range);
    old_v := new Standard_Complex_Vectors.Vector(current.v'range);
    prev_v := new Standard_Complex_Vectors.Vector(current.v'range);
    vv := new Standard_Complex_Vectors.Vector(current.v'range);
    old_t := current.t;
    old_sol.all := current.v;
    prev_t  := old_t;  prev_sol.all  := old_sol.all;
    prev_t0 := prev_t; prev_sol0.all := prev_sol.all;
    prev_t1 := prev_t; prev_sol1.all := prev_sol.all;
    prev_t2 := prev_t; prev_sol2.all := prev_sol.all;
    prev_v.all := (prev_v'range => Create(0.0));
  end Init_Solution_Data;

  procedure Init ( s : in Link_to_Solution ) is
  begin
    current := s;
    point := Shallow_Create(current);
    step := Continuation_Parameters.max_path_step_size;
    nsuccess := 0;
    trial := 0;
    success := true;
    Clear_Solution_Data;
    Init_Solution_Data;
  end Init;

  procedure Init ( p,q : in Link_to_Poly_sys; fixed_gamma : in boolean ) is

    gamma : Complex_Number;

  begin
    if fixed_gamma
     then gamma := Create(0.57670012968461137, 0.8169559109411918);
     else gamma := Standard_Random_Numbers.Random1;
    end if;
    Init(p,q,gamma,2);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_sys; fixed_gamma : in boolean;
                   s : in Link_to_Solution ) is

    gamma : Complex_Number;

  begin
    if fixed_gamma
     then gamma := Create(0.57670012968461137, 0.8169559109411918);
     else gamma := Standard_Random_Numbers.Random1;
    end if;
    Init(p,q,s,gamma,2);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_sys;
                   gamma : in Complex_Number; k : in natural32 ) is
  begin
    Init(p,q,gamma,k,0);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_sys; s : in Link_to_Solution;
                   gamma : in Complex_Number; k : in natural32 ) is
  begin
    Init(p,q,s,gamma,k,0);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_sys;
                   gamma : in Complex_Number; k,cp : in natural32 ) is
  begin
    Standard_Homotopy.Clear;
    Standard_Homotopy.Create(p.all,q.all,k,gamma);
    Continuation_Parameters.Tune(cp);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_sys; s : in Link_to_Solution;
                   gamma : in Complex_Number; k,cp : in natural32 ) is
  begin
    Init(p,q,gamma,k,cp);
    Init(s);
  end Init;

  procedure Init ( h : in Link_to_Poly_Sys; txk : in integer32 ) is
  begin
    Standard_Homotopy.Clear;
    Standard_Homotopy.Create(h.all,txk);
  end Init;

  procedure Init ( h : in Link_to_Poly_Sys; txk : in integer32;
                   s : in Link_to_Solution ) is
  begin
    Standard_Homotopy.Clear;
    Standard_Homotopy.Create(h.all,txk);
    Init(s);
  end Init;

-- PREDICTOR-CORRECTOR STAGE :

  procedure Predictor_Corrector_Stage
              ( target : in Complex_Number;
                p : in Continuation_Parameters.Pred_Pars;
                c : in Continuation_Parameters.Corr_Pars ) is

  -- DESCRIPTION :
  --   Runs one stage of a predictor-corrector method towards the target,
  --   for the given predictor and corrector parameters.

    procedure Predictor is 
      new Single_Predictor
            (Max_Norm,Standard_Homotopy.diff,Standard_Homotopy.diff);
   -- procedure Affine_Corrector is
   --   new Affine_Single_Severe_Normal_Reporting_Corrector
   --         (Max_Norm,Standard_Homotopy.Eval,Standard_Homotopy.diff);
    procedure Affine_Corrector is
      new Affine_Single_Severe_Normal_Silent_Corrector
            (Max_Norm,Standard_Homotopy.Eval,Standard_Homotopy.diff);

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
   -- text_io.put("The predicted t value : "); put(point.sol.t); 
   -- text_io.new_line;
   -- text_io.put_line("The predicted solution vector :");
   -- put_line(point.sol.v);
   -- Affine_Corrector(text_io.standard_output,point,c);
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
    c : constant Continuation_Parameters.Corr_Pars
      := Continuation_Parameters.Create_for_Path;

  begin
    Predictor_Corrector_Stage(target,p,c);
  end Track_along_Path;

  procedure End_Game ( target : in Complex_Number ) is

  -- DESCRIPTION :
  --   Runs one stage of a predictor-corrector method
  --   towards the end of a solution path.

    p : constant Continuation_Parameters.Pred_Pars
      := Continuation_Parameters.Create_End_Game;
    c : constant Continuation_Parameters.Corr_Pars
      := Continuation_Parameters.Create_End_Game;

  begin
    Predictor_Corrector_Stage(target,p,c);
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

    target_t : constant Complex_Number := Create(1.0);

  begin
    return get_next(target_t);
  end get_next;

  function get_next ( target_t : Complex_Number ) return Link_to_Solution is

    dist : constant double_float := Continuation_Parameters.start_end_game;

  begin
    if AbsVal(current.t - target_t) <= dist
     then End_Game(target_t);
     else Track_along_Path(target_t);
    end if;
    current := Shallow_Create(point);
    return current;
  end get_next;

-- DESTRUCTOR :

  procedure Clear is
  begin
    Standard_Homotopy.Clear;
    current := null;
    Standard_Complex_Vectors.Clear(old_sol);
    Standard_Complex_Vectors.Clear(prev_sol);
    Standard_Complex_Vectors.Clear(prev_sol0);
    Standard_Complex_Vectors.Clear(prev_sol1);
    Standard_Complex_Vectors.Clear(prev_sol2);
    Standard_Complex_Vectors.Clear(old_v);
    Standard_Complex_Vectors.Clear(prev_v);
    Standard_Complex_Vectors.Clear(vv);
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
end Standard_Path_Tracker;
