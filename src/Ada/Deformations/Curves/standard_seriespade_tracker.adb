with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_CSeries_Poly_Systems;
with Standard_Homotopy;
with Homotopy_Pade_Approximants;
with Series_and_Homotopies;
with Series_and_Predictors;
with Series_and_Trackers;
with Homotopy_Newton_Steps;

package body Standard_SeriesPade_Tracker is

-- INTERNAL DATA :

  nbeqs : integer32;
  homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters;
  htp : Standard_CSeries_Poly_Systems.Link_to_Poly_Sys;
  poles : Standard_Complex_VecVecs.Link_to_VecVec;
  current : Link_to_Solution;
  current_servec : Standard_Complex_Series_Vectors.Link_to_Vector;
  current_padvec : Standard_Pade_Approximants.Link_to_Pade_Vector;
  current_frp : double_float;
  current_cfp : Complex_Number;

-- CONSTRUCTORS :

  procedure Init ( pars : in Homotopy_Continuation_Parameters.Parameters ) is
  begin
    homconpars := new Homotopy_Continuation_Parameters.Parameters'(pars);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_Sys ) is

    tpow : constant natural32 := 2;
    gamma : constant Complex_Number := homconpars.gamma;
    nbvar : constant integer32
          := integer32(Number_of_Unknowns(p(p'first)));
    numdeg : constant integer32 := integer32(homconpars.numdeg);
    dendeg : constant integer32 := integer32(homconpars.dendeg);

  begin
    Standard_Homotopy.Create(p.all,q.all,tpow,gamma);
    nbeqs := p'last;
   -- setting of homotopy is done when initializing the solution
   -- declare
   --   h : Standard_Complex_Poly_Systems.Poly_Sys(1..nbeqs)
   --     := Standard_Homotopy.Homotopy_System;
   --   s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nbeqs)
   --     := Series_and_Homotopies.Create(h,nbeqs+1,false);
   -- begin
   --   htp := new Standard_CSeries_Poly_Systems.Poly_Sys'(s);
   -- end;
    declare
      servec : Standard_Complex_Series_Vectors.Vector(1..nbvar);
      padvec : constant Standard_Pade_Approximants.Pade_Vector
             := Standard_Pade_Approximants.Allocate(nbvar,numdeg,dendeg);
    begin
      current_servec := new Standard_Complex_Series_Vectors.Vector'(servec);
      current_padvec := new Standard_Pade_Approximants.Pade_Vector'(padvec);
    end;
    declare
      use Homotopy_Pade_Approximants;
      allpoles : constant Standard_Complex_VecVecs.VecVec
               := Allocate_Standard_Poles(nbeqs,dendeg);
    begin
      poles := new Standard_Complex_VecVecs.VecVec'(allpoles);
    end;
  end Init;

  procedure Init ( s : in Link_to_Solution ) is
  begin
    current := s;
    Standard_CSeries_Poly_Systems.Clear(htp);
    declare -- reset the shifted homotopy
      hs : Standard_Complex_Poly_Systems.Poly_Sys(1..nbeqs)
         := Standard_Homotopy.Homotopy_System;
      sh : Standard_CSeries_Poly_Systems.Poly_Sys(1..nbeqs)
         := Series_and_Homotopies.Create(hs,nbeqs+1,false);
    begin
      htp := new Standard_CSeries_Poly_Systems.Poly_Sys'(sh);
    end;
  end Init;

-- PREDICTOR-CORRECTOR STAGE :

  procedure Predict ( fail : out boolean; verbose : in boolean := false ) is

    numdeg : constant integer32 := integer32(homconpars.numdeg);
    dendeg : constant integer32 := integer32(homconpars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := integer32(homconpars.corsteps);
    sol : Standard_Complex_Vectors.Vector(1..current.n) := current.v;
    eva : Standard_Complex_Series_Vectors.Vector(1..nbeqs);
    t,step,predres : double_float := 0.0;
    tolcff : constant double_float := homconpars.epsilon;
    alpha : constant double_float := homconpars.alpha;

  begin
    if verbose then
      Series_and_Predictors.Newton_Prediction
        (standard_output,maxdeg,nit,htp.all,sol,current_servec.all,eva,verbose);
    else
      Series_and_Predictors.Newton_Prediction
        (maxdeg,nit,htp.all,sol,current_servec.all,eva);
    end if;
    Series_and_Predictors.Pade_Approximants
      (current_servec.all,current_padvec.all,poles.all,current_frp,current_cfp);
    if verbose then
      put_line("The poles : "); put_line(poles.all);
      put("Smallest forward pole radius : "); put(current_frp,2); new_line;
      if Standard_Complex_Numbers.REAL_PART(current_cfp) >= 0.0
       then put("Closest forward pole :"); put(current_cfp); new_line;
      end if;
      step := Series_and_Predictors.Set_Step_Size
                (standard_output,eva,tolcff,alpha,verbose);
    else
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
    end if;
    step := homconpars.sbeta*step;
    Standard_Complex_Series_Vectors.Clear(eva);
    if current_frp > 0.0 then
      step := Series_and_Predictors.Cap_Step_Size
                (step,current_frp,homconpars.pbeta);
    end if;
    t := Standard_Complex_Numbers.REAL_PART(current.t);
    Series_and_Trackers.Set_Step(t,step,homconpars.maxsize,1.0);
    if verbose
     then put("Step size :"); put(step,2); put("  t ="); put(t,2);
    end if;
    loop
      sol := Series_and_Predictors.Predicted_Solution
               (current_padvec.all,step);
      predres := Series_and_Trackers.Residual_Prediction(htp.all,sol,step);
      if verbose
       then put("  residual :"); put(predres,2); new_line;
      end if;
      exit when (predres <= alpha);
      t := t - step; step := step/2.0; t := t + step;
      if verbose
       then put("Step size :"); put(step,2); put("  t ="); put(t,2);
      end if;
      exit when (step < homconpars.minsize);
    end loop;
    current.t := Standard_Complex_Numbers.Create(t);
    current.v := sol;
    if t = 1.0 
     then fail := false;
     else fail := (step < homconpars.minsize);
    end if;
    Standard_Complex_Series_Vectors.Clear(eva);
    Series_and_Homotopies.Shift(htp.all,-step);
  end Predict;

  procedure Correct ( fail : out boolean; verbose : in boolean := false ) is

    t : constant double_float := REAL_PART(current.t);
    nbrit : natural32 := 0;

  begin
    if verbose then
      Homotopy_Newton_Steps.Correct
        (standard_output,nbeqs,t,homconpars.tolres,homconpars.corsteps,nbrit,
         current.v,current.err,current.rco,current.res,fail,true);
    else
      Homotopy_Newton_Steps.Correct
        (nbeqs,t,homconpars.tolres,homconpars.corsteps,nbrit,
         current.v,current.err,current.rco,current.res,fail);
    end if;
  end Correct;

  procedure Predict_and_Correct
              ( fail : out boolean; verbose : in boolean := false ) is
  begin
    Predict(fail,verbose);
    if not fail
     then Correct(fail,verbose);
    end if;
  end Predict_and_Correct;

-- SELECTORS :

  function Get_Parameters
    return Homotopy_Continuation_Parameters.Link_to_Parameters is

  begin
    return homconpars;
  end Get_Parameters;

  function Get_Current_Solution return Link_to_Solution is
  begin
    return current;
  end Get_Current_Solution;

  function Get_Current_Series_Vector
    return Standard_Complex_Series_Vectors.Link_to_Vector is
  begin
    return current_servec;
  end Get_Current_Series_Vector;

  function Get_Current_Pade_Vector
    return Standard_Pade_Approximants.Link_to_Pade_Vector is
  begin
    return current_padvec;
  end Get_Current_Pade_Vector;

  function Get_Current_Pole_Radius return double_float is
  begin
    return current_frp;
  end Get_Current_Pole_Radius;

  function Get_Current_Closest_Pole return Complex_Number is
  begin
    return current_cfp;
  end Get_Current_Closest_Pole;

-- DESTRUCTOR :

  procedure Clear is
  begin
    Homotopy_Continuation_Parameters.Clear(homconpars);
    Standard_CSeries_Poly_Systems.Clear(htp);
    Standard_Complex_VecVecs.Deep_Clear(poles);
    Standard_Complex_Series_Vectors.Clear(current_servec);
    Standard_Pade_Approximants.Clear(current_padvec);
  end Clear;

end Standard_SeriesPade_Tracker;
