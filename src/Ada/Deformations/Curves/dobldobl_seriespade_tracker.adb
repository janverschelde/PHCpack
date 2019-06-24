with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_Homotopy;
with DoblDobl_Pade_Approximants_io;
with Homotopy_Pade_Approximants;
with Series_and_Homotopies;
with Series_and_Predictors;
with Series_and_Trackers;
with Standard_Pade_Trackers;
with Homotopy_Newton_Steps;
with Homotopy_Mixed_Residuals;

package body DoblDobl_SeriesPade_Tracker is

-- INTERNAL DATA :

  nbeqs : integer32; -- number of equations
  nbvar : integer32; -- number of variables
  idxpar : integer32; -- index of the continuation parameter, 0 if artificial
  homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters;
  htp : DoblDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
  abh : DoblDobl_Complex_Poly_SysFun.Link_to_Eval_Poly_Sys;
  current_poles : DoblDobl_Complex_VecVecs.Link_to_VecVec;
  current : Link_to_Solution;
  current_servec : DoblDobl_Complex_Series_Vectors.Link_to_Vector;
  current_padvec : DoblDobl_Pade_Approximants.Link_to_Pade_Vector;
  current_frp : double_double;
  current_cfp : DoblDobl_Complex_Numbers.Complex_Number;
  current_step : double_float;

-- CONSTRUCTORS :

  procedure Init ( pars : in Homotopy_Continuation_Parameters.Parameters ) is
  begin
    homconpars := new Homotopy_Continuation_Parameters.Parameters'(pars);
  end Init;

  procedure Initialize_Series_and_Approximants is

  -- DESCRIPTION :
  --   Allocates space for power series and Pade approximants.

    numdeg : constant integer32 := integer32(homconpars.numdeg);
    dendeg : constant integer32 := integer32(homconpars.dendeg);
    servec : DoblDobl_Complex_Series_Vectors.Vector(1..nbvar);
    padvec : constant DoblDobl_Pade_Approximants.Pade_Vector
           := DoblDobl_Pade_Approximants.Allocate(nbvar,numdeg,dendeg);

    use Homotopy_Pade_Approximants;
    allpoles : constant DoblDobl_Complex_VecVecs.VecVec
             := Allocate_DoblDobl_Poles(nbeqs,dendeg);

  begin
    current_servec := new DoblDobl_Complex_Series_Vectors.Vector'(servec);
    current_padvec := new DoblDobl_Pade_Approximants.Pade_Vector'(padvec);
    current_poles := new DoblDobl_Complex_VecVecs.VecVec'(allpoles);
  end Initialize_Series_and_Approximants;

  procedure Init ( p,q : in Link_to_Poly_Sys ) is

    tpow : constant natural32 := 2;
    d_gamma : constant Standard_Complex_Numbers.Complex_Number
            := homconpars.gamma;
    dd_gamma : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(d_gamma);

  begin
    idxpar := 0;
    DoblDobl_Homotopy.Create(p.all,q.all,tpow,dd_gamma);
    abh := new DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys'
                 (Homotopy_Mixed_Residuals.DoblDobl_AbsVal_Homotopy);
    nbeqs := p'last;
    nbvar := integer32(Number_of_Unknowns(p(p'first)));
   -- definition of series homotopy is done with Init of the solution
    Initialize_Series_and_Approximants;
  end Init;

  procedure Init ( h : in Link_to_Poly_Sys; idx : in integer32 ) is
  begin
    idxpar := idx;
    DoblDobl_Homotopy.Create(h.all,idx);
    abh := new DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys'
                 (Homotopy_Mixed_Residuals.DoblDobl_AbsVal_Homotopy);
    nbeqs := h'last;
    nbvar := integer32(Number_of_Unknowns(h(h'first))) - 1;
    Initialize_Series_and_Approximants;
  end Init;

  procedure Init ( s : in Link_to_Solution ) is

    conpar : integer32;

  begin
    if idxpar = 0
     then conpar := nbeqs + 1;
     else conpar := idxpar;
    end if;
    current := s;
    DoblDobl_CSeries_Poly_Systems.Clear(htp);
    declare -- reset the shifted homotopy
      hs : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nbeqs)
         := DoblDobl_Homotopy.Homotopy_System;
      sh : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nbeqs)
         := Series_and_Homotopies.Create(hs,conpar,false);
    begin
      htp := new DoblDobl_CSeries_Poly_Systems.Poly_Sys'(sh);
    end;
  end Init;

-- PREDICTOR-CORRECTOR STAGE :

  procedure Predict ( fail : out boolean; verbose : in boolean := false ) is

    numdeg : constant integer32 := integer32(homconpars.numdeg);
    dendeg : constant integer32 := integer32(homconpars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2; -- + 1; -- + 2;
    nit : constant integer32 := integer32(homconpars.corsteps+2);
    sol : DoblDobl_Complex_Vectors.Vector(1..current.n) := current.v;
    eva : DoblDobl_Complex_Series_Vectors.Vector(1..nbeqs);
    t,predres : double_float;
    dd_t,dd_step : double_double;
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
      (current_servec.all,current_padvec.all,current_poles.all,
       current_frp,current_cfp);
    if verbose then
      put_line("The Pade vector : ");
      for i in current_padvec'range loop
        put_line(DoblDobl_Pade_Approximants_io.Write(current_padvec(i)));
      end loop;
      put_line("The poles : "); put_line(current_poles.all);
      put_line("The poles :"); put_line(current_poles.all);
      put("Smallest pole radius : "); put(current_frp,3); new_line;
      if DoblDobl_Complex_Numbers.REAL_PART(current_cfp) >= 0.0
       then put("Closest pole : "); put(current_cfp); new_line;
      end if;
      current_step := Series_and_Predictors.Set_Step_Size
                        (standard_output,eva,tolcff,alpha,verbose);
    else
      current_step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
    end if;
    current_step := homconpars.sbeta*current_step;
    DoblDobl_Complex_Series_Vectors.Clear(eva);
    if current_frp > 0.0 then
      current_step := Series_and_Predictors.Cap_Step_Size
                        (current_step,hi_part(current_frp),homconpars.pbeta);
    end if;
    dd_t := DoblDobl_Complex_Numbers.REAL_PART(current.t);
    t := hi_part(dd_t);
    Standard_Pade_Trackers.Set_Step(t,current_step,homconpars.maxsize,1.0);
    if verbose
     then put("Step size :"); put(current_step,2); put("  t ="); put(t,2);
    end if;
    loop
      dd_step := Double_Double_Numbers.Create(current_step);
      sol := Series_and_Predictors.Predicted_Solution
               (current_padvec.all,dd_step);
      predres := Series_and_Trackers.Residual_Prediction(sol,t);
      if verbose
       then put("  residual :"); put(predres,2); new_line;
      end if;
      exit when (predres <= alpha);
      t := t - current_step; current_step := current_step/2.0;
      t := t + current_step;
      if verbose
       then put("Step size :"); put(current_step,2); put("  t ="); put(t,2);
      end if;
      exit when (current_step < homconpars.minsize);
    end loop;
    dd_t := Double_Double_Numbers.Create(t);
    dd_step := Double_Double_Numbers.Create(current_step);
    current.t := DoblDobl_Complex_Numbers.Create(dd_t);
    current.v := sol;
    if t = 1.0 
     then fail := false;
     else fail := (current_step < homconpars.minsize);
    end if;
    DoblDobl_Complex_Series_Vectors.Clear(eva);
    Series_and_Homotopies.Shift(htp.all,-dd_step);
  end Predict;

  procedure Correct ( fail : out boolean; verbose : in boolean := false ) is

    t : constant double_float
      := hi_part(DoblDobl_Complex_Numbers.REAL_PART(current.t));
    nbrit : natural32 := 0;
    extra : constant natural32 := homconpars.corsteps;
    err,rco,res : double_float;

  begin
    if verbose then
      Homotopy_Newton_Steps.Correct
        (standard_output,abh.all,t,homconpars.tolres,homconpars.corsteps,nbrit,
         current.v,err,rco,res,fail,extra,true);
    else
      Homotopy_Newton_Steps.Correct
        (abh.all,t,homconpars.tolres,homconpars.corsteps,nbrit,
         current.v,err,rco,res,fail,extra);
    end if;
    current.err := Double_Double_Numbers.Create(err);
    current.rco := Double_Double_Numbers.Create(rco);
    current.res := Double_Double_Numbers.Create(res);
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
    return DoblDobl_Complex_Series_Vectors.Link_to_Vector is
  begin
    return current_servec;
  end Get_Current_Series_Vector;

  function Get_Current_Pade_Vector
    return DoblDobl_Pade_Approximants.Link_to_Pade_Vector is
  begin
    return current_padvec;
  end Get_Current_Pade_Vector;

  function Get_Current_Poles return DoblDobl_Complex_VecVecs.Link_to_VecVec is
  begin
    return current_poles;
  end Get_Current_Poles;

  function Get_Current_Pole_Radius return double_double is
  begin
    return current_frp;
  end Get_Current_Pole_Radius;

  function Get_Current_Closest_Pole
             return DoblDobl_Complex_Numbers.Complex_Number is
  begin
    return current_cfp;
  end Get_Current_Closest_Pole;

  function Get_Current_Step_Size return double_float is
  begin
    return current_step;
  end Get_Current_Step_Size;

  function Get_Current_t_Value return double_float is
  begin
    return hi_part(DoblDobl_Complex_Numbers.REAL_PART(current.t));
  end Get_Current_t_Value;

-- DESTRUCTOR :

  procedure Clear is
  begin
    Homotopy_Continuation_Parameters.Clear(homconpars);
    DoblDobl_CSeries_Poly_Systems.Clear(htp);
    DoblDobl_Complex_Poly_SysFun.Clear(abh);
    DoblDobl_Complex_VecVecs.Deep_Clear(current_poles);
    DoblDobl_Complex_Series_Vectors.Clear(current_servec);
    DoblDobl_Pade_Approximants.Clear(current_padvec);
  end Clear;

end DoblDobl_SeriesPade_Tracker;
