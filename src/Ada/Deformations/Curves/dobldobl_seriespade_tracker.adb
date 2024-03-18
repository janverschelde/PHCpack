with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with Symbol_Table;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Hessians;
with DoblDobl_Complex_Series_VecVecs;
with DoblDobl_CSeries_Vector_Functions;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_SysFun;
with DoblDobl_CSeries_Jaco_Matrices;
with DoblDobl_Homotopy;
with DoblDobl_Coefficient_Homotopy;
with Projective_Transformations;
with DoblDobl_Pade_Approximants_io;
with Homotopy_Pade_Approximants;
with Series_and_Homotopies;
with Series_and_Predictors;
with Standard_Pade_Trackers;
with DoblDobl_Pade_Trackers;
with Homotopy_Series_Readers;
with Homotopy_Newton_Steps;
with Homotopy_Mixed_Residuals;
with Homotopy_Coefficient_Scaling;
with Singular_Values_of_Hessians;

package body DoblDobl_SeriesPade_Tracker is

-- INTERNAL DATA :

  nbeqs : integer32; -- number of equations
  nbvar : integer32; -- number of variables
  idxpar : integer32; -- index of the continuation parameter, 0 if artificial
  homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters;
  htp : DoblDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
  fhm : DoblDobl_CSeries_Poly_SysFun.Link_to_Eval_Coeff_Poly_Sys;
  fcf : DoblDobl_Complex_Series_VecVecs.Link_to_VecVec;
  ejm : DoblDobl_CSeries_Jaco_Matrices.Link_to_Eval_Coeff_Jaco_Mat;
  mlt : DoblDobl_CSeries_Jaco_Matrices.Link_to_Mult_Factors;
  abh : DoblDobl_Complex_Poly_SysFun.Link_to_Eval_Poly_Sys;
  jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
  hs : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
  current_poles : DoblDobl_Complex_VecVecs.Link_to_VecVec;
  predicted : Link_to_Solution; -- predicted solution
  current : Link_to_Solution;   -- current solution
  current_servec : DoblDobl_Complex_Series_Vectors.Link_to_Vector;
  current_padvec : DoblDobl_Pade_Approximants.Link_to_Pade_Vector;
  current_frp,solnrm,eta : double_double;
  current_cfp : DoblDobl_Complex_Numbers.Complex_Number;
  series_step,pole_step,hessian_step,current_step : double_float;
  cntsstp,cntdstp,cntpstp : natural32;
  homcoord : boolean := false; -- homogeneous coordinates or not

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

  procedure Init ( p,q : in Link_to_Poly_Sys; homogeneous : in boolean ) is

    tpow : constant natural32 := 1; -- := 2; may introduce singularities
    d_gamma : constant Standard_Complex_Numbers.Complex_Number
            := homconpars.gamma;
    dd_gamma : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(d_gamma);
    vp : Link_to_Poly_Sys := p;
    vq : Link_to_Poly_Sys := q;

    use Singular_Values_of_Hessians;

  begin
    idxpar := 0;
    homcoord := homogeneous;
    if not homogeneous then
      DoblDobl_Homotopy.Create(p.all,q.all,tpow,dd_gamma);
    else
      Homotopy_Series_Readers.DoblDobl_Projective_Transformation(vp,vq);
      Symbol_Table.Enlarge(1);
      Symbol_Table.Add_String("Z0");
      DoblDobl_Homotopy.Create(vp.all,vq.all,1,dd_gamma);
      DoblDobl_Coefficient_Homotopy.Create(vq.all,vp.all,1,dd_gamma);
    end if;
    abh := new DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys'
                 (Homotopy_Mixed_Residuals.DoblDobl_AbsVal_Homotopy);
    if homogeneous then
      nbeqs := vp'last;
      nbvar := integer32(Number_of_Unknowns(vp(vp'first)));
    else
      nbeqs := p'last;
      nbvar := integer32(Number_of_Unknowns(p(p'first)));
    end if;
   -- definition of series homotopy is done with Init of the solution
    Initialize_Series_and_Approximants;
    DoblDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
  end Init;

  procedure Init ( h : in Link_to_Poly_Sys; idx : in integer32 ) is

    use Singular_Values_of_Hessians;

  begin
    idxpar := idx;
    DoblDobl_Homotopy.Create(h.all,idx);
    abh := new DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys'
                 (Homotopy_Mixed_Residuals.DoblDobl_AbsVal_Homotopy);
    nbeqs := h'last;
    nbvar := integer32(Number_of_Unknowns(h(h'first))) - 1;
    Initialize_Series_and_Approximants;
    DoblDobl_Jacobian_Hessians_of_Homotopy(idx,jm,hs);
  end Init;

  procedure Init ( s : in Link_to_Solution ) is

    conpar : integer32;

  begin
    if idxpar = 0
     then conpar := nbeqs + 1;
     else conpar := idxpar;
    end if;
    if not homcoord then
      current := s;
      predicted := new Solution'(s.all);
    else
      declare
        p1s : constant Solution(s.n+1)
            := Projective_Transformations.Projective_Transformation(s.all);
      begin
        current := new Solution'(p1s);
        predicted := new Solution'(p1s);
      end;
    end if;
    DoblDobl_CSeries_Poly_Systems.Clear(htp);
    DoblDobl_CSeries_Poly_SysFun.Clear(fhm);
    DoblDobl_Complex_Series_VecVecs.Deep_Clear(fcf);
    DoblDobl_CSeries_Jaco_Matrices.Clear(ejm);
    DoblDobl_CSeries_Jaco_Matrices.Clear(mlt);
    declare -- reset the shifted homotopy
      hs : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nbeqs)
         := DoblDobl_Homotopy.Homotopy_System;
      sh : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nbeqs)
         := Series_and_Homotopies.Create(hs,conpar,false);
      fs : constant DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(sh'range)
         := DoblDobl_CSeries_Poly_SysFun.Create(sh);
      fc : constant DoblDobl_Complex_Series_VecVecs.VecVec(sh'range)
         := DoblDobl_CSeries_Poly_SysFun.Coeff(sh);
      nv : constant integer32 := nbvar;
      dm : DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat(sh'range,1..nv);
      mt : DoblDobl_CSeries_Jaco_Matrices.Mult_Factors(sh'range,1..nv);
    begin
      htp := new DoblDobl_CSeries_Poly_Systems.Poly_Sys'(sh);
      fhm := new DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys'(fs);
      fcf := new DoblDobl_Complex_Series_VecVecs.VecVec'(fc);
      DoblDobl_CSeries_Jaco_Matrices.Create(sh,dm,mt);
      ejm := new DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat'(dm);
      mlt := new DoblDobl_CSeries_Jaco_Matrices.Mult_Factors'(mt);
    end;
    cntsstp := 0; cntdstp := 0; cntpstp := 0;
  end Init;

-- PREDICTOR-CORRECTOR STAGE :

  procedure Step_Control ( verbose : in boolean := false ) is

    numdeg : constant integer32 := integer32(homconpars.numdeg);
    dendeg : constant integer32 := integer32(homconpars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32
        := Standard_Pade_Trackers.Maximum(5,maxdeg/2) + 1;
    eva : DoblDobl_Complex_Series_Vectors.Vector(1..nbeqs);
    dd_t : double_double := DoblDobl_Complex_Numbers.REAL_PART(current.t);
    t : double_float := hi_part(dd_t);
   -- tolcff : constant double_float := homconpars.epsilon;
   -- alpha : constant double_float := homconpars.alpha;

  begin
    if verbose then
     -- Series_and_Predictors.Newton_Prediction
     --   (standard_output,maxdeg,nit,htp.all,current.v,current_servec.all,
     --    eva,verbose);
      Series_and_Predictors.Newton_Prediction
        (standard_output,maxdeg,nit,fhm.all,fcf.all,ejm.all,mlt.all,
         current.v,current_servec.all,eva,true);
     -- series_step := Series_and_Predictors.Set_Step_Size
     --                  (standard_output,eva,tolcff,alpha,verbose);
    else
     -- Series_and_Predictors.Newton_Prediction
     --   (maxdeg,nit,htp.all,current.v,current_servec.all,eva);
      Series_and_Predictors.Newton_Prediction
        (maxdeg,nit,fhm.all,fcf.all,ejm.all,mlt.all,
         current.v,current_servec.all,eva);
     -- series_step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
    end if;
    series_step := 1.0; -- disable series step -- homconpars.sbeta*series_step;
   -- if verbose
   --  then put("series step : "); put(series_step,2); new_line;
   -- end if;
    Series_and_Predictors.Pade_Approximants
      (current_servec.all,current_padvec.all,current_poles.all,
       current_frp,current_cfp);
    pole_step := homconpars.pbeta*hi_part(current_frp);
    if verbose then
      put_line("The Pade vector : ");
      for i in current_padvec'range loop
        put_line(DoblDobl_Pade_Approximants_io.Write(current_padvec(i)));
      end loop;
      put_line("The poles : "); put_line(current_poles.all);
      put("pole step : "); put(pole_step,2);
      put("  smallest pole radius : "); put(current_frp,2); new_line;
      put("closest pole : "); put(current_cfp); new_line;
    end if;
    if not verbose then
      hessian_step := Series_and_Predictors.Step_Distance
                        (maxdeg,homconpars.cbeta,dd_t,jm,hs,current.v,
                         current_servec.all,current_padvec.all);
      Standard_Pade_Trackers.Minimum_Step_Size
        (series_step,hessian_step,pole_step,current_step,
         cntsstp,cntdstp,cntpstp);
    else
      declare -- must extend the solution vector with the value for t
        use DoblDobl_Complex_Vectors,Singular_Values_of_Hessians;
        solxt : Vector(current.v'first..current.v'last+1);
      begin
        solxt(current.v'range) := current.v;
        solxt(solxt'last) := current.t;
        eta := DoblDobl_Distance(jm.all,hs.all,solxt);
      end;
      solnrm := Homotopy_Pade_Approximants.Solution_Error_Norm
                  (current_servec.all,current_padvec.all);
      hessian_step := Series_and_Predictors.Step_Distance
                        (maxdeg,homconpars.cbeta,hi_part(eta),hi_part(solnrm));
      put("Hessian step : "); put(hessian_step,2);
      put("  eta : "); put(eta,2);
      put("  nrm : "); put(solnrm,2); new_line;
      Standard_Pade_Trackers.Minimum_Step_Size
        (standard_output,series_step,hessian_step,
         pole_step,current_step,cntsstp,cntdstp,cntpstp);
    end if;
    Standard_Pade_Trackers.Set_Step(t,current_step,homconpars.maxsize,1.0);
    dd_t := Double_Double_Numbers.Create(t);
    current.t := DoblDobl_Complex_Numbers.Create(dd_t);
    if verbose
     then put("Step size :"); put(current_step,2); put("  t ="); put(t,2);
    end if;
    DoblDobl_Complex_Series_Vectors.Clear(eva);
  end Step_Control;

  procedure Predictor_Feedback_Loop
              ( fail : out boolean; verbose : in boolean := false ) is

    sol : DoblDobl_Complex_Vectors.Vector(1..current.n) := current.v;
    t : double_float
      := hi_part(DoblDobl_Complex_Numbers.REAL_PART(current.t));
    predres : double_float;
    dd_step,dd_t : double_double;
    alpha : constant double_float := homconpars.alpha;

  begin
    fail := false;
    loop
      dd_step := Double_Double_Numbers.Create(current_step);
      sol := Series_and_Predictors.Predicted_Solution
               (current_padvec.all,dd_step);
      predres := DoblDobl_Pade_Trackers.Residual_Prediction(abh.all,sol,t);
      if verbose
       then put("  predictor residual :"); put(predres,2); new_line;
      end if;
      exit when (predres <= alpha);
      t := t - current_step; current_step := current_step/2.0;
      t := t + current_step;
      if verbose
       then put("Step size :"); put(current_step,2); put("  t ="); put(t,2);
      end if;
      if current_step < homconpars.minsize
       then fail := true; exit;
      end if;
    end loop;
    dd_t := Double_Double_Numbers.Create(t);
    predicted.t := DoblDobl_Complex_Numbers.Create(dd_t);
    predicted.v := sol;
    predicted.res := Double_Double_Numbers.Create(predres);
    current.t := DoblDobl_Complex_Numbers.Create(dd_t);
    current.v := sol;
  end Predictor_Feedback_Loop;

  procedure Predict ( fail : out boolean; verbose : in boolean := false ) is

    t : double_double := DoblDobl_Complex_Numbers.REAL_PART(current.t);
    dd_step : double_double;
    d_gamma : constant Standard_Complex_Numbers.Complex_Number
            := homconpars.gamma;
    dd_gamma : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(d_gamma);

  begin
   -- if verbose then
   --   Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
   --     (standard_output,fhm.all,fcf.all,current.v,t,dd_gamma,true);
   -- else -- verbose with scaling is for debugging purposes only
    Homotopy_Coefficient_Scaling.Scale_Solution_Coefficients
      (fcf.all,current.v,t,dd_gamma);
   -- end if;
    Step_Control(verbose);
    Predictor_Feedback_Loop(fail,verbose);
    t := DoblDobl_Complex_Numbers.REAL_PART(current.t);
    if hi_part(t) = 1.0 
     then fail := false;
     else fail := (current_step < homconpars.minsize);
    end if;
    dd_step := Double_Double_Numbers.Create(current_step);
    Series_and_Homotopies.Shift(htp.all,-dd_step);
  end Predict;

  procedure Correct ( fail : out boolean; verbose : in boolean := false ) is

    t : constant double_float
      := hi_part(DoblDobl_Complex_Numbers.REAL_PART(current.t));
    nbrit : natural32 := 0;
    extra : constant natural32 := 1;
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

    dd_step : double_double;

  begin
    Predict(fail,verbose);
    if not fail
     then Correct(fail,verbose);
    end if;
    dd_step := create(current_step);
    DoblDobl_CSeries_Vector_Functions.Shift(fcf.all,-dd_step);
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

  function Get_Predicted_Solution return Link_to_Solution is
  begin
    return predicted;
  end Get_Predicted_Solution;

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

  function Get_Current_Series_Step return double_float is
  begin
    return series_step;
  end Get_Current_Series_Step;

  function Get_Current_Pole_Step return double_float is
  begin
    return pole_step;
  end Get_Current_Pole_Step;

  function Get_Current_Estimated_Distance return double_double is
  begin
    return eta;
  end Get_Current_Estimated_Distance;

  function Get_Current_Hessian_Step return double_float is
  begin
    return hessian_step;
  end Get_Current_Hessian_Step;

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
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Hessians.Clear(hs);
    DoblDobl_Complex_VecVecs.Deep_Clear(current_poles);
    DoblDobl_Complex_Series_Vectors.Clear(current_servec);
    DoblDobl_Pade_Approximants.Clear(current_padvec);
    DoblDobl_Complex_Solutions.Clear(predicted);
  end Clear;

end DoblDobl_SeriesPade_Tracker;
