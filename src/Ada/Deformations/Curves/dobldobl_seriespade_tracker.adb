with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_Homotopy;
with Homotopy_Pade_Approximants;
with Series_and_Homotopies;
with Series_and_Predictors;
with Series_and_Trackers;
with Homotopy_Newton_Steps;

package body DoblDobl_SeriesPade_Tracker is

-- INTERNAL DATA :

  nbeqs : integer32;
  homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters;
  htp : DoblDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
  poles : DoblDobl_Complex_VecVecs.Link_to_VecVec;
  current : Link_to_Solution;
  current_servec : DoblDobl_Complex_Series_Vectors.Link_to_Vector;
  current_padvec : DoblDobl_Pade_Approximants.Link_to_Pade_Vector;

-- CONSTRUCTORS :

  procedure Init ( pars : in Homotopy_Continuation_Parameters.Parameters ) is
  begin
    homconpars := new Homotopy_Continuation_Parameters.Parameters'(pars);
  end Init;

  procedure Init ( p,q : in Link_to_Poly_Sys ) is

    tpow : constant natural32 := 2;
    d_gamma : constant Standard_Complex_Numbers.Complex_Number
            := homconpars.gamma;
    dd_gamma : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(d_gamma);
    nbvar : constant integer32
          := integer32(Number_of_Unknowns(p(p'first)));
    numdeg : constant integer32 := integer32(homconpars.numdeg);
    dendeg : constant integer32 := integer32(homconpars.dendeg);

  begin
    DoblDobl_Homotopy.Create(p.all,q.all,tpow,dd_gamma);
    nbeqs := p'last;
   -- declare
   --   h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nbeqs)
   --     := DoblDobl_Homotopy.Homotopy_System;
   --   s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nbeqs)
   --     := Series_and_Homotopies.Create(h,nbeqs+1,false);
   -- begin
   --   htp := new DoblDobl_CSeries_Poly_Systems.Poly_Sys'(s);
   -- end;
    declare
      servec : DoblDobl_Complex_Series_Vectors.Vector(1..nbvar);
      padvec : constant DoblDobl_Pade_Approximants.Pade_Vector
             := DoblDobl_Pade_Approximants.Allocate(nbvar,numdeg,dendeg);
    begin
      current_servec := new DoblDobl_Complex_Series_Vectors.Vector'(servec);
      current_padvec := new DoblDobl_Pade_Approximants.Pade_Vector'(padvec);
    end;
    declare
      use Homotopy_Pade_Approximants;
      allpoles : constant DoblDobl_Complex_VecVecs.VecVec
               := Allocate_DoblDobl_Poles(nbeqs,dendeg);
    begin
      poles := new DoblDobl_Complex_VecVecs.VecVec'(allpoles);
    end;
  end Init;

  procedure Init ( s : in Link_to_Solution ) is
  begin
    current := s;
    DoblDobl_CSeries_Poly_Systems.Clear(htp);
    declare -- reset the shifted homotopy
      hs : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nbeqs)
         := DoblDobl_Homotopy.Homotopy_System;
      sh : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nbeqs)
         := Series_and_Homotopies.Create(hs,nbeqs+1,false);
    begin
      htp := new DoblDobl_CSeries_Poly_Systems.Poly_Sys'(sh);
    end;
  end Init;

-- PREDICTOR-CORRECTOR STAGE :

  procedure Predict ( fail : out boolean; verbose : in boolean := false ) is

    numdeg : constant integer32 := integer32(homconpars.numdeg);
    dendeg : constant integer32 := integer32(homconpars.dendeg);
    maxdeg : constant integer32 := numdeg + dendeg + 2;
    nit : constant integer32 := integer32(homconpars.corsteps);
    sol : DoblDobl_Complex_Vectors.Vector(1..current.n) := current.v;
    eva : DoblDobl_Complex_Series_Vectors.Vector(1..nbeqs);
    t,step,predres : double_float;
    dd_t,dd_step,frp : double_double;
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
      (current_servec.all,current_padvec.all,poles.all,frp);
    if verbose then
      put("The smallest forward pole radius : "); put(frp,2); new_line;
      step := Series_and_Predictors.Set_Step_Size
               (standard_output,eva,tolcff,alpha,verbose);
    else
      step := Series_and_Predictors.Set_Step_Size(eva,tolcff,alpha);
    end if;
    step := homconpars.sbeta*step;
    DoblDobl_Complex_Series_Vectors.Clear(eva);
    if frp > 0.0 then
      step := Series_and_Predictors.Cap_Step_Size
                (step,hi_part(frp),homconpars.pbeta);
    end if;
    dd_t := DoblDobl_Complex_Numbers.REAL_PART(current.t);
    t := hi_part(dd_t);
    Series_and_Trackers.Set_Step(t,step,homconpars.maxsize,1.0);
    if verbose
     then put("Step size :"); put(step,2); put("  t ="); put(t,2);
    end if;
    loop
      dd_step := Double_Double_Numbers.Create(step);
      sol := Series_and_Predictors.Predicted_Solution
               (current_padvec.all,dd_step);
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
    dd_t := Double_Double_Numbers.Create(t);
    current.t := DoblDobl_Complex_Numbers.Create(dd_t);
    current.v := sol;
    if t = 1.0 
     then fail := false;
     else fail := (step < homconpars.minsize);
    end if;
    DoblDobl_Complex_Series_Vectors.Clear(eva);
    Series_and_Homotopies.Shift(htp.all,-dd_step);
  end Predict;

  procedure Correct ( fail : out boolean; verbose : in boolean := false ) is

    t : constant double_float
      := hi_part(DoblDobl_Complex_Numbers.REAL_PART(current.t));
    nbrit : natural32 := 0;
    err,rco,res : double_float;

  begin
    if verbose then
      Homotopy_Newton_Steps.Correct
        (standard_output,nbeqs,t,homconpars.tolres,homconpars.corsteps,nbrit,
         current.v,err,rco,res,fail,true);
    else
      Homotopy_Newton_Steps.Correct
        (nbeqs,t,homconpars.tolres,homconpars.corsteps,nbrit,
         current.v,err,rco,res,fail);
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

-- DESTRUCTOR :

  procedure Clear is
  begin
    Homotopy_Continuation_Parameters.Clear(homconpars);
    DoblDobl_CSeries_Poly_Systems.Clear(htp);
    DoblDobl_Complex_VecVecs.Deep_Clear(poles);
    DoblDobl_Complex_Series_Vectors.Clear(current_servec);
    DoblDobl_Pade_Approximants.Clear(current_padvec);
  end Clear;

end DoblDobl_SeriesPade_Tracker;
