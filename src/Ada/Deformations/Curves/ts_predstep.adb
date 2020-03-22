with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Hessians;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Hessians;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Hessians;
with QuadDobl_Complex_Solutions;
with Solution_Drops;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Complex_Series_Vectors;
with Standard_CSeries_Poly_Systems;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_CSeries_Poly_Systems;
with Series_and_Homotopies;
with Series_and_Predictors;
with Standard_Pade_Approximants;
with DoblDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants;
with Homotopy_Pade_Approximants;
with Test_Series_Predictors;
with Singular_Values_of_Hessians;

procedure ts_predstep is

-- DESCRIPTION :
--   Test on the three ways to predict the step size of a path tracker.

  procedure Standard_Test_Prediction
              ( jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tests the predictor on the one solution in sol for the homotopy hom,
  --   in standard double precision.

    nit : constant integer32 := 4;
    srv : Standard_Complex_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Complex_Series_Vectors.Vector(hom'range);
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range);
    poles : Standard_Complex_VecVecs.VecVec(pv'range);
    cpl : double_float;
    degnum,degden,maxdeg : integer32 := 0;
    verbose : constant boolean := false;
    eta,nrm,step1,step2 : double_float;
    beta1 : constant double_float := 0.5;
    beta2 : constant double_float := 0.005;

  begin
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    maxdeg := degnum + degden + 2;
    Series_and_Predictors.Newton_Prediction(maxdeg,nit,hom,sol.v,srv,eva);
    pv := Standard_Pade_Approximants.Create(degnum,degden,srv,verbose);
    poles := Homotopy_Pade_Approximants.Standard_Poles(pv);
    cpl := Homotopy_Pade_Approximants.Closest_Pole(poles);
    put("The smallest pole radius : "); put(cpl,3); new_line;
    step1 := beta1*cpl;
    put("step : "); put(step1,2); new_line;
    eta := Singular_Values_of_Hessians.Standard_Distance(jm.all,hs.all,sol);
    put(" eta : "); put(eta,2); new_line;
    nrm := Homotopy_Pade_Approximants.Solution_Error_Norm(srv,pv);
    put(" nrm : "); put(nrm,2); new_line;
    step2 := Series_and_Predictors.Step_Distance(maxdeg,beta2,eta,nrm);
    put("step : "); put(step2,2); new_line;
  end Standard_Test_Prediction;

  procedure DoblDobl_Test_Prediction
              ( jm : in DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tests the predictor on the one solution in sol for the homotopy hom,
  --   in double double precision.

    nit : constant integer32 := 4;
    srv : DoblDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : DoblDobl_Complex_Series_Vectors.Vector(hom'range);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range);
    cpl : double_double;
    maxdeg,degnum,degden : integer32 := 0;
    verbose : constant boolean := false;
    eta,nrm : double_double;
    step1,step2 : double_float;
    beta1 : constant double_float := 0.5;
    beta2 : constant double_float := 0.005;

  begin
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    maxdeg := degnum + degden + 2;
    Series_and_Predictors.Newton_Prediction(maxdeg,nit,hom,sol.v,srv,eva);
    pv := DoblDobl_Pade_Approximants.Create(degnum,degden,srv,verbose);
    poles := Homotopy_Pade_Approximants.DoblDobl_Poles(pv);
    cpl := Homotopy_Pade_Approximants.Closest_Pole(poles);
    put("The smallest pole radius : "); put(cpl,3); new_line;
    step1 := beta1*hi_part(cpl);
    put("step : "); put(step1,2); new_line;
    eta := Singular_Values_of_Hessians.DoblDobl_Distance(jm.all,hs.all,sol);
    put(" eta : "); put(eta,2); new_line;
    nrm := Homotopy_Pade_Approximants.Solution_Error_Norm(srv,pv);
    put(" nrm : "); put(nrm,2); new_line;
    step2 := Series_and_Predictors.Step_Distance
               (maxdeg,beta2,hi_part(eta),hi_part(nrm));
    put("step : "); put(step2,2); new_line;
  end DoblDobl_Test_Prediction;

  procedure QuadDobl_Test_Prediction
              ( jm : in QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tests the predictor on the one solution in sol for the homotopy hom,
  --   in quad double precision.

    nit : constant integer32 := 4;
    srv : QuadDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : QuadDobl_Complex_Series_Vectors.Vector(hom'range);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range);
    cpl : quad_double;
    maxdeg,degnum,degden : integer32 := 0;
    verbose : constant boolean := false;
    eta,nrm : quad_double;
    step1,step2 : double_float;
    beta1 : constant double_float := 0.5;
    beta2 : constant double_float := 0.005;

  begin
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    maxdeg := degnum + degden + 2;
    Series_and_Predictors.Newton_Prediction(maxdeg,nit,hom,sol.v,srv,eva);
    pv := QuadDobl_Pade_Approximants.Create(degnum,degden,srv,verbose);
    poles := Homotopy_Pade_Approximants.QuadDobl_Poles(pv);
    cpl := Homotopy_Pade_Approximants.Closest_Pole(poles);
    put("The smallest pole radius : "); put(cpl,3); new_line;
    step1 := beta1*hihi_part(cpl);
    put("step : "); put(step1,2); new_line;
    eta := Singular_Values_of_Hessians.QuadDobl_Distance(jm.all,hs.all,sol);
    put(" eta : "); put(eta,2); new_line;
    nrm := Homotopy_Pade_Approximants.Solution_Error_Norm(srv,pv);
    put(" nrm : "); put(nrm,2); new_line;
    step2 := Series_and_Predictors.Step_Distance
               (maxdeg,beta2,hihi_part(eta),hihi_part(nrm));
    put("step : "); put(step2,2); new_line;
  end QuadDobl_Test_Prediction;

  procedure Standard_Test_Prediction
              ( nq,idxpar : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The Standard_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    h : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : constant Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,idxpar);
    len : constant integer32
        := integer32(Standard_Complex_Solutions.Length_Of(sols));
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    sol : Standard_Complex_Solutions.Link_to_Solution;
    jm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : Standard_Complex_Hessians.Link_to_Array_of_Hessians;

    use Singular_Values_of_Hessians;

  begin
    if idxpar = 0
     then Standard_Jacobian_Hessians_of_Homotopy(jm,hs);
     else Standard_Jacobian_Hessians_of_Homotopy(idxpar,jm,hs);
    end if;
    put("Running the predictor on ");
    put(len,1); put_line(" solutions ...");
    for k in 1..len loop
      put("At solution "); put(k,1); put_line("...");
      sol := Standard_Complex_Solutions.Head_Of(tmp);
      Standard_Test_Prediction(jm,hs,s,sol.all);
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Standard_Test_Prediction;

  procedure DoblDobl_Test_Prediction
              ( nq,idxpar : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The DoblDobl_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    h : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,idxpar);
    len : constant integer32
        := integer32(DoblDobl_Complex_Solutions.Length_Of(sols));
    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    sol : DoblDobl_Complex_Solutions.Link_to_Solution;
    jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;

    use Singular_Values_of_Hessians;

  begin
    if idxpar = 0
     then DoblDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
     else DoblDobl_Jacobian_Hessians_of_Homotopy(idxpar,jm,hs);
    end if;
    put("Running the predictor on ");
    put(len,1); put_line(" solutions ...");
    for k in 1..len loop
      put("At solution "); put(k,1); put_line("...");
      sol := DoblDobl_Complex_Solutions.Head_Of(tmp);
      DoblDobl_Test_Prediction(jm,hs,s,sol.all);
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end DoblDobl_Test_Prediction;

  procedure QuadDobl_Test_Prediction
              ( nq,idxpar : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The QuadDobl_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    h : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,idxpar);
    len : constant integer32
        := integer32(QuadDobl_Complex_Solutions.Length_Of(sols));
    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    sol : QuadDobl_Complex_Solutions.Link_to_Solution;
    jm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hs : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;

    use Singular_Values_of_Hessians;

  begin
    if idxpar = 0
     then QuadDobl_Jacobian_Hessians_of_Homotopy(jm,hs);
     else QuadDobl_Jacobian_Hessians_of_Homotopy(idxpar,jm,hs);
    end if;
    put("Running the predictor on ");
    put(len,1); put_line(" solutions ...");
    for k in 1..len loop
      put("At solution "); put(k,1); put_line("...");
      sol := QuadDobl_Complex_Solutions.Head_Of(tmp);
      QuadDobl_Test_Prediction(jm,hs,s,sol.all);
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end QuadDobl_Test_Prediction;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in standard double precision.

    nbeq,idxpar : integer32;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.Standard_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      Standard_Test_Prediction(nbeq,nbeq+1,sols);
    else
      declare
        dropsols : constant Standard_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        Standard_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in double double precision.

    nbeq,idxpar : integer32;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.DoblDobl_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      DoblDobl_Test_Prediction(nbeq,nbeq+1,sols);
    else
      declare
        dropsols : constant DoblDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        DoblDobl_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in quad double precision.

    nbeq,idxpar : integer32;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.QuadDobl_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      QuadDobl_Test_Prediction(nbeq,nbeq+1,sols);
    else
      declare
        dropsols : constant QuadDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        QuadDobl_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the working precision and launches the test.

    ans : constant character := Prompt_for_Precision;

  begin
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_predstep;
