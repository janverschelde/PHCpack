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
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Double_Double_Vectors;
with Double_Double_Vectors_io;           use Double_Double_Vectors_io;
with Quad_Double_Vectors;
with Quad_Double_Vectors_io;             use Quad_Double_Vectors_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;
 use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Poly_Systems;
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
with Test_Pade_Predictors;

procedure ts_padepred is

-- DESCRIPTION :
--   Test on the application of Pade approximants to predict solutions.

  procedure Standard_Test_Prediction
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tests the predictor on the solution in sol for the homotopy hom,
  --   in standard double precision.

   -- neq : constant integer32 := hom'last;
    maxdeg,nit : integer32 := 4;
    srv : Standard_Complex_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Complex_Series_Vectors.Vector(hom'range);
    tolcff : constant double_float := 1.0e-12;
    tolres,step : double_float := 0.0;
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range);
    poles : Standard_Complex_VecVecs.VecVec(pv'range);
    degnum,degden : integer32 := 0;
    verbose : constant boolean := true;
    frp : double_float;

  begin
    new_line;
    put("Give the number of Newton iterations : "); get(nit);
    put("Give the maximal degree of the series : "); get(maxdeg);
   -- Test_Pade_Predictors.Standard_Test_Pade_Prediction(neq,sol,pv);
    if verbose then
      Series_and_Predictors.Newton_Prediction
        (standard_output,maxdeg,nit,hom,sol.v,srv,eva,true);
    else
      Series_and_Predictors.Newton_Prediction(maxdeg,nit,hom,sol.v,srv,eva);
    end if;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    pv := Standard_Pade_Approximants.Create(degnum,degden,srv,verbose);
    poles := Homotopy_Pade_Approximants.Standard_Poles(pv);
    frp := Homotopy_Pade_Approximants.Closest_Pole(poles);
    put("The closest pole : "); put(frp,3); new_line;
    Test_Pade_Predictors.Standard_Step_Prediction(hom,srv,eva,pv);
    new_line;
    put_line("Setting the step size based on the power series ...");
    put("Give the tolerance on the residual : "); get(tolres);
    step := Series_and_Predictors.Set_Step_Size
              (standard_output,eva,tolcff,tolres,true);
    put("The computed step size : "); put(step,3); new_line;
    Test_Series_Predictors.Standard_Check_Prediction(hom,srv,eva,step);
    Test_Pade_Predictors.Standard_Check_Prediction(hom,pv,step);
  end Standard_Test_Prediction;

  procedure DoblDobl_Test_Prediction
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tests the predictor on the solution in sol for the homotopy hom,
  --   in double double precision.

   -- neq : constant integer32 := hom'last;
    maxdeg,nit : integer32 := 4;
    srv : DoblDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : DoblDobl_Complex_Series_Vectors.Vector(hom'range);
    tolcff : constant double_float := 1.0e-12;
    tolres,step : double_float := 0.0;
    dd_step : double_double;
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range);
    degnum,degden : integer32 := 0;
    verbose : constant boolean := false;
    frp : double_double;

  begin
    new_line;
    put("Give the number of Newton iterations : "); get(nit);
    put("Give the maximal degree of the series : "); get(maxdeg);
   -- Test_Pade_Predictors.DoblDobl_Test_Pade_Prediction(neq,sol,pv);
    Series_and_Predictors.Newton_Prediction(maxdeg,nit,hom,sol.v,srv,eva);
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    pv := DoblDobl_Pade_Approximants.Create(degnum,degden,srv,verbose);
    poles := Homotopy_Pade_Approximants.DoblDobl_Poles(pv);
    frp := Homotopy_Pade_Approximants.Closest_Pole(poles);
    put("The smallest forward pole : "); put(frp,3); new_line;
    Test_Pade_Predictors.DoblDobl_Step_Prediction(hom,srv,eva,pv);
    new_line;
    put_line("Setting the step size based on the power series ...");
    put("Give the tolerance on the residual : "); get(tolres);
    step := Series_and_Predictors.Set_Step_Size
              (standard_output,eva,tolcff,tolres,true);
    put("The computed step size : "); put(step,3); new_line;
    dd_step := create(step);
    Test_Series_Predictors.DoblDobl_Check_Prediction(hom,srv,eva,dd_step);
    Test_Pade_Predictors.DoblDobl_Check_Prediction(hom,pv,dd_step);
  end DoblDobl_Test_Prediction;

  procedure QuadDobl_Test_Prediction
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tests the predictor on the solution in sol for the homotopy hom,
  --   in quad double precision.

   -- neq : constant integer32 := hom'last;
    maxdeg,nit : integer32 := 4;
    srv : QuadDobl_Complex_Series_Vectors.Vector(1..sol.n);
    eva : QuadDobl_Complex_Series_Vectors.Vector(hom'range);
    tolcff : constant double_float := 1.0e-12;
    tolres,step : double_float := 0.0;
    qd_step : quad_double;
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range);
    degnum,degden : integer32 := 0;
    verbose : constant boolean := false;
    frp : quad_double;

  begin
    new_line;
    put("Give the number of Newton iterations : "); get(nit);
    put("Give the maximal degree of the series : "); get(maxdeg);
   -- Test_Pade_Predictors.QuadDobl_Test_Pade_Prediction(neq,sol,pv);
    Series_and_Predictors.Newton_Prediction(maxdeg,nit,hom,sol.v,srv,eva);
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    pv := QuadDobl_Pade_Approximants.Create(degnum,degden,srv,verbose);
    poles := Homotopy_Pade_Approximants.QuadDobl_Poles(pv);
    frp := Homotopy_Pade_Approximants.Closest_Pole(poles);
    put("The smallest forward pole : "); put(frp,3); new_line;
    Test_Pade_Predictors.QuadDobl_Step_Prediction(hom,srv,eva,pv);
    new_line;
    put_line("Setting the step size based on the power series ...");
    put("Give the tolerance on the residual : "); get(tolres);
    step := Series_and_Predictors.Set_Step_Size
              (standard_output,eva,tolcff,tolres,true);
    put("The computed step size : "); put(step,3); new_line;
    qd_step := create(step);
    Test_Series_Predictors.QuadDobl_Check_Prediction(hom,srv,eva,qd_step);
    Test_Pade_Predictors.QuadDobl_Check_Prediction(hom,pv,qd_step);
  end QuadDobl_Test_Prediction;

  procedure Standard_Forward_Pole_Radius
              ( nq : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Computes the forward pole radius for the solutions in sols,
  --   for a homotopy with nq equations, in standard double precision.

  -- REQUIRED :
  --   Standard_Homotopy has been initialized.

    len : constant integer32
        := integer32(Standard_Complex_Solutions.Length_Of(sols));
    rad : Standard_Floating_Vectors.Vector(1..len);
    numdeg,dendeg : integer32 := 0;
    fpr : double_float;

  begin
    new_line;
    put("Give the degree of the numerator : "); get(numdeg);
    put("Give the degree of the denominator : "); get(dendeg);
    Test_Pade_Predictors.Forward_Pole_Radius(nq,numdeg,dendeg,sols,rad,fpr);
    put_line("The forward pole radii : "); put_line(rad);
    put("The forward pole radius : "); put(fpr,3); new_line;
  end Standard_Forward_Pole_Radius;

  procedure DoblDobl_Forward_Pole_Radius
              ( nq : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Computes the forward pole radius for the solutions in sols,
  --   for a homotopy with nq equations, in double double precision.

  -- REQUIRED :
  --   DoblDobl_Homotopy has been initialized.

    len : constant integer32
        := integer32(DoblDobl_Complex_Solutions.Length_Of(sols));
    rad : Double_Double_Vectors.Vector(1..len);
    numdeg,dendeg : integer32 := 0;
    fpr : double_double;

  begin
    new_line;
    put("Give the degree of the numerator : "); get(numdeg);
    put("Give the degree of the denominator : "); get(dendeg);
    Test_Pade_Predictors.Forward_Pole_Radius(nq,numdeg,dendeg,sols,rad,fpr);
    put_line("The forward pole radii : "); put_line(rad);
    put("The forward pole radius : "); put(fpr,3); new_line;
  end DoblDobl_Forward_Pole_Radius;

  procedure QuadDobl_Forward_Pole_Radius
              ( nq : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Computes the forward pole radius for the solutions in sols,
  --   for a homotopy with nq equations, in quad double precision.

  -- REQUIRED :
  --   QuadDobl_Homotopy has been initialized.

    len : constant integer32
        := integer32(QuadDobl_Complex_Solutions.Length_Of(sols));
    rad : Quad_Double_Vectors.Vector(1..len);
    numdeg,dendeg : integer32 := 0;
    fpr : quad_double;

  begin
    new_line;
    put("Give the degree of the numerator : "); get(numdeg);
    put("Give the degree of the denominator : "); get(dendeg);
    Test_Pade_Predictors.Forward_Pole_Radius(nq,numdeg,dendeg,sols,rad,fpr);
    put_line("The forward pole radii : "); put_line(rad);
    put("The forward pole radius : "); put(fpr,3); new_line;
  end QuadDobl_Forward_Pole_Radius;

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
    ans : character;

  begin
    put_line("The homotopy system : "); put(h);
    put("Compute the global forward pole radius ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Forward_Pole_Radius(nq,sols);
    end if;
    put("Running the predictor on ");
    put(len,1); put_line(" solutions ...");
    for k in 1..len loop
      put("At solution "); put(k,1); put_line("...");
      sol := Standard_Complex_Solutions.Head_Of(tmp);
      Standard_Test_Prediction(s,sol.all);
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
    ans : character;

  begin
    put("Compute the global forward pole radius ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then DoblDobl_Forward_Pole_Radius(nq,sols);
    end if;
    put("Running the predictor on ");
    put(len,1); put_line(" solutions ...");
    for k in 1..len loop
      put("At solution "); put(k,1); put_line("...");
      sol := DoblDobl_Complex_Solutions.Head_Of(tmp);
      DoblDobl_Test_Prediction(s,sol.all);
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
    ans : character;

  begin
    put("Compute the global forward pole radius ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then QuadDobl_Forward_Pole_Radius(nq,sols);
    end if;
    put("Running the predictor on ");
    put(len,1); put_line(" solutions ...");
    for k in 1..len loop
      put("At solution "); put(k,1); put_line("...");
      sol := QuadDobl_Complex_Solutions.Head_Of(tmp);
      QuadDobl_Test_Prediction(s,sol.all);
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
  --   Prompts for the working precision and launches the test.

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
end ts_padepred;
