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
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with Solution_Drops;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Dense_Series_Vectors;
with Standard_Series_Poly_Systems;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Series_Poly_Systems;
with Series_and_Homotopies;
with Series_and_Predictors;
with Standard_Pade_Approximants;
with DoblDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants;
with Test_Series_Predictors;
with Test_Pade_Predictors;

procedure ts_padepred is

-- DESCRIPTION :
--   Test on the application of Pade approximants to predict solutions.

  procedure Standard_Test_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tests the predictor on the solution in sol for the homotopy hom,
  --   in standard double precision.

    nit : integer32 := 4;
    srv : Standard_Dense_Series_Vectors.Vector(1..sol.n);
    eva : Standard_Dense_Series_Vectors.Vector(hom'range);
    tolcff : constant double_float := 1.0e-12;
    tolres,step : double_float := 0.0;
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range);

  begin
    new_line;
    put("Give the number of Newton iterations : "); get(nit);
    Test_Pade_Predictors.Standard_Test_Pade_Prediction(hom,sol,nit,pv);
    Series_and_Predictors.Newton_Prediction(nit,hom,sol.v,srv,eva);
    Test_Pade_Predictors.Standard_Step_Prediction(hom,srv,eva,pv);
    new_line;
    put_line("Setting the step size based on the power series ...");
    put("Give the tolerance on the residual : "); get(tolres);
    step := Series_and_Predictors.Set_Step_Size
              (standard_output,eva,tolcff,tolres,true);
    put("The computed step size : "); put(step,3); new_line;
    Test_Series_Predictors.Standard_Check_Prediction(hom,srv,eva,step);
    Test_Pade_Predictors.Standard_Check_Prediction(hom,srv,eva,pv,step);
  end Standard_Test_Prediction;

  procedure DoblDobl_Test_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tests the predictor on the solution in sol for the homotopy hom,
  --   in double double precision.

    nit : integer32 := 4;
    srv : DoblDobl_Dense_Series_Vectors.Vector(1..sol.n);
    eva : DoblDobl_Dense_Series_Vectors.Vector(hom'range);
    tolcff : constant double_float := 1.0e-12;
    tolres,step : double_float := 0.0;
    dd_step : double_double;
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range);

  begin
    new_line;
    put("Give the number of Newton iterations : "); get(nit);
    Test_Pade_Predictors.DoblDobl_Test_Pade_Prediction(hom,sol,nit,pv);
    Series_and_Predictors.Newton_Prediction(nit,hom,sol.v,srv,eva);
    Test_Pade_Predictors.DoblDobl_Step_Prediction(hom,srv,eva,pv);
    new_line;
    put_line("Setting the step size based on the power series ...");
    put("Give the tolerance on the residual : "); get(tolres);
    step := Series_and_Predictors.Set_Step_Size
              (standard_output,eva,tolcff,tolres,true);
    put("The computed step size : "); put(step,3); new_line;
    dd_step := create(step);
    Test_Series_Predictors.DoblDobl_Check_Prediction(hom,srv,eva,dd_step);
    Test_Pade_Predictors.DoblDobl_Check_Prediction(hom,srv,eva,pv,dd_step);
  end DoblDobl_Test_Prediction;

  procedure QuadDobl_Test_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Tests the predictor on the solution in sol for the homotopy hom,
  --   in quad double precision.

    nit : integer32 := 4;
    srv : QuadDobl_Dense_Series_Vectors.Vector(1..sol.n);
    eva : QuadDobl_Dense_Series_Vectors.Vector(hom'range);
    tolcff : constant double_float := 1.0e-12;
    tolres,step : double_float := 0.0;
    qd_step : quad_double;
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range);

  begin
    new_line;
    put("Give the number of Newton iterations : "); get(nit);
    Test_Pade_Predictors.QuadDobl_Test_Pade_Prediction(hom,sol,nit,pv);
    Series_and_Predictors.Newton_Prediction(nit,hom,sol.v,srv,eva);
    Test_Pade_Predictors.QuadDobl_Step_Prediction(hom,srv,eva,pv);
    new_line;
    put_line("Setting the step size based on the power series ...");
    put("Give the tolerance on the residual : "); get(tolres);
    step := Series_and_Predictors.Set_Step_Size
              (standard_output,eva,tolcff,tolres,true);
    put("The computed step size : "); put(step,3); new_line;
    qd_step := create(step);
    Test_Series_Predictors.QuadDobl_Check_Prediction(hom,srv,eva,qd_step);
    Test_Pade_Predictors.QuadDobl_Check_Prediction(hom,srv,eva,pv,qd_step);
  end QuadDobl_Test_Prediction;

  procedure Standard_Test_Prediction
              ( nq,idxpar : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The Standard_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_Series_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,idxpar);
    len : constant integer32
        := integer32(Standard_Complex_Solutions.Length_Of(sols));
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    sol : Standard_Complex_Solutions.Link_to_Solution;

  begin
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

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_Series_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,idxpar);
    len : constant integer32
        := integer32(DoblDobl_Complex_Solutions.Length_Of(sols));
    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    sol : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
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

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_Series_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,idxpar);
    len : constant integer32
        := integer32(QuadDobl_Complex_Solutions.Length_Of(sols));
    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    sol : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
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
        dropsols : Standard_Complex_Solutions.Solution_List
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
        dropsols : DoblDobl_Complex_Solutions.Solution_List
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
        dropsols : QuadDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        QuadDobl_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision
  --   and then launches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
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
