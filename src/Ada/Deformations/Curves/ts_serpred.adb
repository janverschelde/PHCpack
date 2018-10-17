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
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Dense_Series;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_Vectors_io;
with Standard_Series_Poly_Systems;
with Standard_Series_Poly_SysFun;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Dense_Series_Vectors_io;
with DoblDobl_Series_Poly_Systems;
with DoblDobl_Series_Poly_SysFun;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_Vectors_io;
with QuadDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_SysFun;
with Series_and_Polynomials;
with Series_and_Polynomials_io;          use Series_and_Polynomials_io;
with Power_Series_Methods;               use Power_Series_Methods;
with Series_and_Solutions;
with Series_and_Homotopies;
with Series_and_Predictors;
with Homotopy_Series_Readers;
with Standard_Pade_Approximants;
with Standard_Pade_Approximants_io;
with DoblDobl_Pade_Approximants;
with DoblDobl_Pade_Approximants_io;
with QuadDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants_io;
with Homotopy_Pade_Approximants;

procedure ts_serpred is

-- DESCRIPTION :
--   Test on the application of power series to predict solutions.

  procedure Standard_Check_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Dense_Series_Vectors.Vector;
                step : in double_float ) is

  -- DESCRIPTION :
  --   Checks the quality of the predicted step.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv;
  --   step     a step size.

    pred_err : double_float;
    pred_sol : Standard_Complex_Vectors.Vector(srv'range);
    phm : Standard_Complex_Poly_Systems.Poly_Sys(hom'range);
    val : Standard_Complex_Vectors.Vector(phm'range);
    nrm : double_float;

  begin
    pred_err := Series_and_Predictors.Predicted_Error(eva,step);
    put("The predicted error : "); put(pred_err,3); new_line;
    pred_sol := Series_and_Predictors.Predicted_Solution(srv,step);
    put_line("The predicted solution :"); put_line(pred_sol);
    phm := Series_and_Homotopies.Eval(hom,step);
    val := Standard_Complex_Poly_SysFun.Eval(phm,pred_sol);
    nrm := Standard_Complex_Vector_Norms.Max_Norm(val);
    put("The actual error : "); put(nrm,3); new_line;
    Standard_Complex_Poly_Systems.Clear(phm);
  end Standard_Check_Prediction;

  procedure DoblDobl_Check_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Dense_Series_Vectors.Vector;
                step : in double_double ) is

  -- DESCRIPTION :
  --   Checks the quality of the predicted step.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv;
  --   step     a step size.

    pred_err : double_double;
    pred_sol : DoblDobl_Complex_Vectors.Vector(srv'range);
    phm : DoblDobl_Complex_Poly_Systems.Poly_Sys(hom'range);
    val : DoblDobl_Complex_Vectors.Vector(phm'range);
    nrm : double_double;

  begin
    pred_err := Series_and_Predictors.Predicted_Error(eva,step);
    put("The predicted error : "); put(pred_err,3); new_line;
    pred_sol := Series_and_Predictors.Predicted_Solution(srv,step);
    put_line("The predicted solution :"); put_line(pred_sol);
    phm := Series_and_Homotopies.Eval(hom,step);
    val := DoblDobl_Complex_Poly_SysFun.Eval(phm,pred_sol);
    nrm := DoblDobl_Complex_Vector_Norms.Max_Norm(val);
    put("The actual error : "); put(nrm,3); new_line;
    DoblDobl_Complex_Poly_Systems.Clear(phm);
  end DoblDobl_Check_Prediction;

  procedure QuadDobl_Check_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Dense_Series_Vectors.Vector;
                step : in quad_double ) is

  -- DESCRIPTION :
  --   Checks the quality of the predicted step.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv;
  --   step     a step size.

    pred_err : quad_double;
    pred_sol : QuadDobl_Complex_Vectors.Vector(srv'range);
    phm : QuadDobl_Complex_Poly_Systems.Poly_Sys(hom'range);
    val : QuadDobl_Complex_Vectors.Vector(phm'range);
    nrm : quad_double;

  begin
    pred_err := Series_and_Predictors.Predicted_Error(eva,step);
    put("The predicted error : "); put(pred_err,3); new_line;
    pred_sol := Series_and_Predictors.Predicted_Solution(srv,step);
    put_line("The predicted solution :"); put_line(pred_sol);
    phm := Series_and_Homotopies.Eval(hom,step);
    val := QuadDobl_Complex_Poly_SysFun.Eval(phm,pred_sol);
    nrm := QuadDobl_Complex_Vector_Norms.Max_Norm(val);
    put("The actual error : "); put(nrm,3); new_line;
    QuadDobl_Complex_Poly_Systems.Clear(phm);
  end QuadDobl_Check_Prediction;

  procedure Standard_Step_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Tests the determination of the step size interactively
  --   to gain computational experience with its effectiveness.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv.

    step : double_float := 0.0;
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    loop
      put("Give the step size : "); get(step);
      Standard_Check_Prediction(hom,srv,eva,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Step_Prediction;

  procedure DoblDobl_Step_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Tests the determination of the step size interactively
  --   to gain computational experience with its effectiveness.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv.

    step : double_double := create(0.0);
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    loop
      put("Give the step size : "); get(step);
      DoblDobl_Check_Prediction(hom,srv,eva,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Step_Prediction;

  procedure QuadDobl_Step_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Tests the determination of the step size interactively
  --   to gain computational experience with its effectiveness.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv.
 
    step : quad_double := create(0.0);
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    loop
      put("Give the step size : "); get(step);
      QuadDobl_Check_Prediction(hom,srv,eva,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end QuadDobl_Step_Prediction;

  procedure Standard_Test_Pade_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Solutions.Solution;
                nit : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the Pade predictor on the solution in sol for the 
  --   homotopy hom, in standard double precision.

    degnum,degden : integer32 := 0;
    nbt : natural32;
    nbeq : constant integer32 := hom'last;
    srv : Standard_Dense_Series_Vectors.Vector(sol.v'range);
    eva : Standard_Dense_Series_Vectors.Vector(1..nbeq);
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range);
    poles : Standard_Complex_VecVecs.VecVec(pv'range);

  begin
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    nbt := natural32(degnum+degden+1);
    Homotopy_Pade_Approximants.Standard_Pade_Approximant
      (sol.v,nbeq+1,nbeq,degnum,degden,natural32(nit),srv,eva,pv,true);
    put_line("The solution series :");
    Standard_Dense_Series_Vectors_io.put(srv);
    put_line("The evaluated solution series :");
    Standard_Dense_Series_Vectors_io.put(eva);
    put_line("The Pade approximant :");
    for i in pv'range loop
      put("numerator for component "); put(i,1); put_line(" :");
      put_line(Standard_Pade_Approximants.Numerator_Coefficients(pv(i)));
      put("denominator for component "); put(i,1); put_line(" :");
      put_line(Standard_Pade_Approximants.Denominator_Coefficients(pv(i)));
      put_line(Standard_Pade_Approximants_io.Write(pv(i)));
    end loop;
    poles := Homotopy_Pade_Approximants.Standard_Poles(pv);
    put_line("The poles : "); put_line(poles);
    Standard_Pade_Approximants.Clear(pv);
  end Standard_Test_Pade_Prediction;

  procedure DoblDobl_Test_Pade_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Solutions.Solution;
                nit : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the Pade predictor on the solution in sol for the 
  --   homotopy hom, in double double precision.

    degnum,degden : integer32 := 0;
    nbt : natural32;
    nbeq : constant integer32 := hom'last;
    srv : DoblDobl_Dense_Series_Vectors.Vector(sol.v'range);
    eva : DoblDobl_Dense_Series_Vectors.Vector(1..nbeq);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range);

  begin
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    nbt := natural32(degnum+degden+1);
    Homotopy_Pade_Approximants.DoblDobl_Pade_Approximant
      (sol.v,nbeq+1,nbeq,degnum,degden,natural32(nit),srv,eva,pv);
    put_line("The solution series :");
    DoblDobl_Dense_Series_Vectors_io.put(srv);
    put_line("The evaluated solution series :");
    DoblDobl_Dense_Series_Vectors_io.put(eva);
    put_line("The Pade approximant :");
    for i in pv'range loop
      put_line(DoblDobl_Pade_Approximants_io.Write(pv(i)));
    end loop;
    poles := Homotopy_Pade_Approximants.DoblDobl_Poles(pv);
    put_line("The poles : "); put_line(poles);
    DoblDobl_Pade_Approximants.Clear(pv);
  end DoblDobl_Test_Pade_Prediction;

  procedure QuadDobl_Test_Pade_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Solutions.Solution;
                nit : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the Pade predictor on the solution in sol for the 
  --   homotopy hom, in quad double precision.

    degnum,degden : integer32 := 0;
    nbt : natural32;
    nbeq : constant integer32 := hom'last;
    srv : QuadDobl_Dense_Series_Vectors.Vector(sol.v'range);
    eva : QuadDobl_Dense_Series_Vectors.Vector(1..nbeq);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range);

  begin
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    nbt := natural32(degnum+degden+1);
    Homotopy_Pade_Approximants.QuadDobl_Pade_Approximant
      (sol.v,nbeq+1,nbeq,degnum,degden,natural32(nit),srv,eva,pv);
    put_line("The solution series :");
    QuadDobl_Dense_Series_Vectors_io.put(srv);
    put_line("The evaluated solution series :");
    QuadDobl_Dense_Series_Vectors_io.put(eva);
    put_line("The Pade approximant :");
    for i in pv'range loop
      put_line(QuadDobl_Pade_Approximants_io.Write(pv(i)));
    end loop;
    poles := Homotopy_Pade_Approximants.QuadDobl_Poles(pv);
    put_line("The poles : "); put_line(poles);
    QuadDobl_Pade_Approximants.Clear(pv);
  end QuadDobl_Test_Pade_Prediction;

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
    ans : character;

  begin
    new_line;
    put("Give the number of Newton iterations : "); get(nit);
    new_line;
    put("Compute a Pade approximation ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Standard_Test_Pade_Prediction(hom,sol,nit);
    else
      Series_and_Predictors.Newton_Prediction(nit,hom,sol.v,srv,eva);
      Standard_Step_Prediction(hom,srv,eva);
      new_line;
      put_line("Setting the step size based on the power series ...");
      put("Give the tolerance on the residual : "); get(tolres);
      step := Series_and_Predictors.Set_Step_Size
                (standard_output,eva,tolcff,tolres,true);
      put("The computed step size : "); put(step,3); new_line;
      Standard_Check_Prediction(hom,srv,eva,step);
    end if;
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
    ans : character;

  begin
    new_line;
    put("Give the number of Newton iterations : "); get(nit);
    skip_line; -- needed for prompting of double double ...
    new_line;
    put("Compute a Pade approximation ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      DoblDobl_Test_Pade_Prediction(hom,sol,nit);
    else
      Series_and_Predictors.Newton_Prediction(nit,hom,sol.v,srv,eva);
      DoblDobl_Step_Prediction(hom,srv,eva);
      new_line;
      put_line("Setting the step size based on the power series ...");
      put("Give the tolerance on the residual : "); get(tolres);
      step := Series_and_Predictors.Set_Step_Size
                (standard_output,eva,tolcff,tolres,true);
      put("The computed step size : "); put(step,3); new_line;
      dd_step := create(step);
      DoblDobl_Check_Prediction(hom,srv,eva,dd_step);
    end if;
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
    ans : character;

  begin
    new_line;
    put("Give the number of Newton iterations : "); get(nit);
    skip_line; -- needed for prompting of quad double ...
    new_line;
    put("Compute a Pade approximation ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      QuadDobl_Test_Pade_Prediction(hom,sol,nit);
    else
      Series_and_Predictors.Newton_Prediction(nit,hom,sol.v,srv,eva);
      QuadDobl_Step_Prediction(hom,srv,eva);
      new_line;
      put_line("Setting the step size based on the power series ...");
      put("Give the tolerance on the residual : "); get(tolres);
      step := Series_and_Predictors.Set_Step_Size
                (standard_output,eva,tolcff,tolres,true);
      put("The computed step size : "); put(step,3); new_line;
      qd_step := create(step);
      QuadDobl_Check_Prediction(hom,srv,eva,qd_step);
    end if;
  end QuadDobl_Test_Prediction;

  procedure Standard_Test_Prediction
              ( nq : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The Standard_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_Series_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1);
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
              ( nq : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The DoblDobl_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_Series_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1);
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
              ( nq : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The QuadDobl_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_Series_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1);
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

    nbeq : integer32;
    sols : Standard_Complex_Solutions.Solution_List;
    ans : character;

  begin
    new_line;
    put("Random gamma ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Homotopy_Series_Readers.Standard_Reader(nbeq,sols,tpow=>1);
    else
      declare
        gamma : constant Standard_Complex_Numbers.Complex_Number
              := Standard_Complex_Numbers.Create(1.0);
      begin
        Homotopy_Series_Readers.Standard_Reader(nbeq,sols,1,gamma);
      end;
    end if;
    new_line;
    Standard_Test_Prediction(nbeq,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in double double precision.

    nbeq : integer32;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    ans : character;

  begin
    new_line;
    put("Random gamma ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols,tpow=>1);
    else
      declare
        one : constant double_double := create(1.0);
        gamma : constant DoblDobl_Complex_Numbers.Complex_Number
              := DoblDobl_Complex_Numbers.Create(one);
      begin
        Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols,1,gamma);
      end;
    end if;
    new_line;
    DoblDobl_Test_Prediction(nbeq,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in quad double precision.

    nbeq : integer32;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    ans : character;

  begin
    new_line;
    put("Random gamma ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols,tpow=>1);
    else
      declare
        one : constant quad_double := create(1.0);
        gamma : constant QuadDobl_Complex_Numbers.Complex_Number
              := QuadDobl_Complex_Numbers.Create(one);
      begin
        Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols,1,gamma);
      end;
    end if;
    new_line;
    QuadDobl_Test_Prediction(nbeq,sols);
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
end ts_serpred;
