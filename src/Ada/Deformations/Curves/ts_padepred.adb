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
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_Vectors_io;
with Standard_Series_Poly_Systems;
with Standard_Series_Poly_SysFun;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Dense_Series_Vectors_io;
with DoblDobl_Series_Poly_Systems;
with DoblDobl_Series_Poly_SysFun;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_Vectors_io;
with QuadDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_SysFun;
with Series_and_Homotopies;
with Series_and_Predictors;
with Standard_Pade_Approximants;
with Standard_Pade_Approximants_io;
with DoblDobl_Pade_Approximants;
with DoblDobl_Pade_Approximants_io;
with QuadDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants_io;
with Homotopy_Pade_Approximants;
with Test_Series_Predictors;
with Test_Pade_Predictors;

procedure ts_padepred is

-- DESCRIPTION :
--   Test on the application of Pade approximants to predict solutions.

  procedure Standard_Step_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Dense_Series_Vectors.Vector;
                pv : in Standard_Pade_Approximants.Pade_Vector ) is

  -- DESCRIPTION :
  --   Tests the determination of the step size interactively
  --   to gain computational experience with its effectiveness.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv;
  --   pv       a vector of Pade approximants.

    step : double_float := 0.0;
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    loop
      put("Give the step size : "); get(step);
      Test_Series_Predictors.Standard_Check_Prediction(hom,srv,eva,step);
      Test_Pade_Predictors.Standard_Check_Prediction(hom,srv,eva,pv,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Step_Prediction;

  procedure DoblDobl_Step_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Dense_Series_Vectors.Vector;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector ) is

  -- DESCRIPTION :
  --   Tests the determination of the step size interactively
  --   to gain computational experience with its effectiveness.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv;
  --   pv       a vector of Pade approximants.

    step : double_double := create(0.0);
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    loop
      put("Give the step size : "); get(step);
      Test_Series_Predictors.DoblDobl_Check_Prediction(hom,srv,eva,step);
      Test_Pade_Predictors.DoblDobl_Check_Prediction(hom,srv,eva,pv,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Step_Prediction;

  procedure QuadDobl_Step_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Dense_Series_Vectors.Vector;
                pv : QuadDobl_Pade_Approximants.Pade_Vector ) is

  -- DESCRIPTION :
  --   Tests the determination of the step size interactively
  --   to gain computational experience with its effectiveness.

  -- ON ENTRY :
  --   hom      homotopy with series coefficients, where the parameter
  --            in the series is the continuation parameter;
  --   srv      series development for the solution;
  --   eva      evaluation of hom at srv;
  --   pv       a vector of Pade approximants.
 
    step : quad_double := create(0.0);
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    loop
      put("Give the step size : "); get(step);
      Test_Series_Predictors.QuadDobl_Check_Prediction(hom,srv,eva,step);
      Test_Pade_Predictors.QuadDobl_Check_Prediction(hom,srv,eva,pv,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end QuadDobl_Step_Prediction;

  procedure Standard_Test_Pade_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Solutions.Solution;
                nit : in integer32;
                pv : out Standard_Pade_Approximants.Pade_Vector ) is

  -- DESCRIPTION :
  --   Tests the Pade predictor on the solution in sol for the 
  --   homotopy hom, in standard double precision.

    degnum,degden : integer32 := 0;
    nbt : natural32;
    nbeq : constant integer32 := hom'last;
    srv : Standard_Dense_Series_Vectors.Vector(sol.v'range);
    eva : Standard_Dense_Series_Vectors.Vector(1..nbeq);
    poles : Standard_Complex_VecVecs.VecVec(pv'range);
    lead,idx : integer32;
    minpole : double_float;

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
    Homotopy_Pade_Approximants.Smallest_Forward_Pole(poles,lead,idx,minpole);
    put("The radius of the smallest pole : "); put(minpole,3); new_line;
  end Standard_Test_Pade_Prediction;

  procedure DoblDobl_Test_Pade_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Solutions.Solution;
                nit : in integer32;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector ) is

  -- DESCRIPTION :
  --   Tests the Pade predictor on the solution in sol for the 
  --   homotopy hom, in double double precision.

    degnum,degden : integer32 := 0;
    nbt : natural32;
    nbeq : constant integer32 := hom'last;
    srv : DoblDobl_Dense_Series_Vectors.Vector(sol.v'range);
    eva : DoblDobl_Dense_Series_Vectors.Vector(1..nbeq);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range);
    lead,idx : integer32;
    minpole : double_double;

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
    Homotopy_Pade_Approximants.Smallest_Forward_Pole(poles,lead,idx,minpole);
    put("The radius of the smallest pole : "); put(minpole,3); new_line;
  end DoblDobl_Test_Pade_Prediction;

  procedure QuadDobl_Test_Pade_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Solutions.Solution;
                nit : in integer32;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector ) is

  -- DESCRIPTION :
  --   Tests the Pade predictor on the solution in sol for the 
  --   homotopy hom, in quad double precision.

    degnum,degden : integer32 := 0;
    nbt : natural32;
    nbeq : constant integer32 := hom'last;
    srv : QuadDobl_Dense_Series_Vectors.Vector(sol.v'range);
    eva : QuadDobl_Dense_Series_Vectors.Vector(1..nbeq);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range);
    lead,idx : integer32;
    minpole : quad_double;

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
    Homotopy_Pade_Approximants.Smallest_Forward_Pole(poles,lead,idx,minpole);
    put("The radius of the smallest pole : "); put(minpole,3); new_line;
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
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range);

  begin
    new_line;
    put("Give the number of Newton iterations : "); get(nit);
    Standard_Test_Pade_Prediction(hom,sol,nit,pv);
    Series_and_Predictors.Newton_Prediction(nit,hom,sol.v,srv,eva);
    Standard_Step_Prediction(hom,srv,eva,pv);
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
    skip_line; -- needed for prompting of double double ...
    DoblDobl_Test_Pade_Prediction(hom,sol,nit,pv);
    Series_and_Predictors.Newton_Prediction(nit,hom,sol.v,srv,eva);
    DoblDobl_Step_Prediction(hom,srv,eva,pv);
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
    skip_line; -- needed for prompting of quad double ...
    QuadDobl_Test_Pade_Prediction(hom,sol,nit,pv);
    Series_and_Predictors.Newton_Prediction(nit,hom,sol.v,srv,eva);
    QuadDobl_Step_Prediction(hom,srv,eva,pv);
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

  begin
    Test_Series_Predictors.Standard_Homotopy_Reader(nbeq,sols);
    new_line;
    Standard_Test_Prediction(nbeq,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in double double precision.

    nbeq : integer32;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.DoblDobl_Homotopy_Reader(nbeq,sols);
    new_line;
    DoblDobl_Test_Prediction(nbeq,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in quad double precision.

    nbeq : integer32;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.QuadDobl_Homotopy_Reader(nbeq,sols);
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
end ts_padepred;
