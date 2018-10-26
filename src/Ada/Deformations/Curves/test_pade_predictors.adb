with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with Standard_Dense_Series_Vectors_io;
with DoblDobl_Dense_Series_Vectors_io;
with QuadDobl_Dense_Series_Vectors_io;
with QuadDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_SysFun;
with Series_and_Homotopies;
with Series_and_Predictors;
with Test_Series_Predictors;
with Standard_Pade_Approximants_io;
with DoblDobl_Pade_Approximants_io;
with QuadDobl_Pade_Approximants_io;
with Homotopy_Pade_Approximants;

package body Test_Pade_Predictors is

  procedure Standard_Check_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Dense_Series_Vectors.Vector;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                step : in double_float ) is

    phm : Standard_Complex_Poly_Systems.Poly_Sys(hom'range);
    valxpv : Standard_Complex_Vectors.Vector(phm'range);
    nrmxpv : double_float;
    xpv : constant Standard_Complex_Vectors.Vector(pv'range)
        := Standard_Pade_Approximants.Eval(pv,step);

  begin
    put_line("The evaluated Pade approximation :"); put_line(xpv);
    phm := Series_and_Homotopies.Eval(hom,step);
    valxpv := Standard_Complex_Poly_SysFun.Eval(phm,xpv);
    nrmxpv := Standard_Complex_Vector_Norms.Max_Norm(valxpv);
    put("The actual error of the Pade predictor : "); put(nrmxpv,3); new_line;
    Standard_Complex_Poly_Systems.Clear(phm);
  end Standard_Check_Prediction;

  procedure DoblDobl_Check_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Dense_Series_Vectors.Vector;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector;
                step : in double_double ) is

    phm : DoblDobl_Complex_Poly_Systems.Poly_Sys(hom'range);
    valxpv : DoblDobl_Complex_Vectors.Vector(phm'range);
    nrmxpv : double_double;
    xpv : constant DoblDobl_Complex_Vectors.Vector(pv'range)
        := DoblDobl_Pade_Approximants.Eval(pv,step);

  begin
    put_line("The evaluated Pade approximant :"); put_line(xpv);
    phm := Series_and_Homotopies.Eval(hom,step);
    valxpv := DoblDobl_Complex_Poly_SysFun.Eval(phm,xpv);
    nrmxpv := DoblDobl_Complex_Vector_Norms.Max_Norm(valxpv);
    put("The actual error of the Pade predictor : "); put(nrmxpv,3); new_line;
    DoblDobl_Complex_Poly_Systems.Clear(phm);
  end DoblDobl_Check_Prediction;

  procedure QuadDobl_Check_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Dense_Series_Vectors.Vector;
                pv : in QuadDobl_Pade_Approximants.Pade_Vector;
                step : in quad_double ) is

    phm : QuadDobl_Complex_Poly_Systems.Poly_Sys(hom'range);
    valxpv : QuadDobl_Complex_Vectors.Vector(phm'range);
    nrmxpv : quad_double;
    xpv : constant QuadDobl_Complex_Vectors.Vector(pv'range)
        := QuadDobl_Pade_Approximants.Eval(pv,step);

  begin
    put_line("The evaluated Pade approximant :"); put_line(xpv);
    phm := Series_and_Homotopies.Eval(hom,step);
    valxpv := QuadDobl_Complex_Poly_SysFun.Eval(phm,xpv);
    nrmxpv := QuadDobl_Complex_Vector_Norms.Max_Norm(valxpv);
    put("The actual error of the Pade predictor : "); put(nrmxpv,3); new_line;
    QuadDobl_Complex_Poly_Systems.Clear(phm);
  end QuadDobl_Check_Prediction;

  procedure Standard_Step_Prediction
              ( hom : in Standard_Series_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Dense_Series_Vectors.Vector;
                pv : in Standard_Pade_Approximants.Pade_Vector ) is

    step : double_float := 0.0;
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    loop
      put("Give the step size : "); get(step);
      Test_Series_Predictors.Standard_Check_Prediction(hom,srv,eva,step);
      Standard_Check_Prediction(hom,srv,eva,pv,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Step_Prediction;

  procedure DoblDobl_Step_Prediction
              ( hom : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Dense_Series_Vectors.Vector;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector ) is

    step : double_double := create(0.0);
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    skip_line; -- needed for the prompting of the double double step
    loop
      put("Give the step size : "); get(step);
      Test_Series_Predictors.DoblDobl_Check_Prediction(hom,srv,eva,step);
      DoblDobl_Check_Prediction(hom,srv,eva,pv,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Step_Prediction;

  procedure QuadDobl_Step_Prediction
              ( hom : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Dense_Series_Vectors.Vector;
                pv : QuadDobl_Pade_Approximants.Pade_Vector ) is

    step : quad_double := create(0.0);
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    skip_line; -- needed for prompting of the quad double step
    loop
      put("Give the step size : "); get(step);
      Test_Series_Predictors.QuadDobl_Check_Prediction(hom,srv,eva,step);
      QuadDobl_Check_Prediction(hom,srv,eva,pv,step);
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

end Test_Pade_Predictors;
