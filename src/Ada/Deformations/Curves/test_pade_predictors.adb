with text_io;                            use text_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_SysFun;
with Series_and_Homotopies;
with Series_and_Predictors;

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

end Test_Pade_Predictors;
