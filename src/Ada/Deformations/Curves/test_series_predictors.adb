with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
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
with Homotopy_Series_Readers;
with Series_and_Homotopies;
with Series_and_Predictors;

package body Test_Series_Predictors is

  procedure Standard_Check_Prediction
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Complex_Series_Vectors.Vector;
                step : in double_float ) is

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
    put("The actual error of the series predictor : "); put(nrm,3); new_line;
    Standard_Complex_Poly_Systems.Clear(phm);
  end Standard_Check_Prediction;

  procedure DoblDobl_Check_Prediction
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Complex_Series_Vectors.Vector;
                step : in double_double ) is

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
    put("The actual error of the series predictor : "); put(nrm,3); new_line;
    DoblDobl_Complex_Poly_Systems.Clear(phm);
  end DoblDobl_Check_Prediction;

  procedure QuadDobl_Check_Prediction
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Complex_Series_Vectors.Vector;
                step : in quad_double ) is

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
    put("The actual error of the series predictor : "); put(nrm,3); new_line;
    QuadDobl_Complex_Poly_Systems.Clear(phm);
  end QuadDobl_Check_Prediction;

  procedure Standard_Homotopy_Reader
              ( nbeq,idxpar : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    ans : character;
    nvr : integer32;

  begin
    new_line;
    put("Natural parameter homotopy ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Homotopy_Series_Readers.Standard_Parameter_Reader(nbeq,nvr,idxpar,sols);
    else
      idxpar := 0;
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
    end if;
  end Standard_Homotopy_Reader;

  procedure DoblDobl_Homotopy_Reader
              ( nbeq,idxpar : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    ans : character;
    nvr : integer32;

  begin
    new_line;
    put("Natural parameter homotopy ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Homotopy_Series_Readers.DoblDobl_Parameter_Reader(nbeq,nvr,idxpar,sols);
    else
      idxpar := 0;
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
    end if;
  end DoblDobl_Homotopy_Reader;

  procedure QuadDobl_Homotopy_Reader
              ( nbeq,idxpar : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    ans : character;
    nvr : integer32;

  begin
    new_line;
    put("Natural parameter homotopy ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Homotopy_Series_Readers.QuadDobl_Parameter_Reader(nbeq,nvr,idxpar,sols);
    else
      idxpar := 0;
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
    end if;
  end QuadDobl_Homotopy_Reader;

end Test_Series_Predictors;
