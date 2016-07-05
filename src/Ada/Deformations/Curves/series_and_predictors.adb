with Standard_Complex_Vector_Norms;
with Standard_Series_Vector_Functions;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Series_Vector_Functions;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Series_Vector_Functions;

package body Series_and_Predictors is

  function Predicted_Error
             ( evls : Standard_Dense_Series_Vectors.Vector;
               step : double_float ) return double_float is

    eva : constant Standard_Complex_Vectors.Vector
        := Standard_Series_Vector_Functions.Eval(evls,step);
    res : constant double_float
        := Standard_Complex_Vector_Norms.Max_Norm(eva);

  begin
    return res;
  end Predicted_Error;

  function Predicted_Error
             ( evls : DoblDobl_Dense_Series_Vectors.Vector;
               step : double_double ) return double_double is

    eva : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_Series_Vector_Functions.Eval(evls,step);
    res : constant double_double
        := DoblDobl_Complex_Vector_Norms.Max_Norm(eva);

  begin
    return res;
  end Predicted_Error;

  function Predicted_Error
             ( evls : QuadDobl_Dense_Series_Vectors.Vector;
               step : quad_double ) return quad_double is

    eva : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_Series_Vector_Functions.Eval(evls,step);
    res : constant quad_double
        := QuadDobl_Complex_Vector_Norms.Max_Norm(eva);

  begin
    return res;
  end Predicted_Error;

  function Predicted_Solution
             ( srv : Standard_Dense_Series_Vectors.Vector;
               step : double_float )
             return Standard_Complex_Vectors.Vector is

    res : constant Standard_Complex_Vectors.Vector
        := Standard_Series_Vector_Functions.Eval(srv,step);

  begin
    return res;
  end Predicted_Solution;

  function Predicted_Solution
             ( srv : DoblDobl_Dense_Series_Vectors.Vector;
               step : double_double )
             return DoblDobl_Complex_Vectors.Vector is

    res : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_Series_Vector_Functions.Eval(srv,step);

  begin
    return res;
  end Predicted_Solution;

  function Predicted_Solution
             ( srv : QuadDobl_Dense_Series_Vectors.Vector;
               step : quad_double )
             return QuadDobl_Complex_Vectors.Vector is

    res : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_Series_Vector_Functions.Eval(srv,step);

  begin
    return res;
  end Predicted_Solution;

end Series_and_Predictors;
