with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with Standard_Dense_Series;
with Standard_Dense_Series_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_Vectors;

package Series_and_Predictors is

-- DESCRIPTION :
--   Evaluates power series to predict a solution.

  function Predicted_Error
             ( evls : Standard_Dense_Series_Vectors.Vector;
               step : double_float ) return double_float;
  function Predicted_Error
             ( evls : DoblDobl_Dense_Series_Vectors.Vector;
               step : double_double ) return double_double;
  function Predicted_Error
             ( evls : QuadDobl_Dense_Series_Vectors.Vector;
               step : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Given in evls the vector of power series evaluated at
  --   a power series solution and in step the step size,
  --   return the max norm of the series evls evaluated at step,
  --   in double, double double, or quad double precision.

  function Least_Order
             ( s : Standard_Dense_Series.Series; tol : double_float )
             return integer32;
  function Least_Order
             ( s : DoblDobl_Dense_Series.Series; tol : double_float )
             return integer32;
  function Least_Order
             ( s : QuadDobl_Dense_Series.Series; tol : double_float )
             return integer32;

  -- DESCRIPTION :
  --   Returns the smallest integer k in the range 0..s.order
  --   for which AbsVal(s.cff(k)) > tol.
  --   If all coefficients are less than tol, then s.order+1 is returned.

  procedure Least_Order
             ( v : in Standard_Dense_Series_Vectors.Vector;
               tol : in double_float; vk,ord : out integer32 );
  procedure Least_Order
             ( v : in DoblDobl_Dense_Series_Vectors.Vector;
               tol : in double_float; vk,ord : out integer32 );
  procedure Least_Order
             ( v : in QuadDobl_Dense_Series_Vectors.Vector;
               tol : in double_float; vk,ord : out integer32 );

  -- DESCRIPTION :
  --   For all series s in v, returns in ord the smallest Least_Order(s,tol)
  --   and in vk the compoent of vk for which ord = Least_Order(v(vk),tol).
  --   The corresponding coefficient with t^ord is then in v(vk).cff(ord).

  function Set_Step_Size
             ( v : Standard_Dense_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float;
  function Set_Step_Size
             ( v : DoblDobl_Dense_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float;

  -- DESCRIPTION :
  --   Computes a step size so that the residual will be of order tolres.
  --   The tolcff is needed to compute the least order of v.

  function Predicted_Solution
             ( srv : Standard_Dense_Series_Vectors.Vector;
               step : double_float )
             return Standard_Complex_Vectors.Vector;
  function Predicted_Solution
             ( srv : DoblDobl_Dense_Series_Vectors.Vector;
               step : double_double )
             return DoblDobl_Complex_Vectors.Vector;
  function Predicted_Solution
             ( srv : QuadDobl_Dense_Series_Vectors.Vector;
               step : quad_double )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the vector of power series srv at the step,
  --   in double, double double, or quad double precision.

end Series_and_Predictors;
