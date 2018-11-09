with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Series;
with Standard_Complex_Series_Vectors;
with Standard_CSeries_Poly_Systems;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_CSeries_Poly_Systems;

package Series_and_Predictors is

-- DESCRIPTION :
--   Evaluates power series to predict a solution.

  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Vectors.Vector;
                srv : out Standard_Complex_Series_Vectors.Vector;
                eva : out Standard_Complex_Series_Vectors.Vector );
  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector;
                srv : out DoblDobl_Complex_Series_Vectors.Vector;
                eva : out DoblDobl_Complex_Series_Vectors.Vector );
  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Vectors.Vector;
                srv : out QuadDobl_Complex_Series_Vectors.Vector;
                eva : out QuadDobl_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   Applies Newton's method on the homotopy in hom,
  --   starting at the coordinates of the solution as the leading
  --   coefficients in the power series,
  --   in double, double double, or quad double precision.
  --   The versions are silent, do not write extra output.

  -- ON ENTRY :
  --   maxdeg   maximal degree of the series;
  --   nit      number of iterations with Newton's method;
  --   hom      a homotopy with coefficients as power series,
  --            where the series parameter is the continuation parameter;
  --   sol      solution of the homotopy for t = 0, its coordinates
  --            contain the leading coefficients of the power series.

  -- ON RETURN :
  --   srv      vector of power series, to predict the solution
  --            for some small positive value of the continuation parameter,
  --            srv'range = sol'range;
  --   eva      evaluated series vector srv in the homotopy hom;
  --            eva'range = hom'range.

  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Vectors.Vector;
                srv : out Standard_Complex_Series_Vectors.Vector;
                eva : out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector;
                srv : out DoblDobl_Complex_Series_Vectors.Vector;
                eva : out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Vectors.Vector;
                srv : out QuadDobl_Complex_Series_Vectors.Vector;
                eva : out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Applies Newton's method on the homotopy in hom,
  --   starting at the coordinates of the solution as the leading
  --   coefficients in the power series,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write extra diagnostic output to;
  --   maxdeg   maximal degree of the series;
  --   nit      number of iterations with Newton's method;
  --   hom      a homotopy with coefficients as power series,
  --            where the series parameter is the continuation parameter;
  --   sol      solution of the homotopy for t = 0, its coordinates
  --            contain the leading coefficients of the power series;
  --   verbose  if true, then additional output is written to file.

  -- ON RETURN :
  --   srv      vector of power series, to predict the solution
  --            for some small positive value of the continuation parameter,
  --            srv'range = sol'range;
  --   eva      evaluated series vector srv in the homotopy hom;
  --            eva'range = hom'range.

  function Predicted_Error
             ( evls : Standard_Complex_Series_Vectors.Vector;
               step : double_float ) return double_float;
  function Predicted_Error
             ( evls : DoblDobl_Complex_Series_Vectors.Vector;
               step : double_double ) return double_double;
  function Predicted_Error
             ( evls : QuadDobl_Complex_Series_Vectors.Vector;
               step : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Given in evls the vector of power series evaluated at
  --   a power series solution and in step the step size,
  --   return the max norm of the series evls evaluated at step,
  --   in double, double double, or quad double precision.

  procedure Order ( v : in Standard_Complex_Series_Vectors.Vector;
                    tol : in double_float; vk,ord : out integer32 );
  procedure Order ( v : in DoblDobl_Complex_Series_Vectors.Vector;
                    tol : in double_float; vk,ord : out integer32 );
  procedure Order ( v : in QuadDobl_Complex_Series_Vectors.Vector;
                    tol : in double_float; vk,ord : out integer32 );

  -- DESCRIPTION :
  --   For all series s in v, returns in ord the smallest Order(s,tol)
  --   and in vk the compoent of vk for which ord = Order(v(vk),tol).
  --   The corresponding coefficient with t^ord is then in v(vk).cff(ord).

  function Set_Step_Size
             ( v : Standard_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float;
  function Set_Step_Size
             ( v : DoblDobl_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float;
  function Set_Step_Size
             ( v : QuadDobl_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float;

  -- DESCRIPTION :
  --   Computes a step size so that the residual will be of order tolres.
  --   The tolcff is needed to compute the least order of v.
  --   These versions are silent, they do not write extra output.

  function Set_Step_Size
             ( file : file_type;
               v : Standard_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float;
               verbose : boolean := false ) return double_float;
  function Set_Step_Size
             ( file : file_type;
               v : DoblDobl_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float;
               verbose : boolean := false ) return double_float;
  function Set_Step_Size
             ( file : file_type;
               v : QuadDobl_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float;
               verbose : boolean := false ) return double_float;

  -- DESCRIPTION :
  --   Computes a step size so that the residual will be of order tolres.
  --   The tolcff is needed to compute the least order of v.
  --   If verbose, then extra output is written to file.

  function Predicted_Solution
             ( srv : Standard_Complex_Series_Vectors.Vector;
               step : double_float )
             return Standard_Complex_Vectors.Vector;
  function Predicted_Solution
             ( srv : DoblDobl_Complex_Series_Vectors.Vector;
               step : double_double )
             return DoblDobl_Complex_Vectors.Vector;
  function Predicted_Solution
             ( srv : QuadDobl_Complex_Series_Vectors.Vector;
               step : quad_double )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the vector of power series srv at the step,
  --   in double, double double, or quad double precision.

end Series_and_Predictors;
