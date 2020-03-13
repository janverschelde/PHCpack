with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Hessians;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_SysFun;
with Standard_CSeries_Jaco_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Hessians;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_SysFun;
with DoblDobl_CSeries_Jaco_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Hessians;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_SysFun;
with QuadDobl_CSeries_Jaco_Matrices;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Pade_Approximants;
with DoblDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants;

package Series_and_Predictors is

-- DESCRIPTION :
--   Evaluates power series to predict a solution.

  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Vectors.Vector;
                srv : out Standard_Complex_Series_Vectors.Vector;
                eva : out Standard_Complex_Series_Vectors.Vector;
                vrblvl : in integer32 := 0 );
  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector;
                srv : out DoblDobl_Complex_Series_Vectors.Vector;
                eva : out DoblDobl_Complex_Series_Vectors.Vector;
                vrblvl : in integer32 := 0 );
  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Vectors.Vector;
                srv : out QuadDobl_Complex_Series_Vectors.Vector;
                eva : out QuadDobl_Complex_Series_Vectors.Vector;
                vrblvl : in integer32 := 0 );

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
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector;
                srv : out DoblDobl_Complex_Series_Vectors.Vector;
                eva : out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Vectors.Vector;
                srv : out QuadDobl_Complex_Series_Vectors.Vector;
                eva : out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

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

-- ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in Standard_Complex_Vectors.Vector;
                srv : out Standard_Complex_Series_Vectors.Vector;
                eva : out Standard_Complex_Series_Vectors.Vector;
                vrblvl : in integer32 := 0 );
  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                fhm : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in DoblDobl_Complex_Vectors.Vector;
                srv : out DoblDobl_Complex_Series_Vectors.Vector;
                eva : out DoblDobl_Complex_Series_Vectors.Vector;
                vrblvl : in integer32 := 0 );
  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                fhm : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in QuadDobl_Complex_Vectors.Vector;
                srv : out QuadDobl_Complex_Series_Vectors.Vector;
                eva : out QuadDobl_Complex_Series_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies Newton's method on the homotopy in hom,
  --   starting at the coordinates of the solution as the leading
  --   coefficients in the power series,
  --   in double, double double, or quad double precision.
  --   The versions are silent, do not write extra output.

  -- ON ENTRY :
  --   maxdeg   maximal degree of the series;
  --   nit      number of iterations with Newton's method;
  --   fhm      coefficient-parameter homotopy for evaluation,
  --            the series parameter is the continuation parameter;
  --   fcf      coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all partial derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   sol      solution of the homotopy for t = 0, its coordinates
  --            contain the leading coefficients of the power series.

  -- ON RETURN :
  --   srv      vector of power series, to predict the solution
  --            for some small positive value of the continuation parameter,
  --            srv'range = sol'range;
  --   eva      evaluated series vector srv in the homotopy fhm;
  --            eva'range = fhm'range.

  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in Standard_Complex_Vectors.Vector;
                srv : out Standard_Complex_Series_Vectors.Vector;
                eva : out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                fhm : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in DoblDobl_Complex_Vectors.Vector;
                srv : out DoblDobl_Complex_Series_Vectors.Vector;
                eva : out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                fhm : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in QuadDobl_Complex_Vectors.Vector;
                srv : out QuadDobl_Complex_Series_Vectors.Vector;
                eva : out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies Newton's method on the homotopy in hom,
  --   starting at the coordinates of the solution as the leading
  --   coefficients in the power series 
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     to write extra diagnostic output to;
  --   maxdeg   maximal degree of the series;
  --   nit      number of iterations with Newton's method;
  --   fhm      coefficient-parameter homotopy for evaluation,
  --            the series parameter is the continuation parameter;
  --   fcf      coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all partial derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   sol      solution of the homotopy for t = 0, its coordinates
  --            contain the leading coefficients of the power series;
  --   verbose  if true, then additional output is written to file.

  -- ON RETURN :
  --   srv      vector of power series, to predict the solution
  --            for some small positive value of the continuation parameter,
  --            srv'range = sol'range;
  --   eva      evaluated series vector srv in the homotopy fhm;
  --            eva'range = fhm'range.

  procedure Pade_Approximants
              ( srv : in Standard_Complex_Series_Vectors.Vector;
                pv : in out Standard_Pade_Approximants.Pade_Vector;
                poles : in out Standard_Complex_VecVecs.VecVec;
                frp : out double_float;
                cfp : out Standard_Complex_Numbers.Complex_Number;
                verbose : in boolean := false );
  procedure Pade_Approximants
              ( srv : in DoblDobl_Complex_Series_Vectors.Vector;
                pv : in out DoblDobl_Pade_Approximants.Pade_Vector;
                poles : in out DoblDobl_Complex_VecVecs.VecVec;
                frp : out double_double;
                cfp : out DoblDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false );
  procedure Pade_Approximants
              ( srv : in QuadDobl_Complex_Series_Vectors.Vector;
                pv : in out QuadDobl_Pade_Approximants.Pade_Vector;
                poles : in out QuadDobl_Complex_VecVecs.VecVec;
                frp : out quad_double;
                cfp : out QuadDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Given a power series vector, constructs Pade approximants,
  --   computes its poles and the smallest forward pole radius,
  --   in double, double double, or quad double precision.

  -- REQUIRED : the degree of the series in srv is large enough
  --   for the sum of numdeg and dendeg: srv.deg >= numdeg+dendeg.
  --   The numdeg and dendeg are defined in the allocated Pade vector.

  -- ON ENTRY :
  --   srv      vector power series, truncated to a certain degree;
  --   pv       space allocated for Pade vector of numdeg and dendeg;
  --   poles    space allocated to hold all poles;
  --   verbose  if verbose, then extra information is written to screen
  --            during the construction of the Pade approximants.

  -- ON RETURN :
  --   pv       vector of Pade approximants
  --   poles    poles of the Pade approximants;
  --   frp      radius of the closest pole;
  --   cfp      closest pole.

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

  function Cap_Step_Size
             ( step,frp,factor : double_float ) return double_float;

  -- DESCRIPTION :
  --   Caps the step size with the smallest pole radius frp,
  --   multiplied by a factor, which is typically smaller than one.
  --   The formula is min(step,frp*factor).

  function Step_Distance
             ( k : integer32; beta,eta,errnrm : double_float )
             return double_float;
  function Step_Distance
             ( k : integer32; beta,eta,errnrm : double_double )
             return double_double;
  function Step_Distance
             ( k : integer32; beta,eta,errnrm : quad_double )
             return quad_double;

  -- DESCRIPTION :
  --   Returns the k-th root of the ratio (beta*eta)/errnrm,
  --   where k is the maximal degree L+M+2, beta is some small factor,
  --   eta is the estimate for the distance computed with
  --   singular values of Jacobian and Hessians, and
  --   errnrm is the estimated solution error norm.
  --   The returned value is an estimate for the step size,
  --   based on the estimate to the nearest solution.

  function Step_Distance
            ( k : integer32; beta,tval : double_float;
              jm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : Standard_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : Standard_Complex_Vectors.Vector;
              srv : Standard_Complex_Series_Vectors.Vector;
              pv : Standard_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;
  function Step_Distance
            ( file : in file_type;
              k : integer32; beta,tval : double_float;
              jm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : Standard_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : Standard_Complex_Vectors.Vector;
              srv : Standard_Complex_Series_Vectors.Vector;
              pv : Standard_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;
  function Step_Distance
            ( k : integer32; beta : double_float; tval : double_double;
              jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : DoblDobl_Complex_Vectors.Vector;
              srv : DoblDobl_Complex_Series_Vectors.Vector;
              pv : DoblDobl_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;
  function Step_Distance
            ( file : in file_type;
              k : integer32; beta : double_float; tval : double_double;
              jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : DoblDobl_Complex_Vectors.Vector;
              srv : DoblDobl_Complex_Series_Vectors.Vector;
              pv : DoblDobl_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;
  function Step_Distance
            ( k : integer32; beta : double_float; tval : quad_double;
              jm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : QuadDobl_Complex_Vectors.Vector;
              srv : QuadDobl_Complex_Series_Vectors.Vector;
              pv : QuadDobl_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;
  function Step_Distance
            ( file : in file_type;
              k : integer32; beta : double_float; tval : quad_double;
              jm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : QuadDobl_Complex_Vectors.Vector;
              srv : QuadDobl_Complex_Series_Vectors.Vector;
              pv : QuadDobl_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;

  -- DESCRIPTION :
  --   Returns the step size set by the distance to the nearest solution,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file   optional file for output (if verbose);
  --   k      is the maximal degree L+M+2;
  --   beta   some small factor;
  --   tval   value of the the continuation parameter,
  --          corresponding to the current solution sol;
  --   jm     Jacobian matrix of the polynomial homotopy;
  --   hs     Hessians of the polynomials in the homotopy;
  --   sol    the current solution;
  --   srv    series approximation at the current solution;
  --   pv     vector of Pade approximants;
  --   verbose is the verbose flag.

  function Step_Distance
            ( k : integer32; beta : double_float;
              jm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : Standard_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : Standard_Complex_Solutions.Solution;
              srv : Standard_Complex_Series_Vectors.Vector;
              pv : Standard_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;
  function Step_Distance
            ( file : in file_type; k : integer32; beta : double_float;
              jm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : Standard_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : Standard_Complex_Solutions.Solution;
              srv : Standard_Complex_Series_Vectors.Vector;
              pv : Standard_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;
  function Step_Distance
            ( k : integer32; beta : double_float;
              jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : DoblDobl_Complex_Solutions.Solution;
              srv : DoblDobl_Complex_Series_Vectors.Vector;
              pv : DoblDobl_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;
  function Step_Distance
            ( file : in file_type; k : integer32; beta : double_float;
              jm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : DoblDobl_Complex_Solutions.Solution;
              srv : DoblDobl_Complex_Series_Vectors.Vector;
              pv : DoblDobl_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;
  function Step_Distance
            ( k : integer32; beta : double_float;
              jm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : QuadDobl_Complex_Solutions.Solution;
              srv : QuadDobl_Complex_Series_Vectors.Vector;
              pv : QuadDobl_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;
  function Step_Distance
            ( file : in file_type; k : integer32; beta : double_float;
              jm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
              hs : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
              sol : QuadDobl_Complex_Solutions.Solution;
              srv : QuadDobl_Complex_Series_Vectors.Vector;
              pv : QuadDobl_Pade_Approximants.Pade_Vector;
              verbose : boolean := false ) return double_float;

  -- DESCRIPTION :
  --   Returns the step size set by the distance to the nearest solution,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file   optional file for output (if verbose);
  --   k      is the maximal degree L+M+2;
  --   beta   some small factor;
  --   jm     Jacobian matrix of the polynomial homotopy;
  --   hs     Hessians of the polynomials in the homotopy;
  --   sol    the current solution;
  --   srv    series approximation at the current solution;
  --   pv     vector of Pade approximants;
  --   verbose is the verbose flag.

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

  function Predicted_Solution
             ( pv : Standard_Pade_Approximants.Pade_Vector;
               step : double_float )
             return Standard_Complex_Vectors.Vector;
  function Predicted_Solution
             ( pv : DoblDobl_Pade_Approximants.Pade_Vector;
               step : double_double )
             return DoblDobl_Complex_Vectors.Vector;
  function Predicted_Solution
             ( pv : QuadDobl_Pade_Approximants.Pade_Vector;
               step : quad_double )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the vector of Pade approximants pv at the step,
  --   in double, double double, or quad double precision.

end Series_and_Predictors;
