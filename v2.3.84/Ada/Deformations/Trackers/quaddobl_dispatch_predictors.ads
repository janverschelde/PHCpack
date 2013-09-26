with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;
with Continuation_Parameters;            use Continuation_Parameters;
with QuadDobl_Continuation_Data;         use QuadDobl_Continuation_Data;

package QuadDobl_Dispatch_Predictors is

-- DESCRIPTION :
--   This package provides generic predictors, based on the type of
--   predictor in Pred_Pars, the appropriate predictor is invoked.

  generic

    with function Norm ( x : Vector ) return quad_double;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Single_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars; xt : in boolean;
                prev_x,prev_v : in Vector; v : in out Vector;
                prev_t,target : in Complex_Number;
                step,tol : in double_float; trial : in out natural32 );

  -- DESCRIPTION :
  --   Generic predictor for one solution.

  -- ON ENTRY :
  --   s        information about the current solution;
  --   p        parameters for the predictor;
  --   xt       if true then prediction for both x and t is needed,
  --            if false then no prediction for x, only for t;
  --   prev_x   previous solution component (only for secant);
  --   prev_t   previous value for t (only useful for secant);
  --   target   target value for continuation parameter;
  --   step     current step size;
  --   tol      tolerance for floating equalities;
  --   trial    number of consecutive trials (for complex predictor).

  -- ON RETURN :
  --   s        predicted value for solution.

  function Real_Equal ( a,b : Complex_Number; tol : in double_float )
                      return boolean;

  -- DESCRIPTION :
  --   Returns true if the real parts of a and b are within distance tol
  --   from each other.  This is needed as it appears that because of
  --   roundoff accumulation in summing various values of t, the t could
  --   become very close to a target value.

  procedure Single_Quadratic_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars; xt : in boolean;
                x1,x0 : in Vector;
                t1,t0,target : in Complex_Number;
                step,tol : in double_float );

  -- DESCRIPTION :
  --   Quadratic predictor for one solution.

  -- ON ENTRY :
  --   s        information about the current solution;
  --   p        parameters for the predictor;
  --   xt       if true then prediction for both x and t is needed,
  --            if false then no prediction for x, only for t;
  --   x1       the approximation for t = t1;
  --   x0       the approximation for t = t0;
  --   t        the current value of the continuation parameter;
  --   t1       is the previous value for t;
  --   t0       is the value for t, prior to t1;
  --   target   target value for continuation parameter;
  --   step     current step size;
  --   tol      tolerance for floating equalities.

  -- ON RETURN :
  --   s        predicted value for solution.

  procedure Single_Cubic_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars; xt : in boolean;
                x2,x1,x0 : in Vector;
                t2,t1,t0,target : in Complex_Number;
                step,tol : in double_float );

  -- DESCRIPTION :
  --   Quadratic predictor for one solution.

  -- ON ENTRY :
  --   s        information about the current solution;
  --   p        parameters for the predictor;
  --   xt       if true then prediction for both x and t is needed,
  --            if false then no prediction for x, only for t;
  --   x2       the approximation for t = t2;
  --   x1       the approximation for t = t1;
  --   x0       the approximation for t = t0;
  --   t        the current value of the continuation parameter;
  --   t2       is the previous value for t;
  --   t1       is the value for t, prior to t2;
  --   t0       is the value for t, prior to t1;
  --   target   target value for continuation parameter;
  --   step     current step size;
  --   tol      tolerance for floating equalities.

  -- ON RETURN :
  --   s        predicted value for solution.

  procedure Single_Quartic_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars; xt : in boolean;
                x3,x2,x1,x0 : in Vector;
                t3,t2,t1,t0,target : in Complex_Number;
                step,tol : in double_float );

  -- DESCRIPTION :
  --   Quadratic predictor for one solution.

  -- ON ENTRY :
  --   s        information about the current solution;
  --   p        parameters for the predictor;
  --   xt       if true then prediction for both x and t is needed,
  --            if false then no prediction for x, only for t;
  --   x3       the approximation for t = t3;
  --   x2       the approximation for t = t2;
  --   x1       the approximation for t = t1;
  --   x0       the approximation for t = t0;
  --   t        the current value of the continuation parameter;
  --   t3       is the previous value for t;
  --   t2       is the value for t, prior to t3;
  --   t1       is the value for t, prior to t2;
  --   t0       is the value for t, prior to t1;
  --   target   target value for continuation parameter;
  --   step     current step size;
  --   tol      tolerance for floating equalities.

  -- ON RETURN :
  --   s        predicted value for solution.

  generic

    with function Norm ( x : Vector ) return quad_double;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Multiple_Predictor
              ( s : in out Solu_Info_Array; p : in Pred_Pars; xt : in boolean;
                sa : in out Solution_Array; prev_sa : in Solution_Array; 
                t : in out Complex_Number; prev_t,target : in Complex_Number;
                step,tol,dist : in double_float; trial : in natural32 );

  -- DESCRIPTION :
  --   Generic predictor for an array of solutions.

  -- ON ENTRY :
  --   s        array with information of current solutions;
  --   p        parameters for the predictor;
  --   xt       if true then prediction for both x and t is needed,
  --            if false then no prediction for x, only for t; 
  --   sa       the current solutions;
  --   prev_sa  previous solution component (only for secant);
  --   t        current value for continuation parameter;
  --   prev_t   previous value for t (only useful for secant);
  --   target   target value for continuation parameter;
  --   step     current step size;
  --   tol      tolerance for floating equalities;
  --   dist     tolerance for distance between solutions;
  --   trial    number of consecutive trials (for complex predictor).

  -- ON RETURN :
  --   sa       predicted values for solutions;
  --   t        predicted continuation parameter.

end QuadDobl_Dispatch_Predictors;
