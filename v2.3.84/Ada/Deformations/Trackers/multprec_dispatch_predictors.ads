with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;          use Multprec_Complex_Matrices;
with Multprec_Complex_Solutions;         use Multprec_Complex_Solutions;
with Multprec_Continuation_Data;         use Multprec_Continuation_Data;

package Multprec_Dispatch_Predictors is

-- DESCRIPTION :
--   This package provides generic predictors.  Based on the type of
--   predictor in Pred_Pars, the appropriate predictor is invoked.

  generic

    with function Norm ( x : Vector ) return Floating_Number;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Single_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars; xt : in boolean;
                prev_x,prev_v : in Vector; v : in out Vector;
                prev_t,target : in Complex_Number;
                step,tol : in Floating_Number; trial : in out natural32 );

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

  generic

    with function Norm ( x : Vector ) return Floating_Number;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Multiple_Predictor
              ( s : in out Solu_Info_Array; p : in Pred_Pars; xt : in boolean;
                sa : in out Solution_Array; prev_sa : in Solution_Array; 
                t : in out Complex_Number; prev_t,target : in Complex_Number;
                step,tol,dist : in Floating_Number; trial : in natural32 );

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

end Multprec_Dispatch_Predictors;
