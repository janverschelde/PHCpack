with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;

package QuadDobl_Predictors is

-- DESCRIPTION :
--   This package contains several implementations for the predictor 
--   in an increment-and-fix continuation.

--   The predictor provides a prediction both for the continuation parameter t
--   and for the solution(s) x.

--   For the continuation paramter t the following options can be made :
--     Real      : linear prediction, simply adds the step size;
--     Complex   : can make predictions in complex space;
--     Circular  : to perform a circular sample, for winding numbers;
--     Geometric : distances to target form geometric series.

--   For the solution vector x the following options are provided :
--     Secant  : linear extrapolation using differences;
--     Tangent : linear extrapolation using the first derivatives;
--     Hermite : third-order extrapolation using first derivatives.
--   Furthermore, these predictors for x can be applied
--   for one solution (Single) or for an array of solutions (Multiple).

--   By combining these options, the following 13 predictors are provided :

--     Secant_Single_Real_Predictor
--     Secant_Single_Complex_Predictor
--     Secant_Multiple_Real_Predictor
--     Secant_Multiple_Complex_Predictor
--     Tangent_Single_Real_Predictor
--     Tangent_Single_Complex_Predictor
--     Tangent_Multiple_Real_Predictor
--     Tangent_Multiple_Complex_Predictor

--     Secant_Circular_Predictor
--     Secant_Geometric_Predictor
--     Tangent_Circular_Predictor
--     Tangent_Geometric_Predictor

--     Hermite_Single_Real_Predictor

-- The order in which these predictors are listed depends on their mutual
-- resemblances of the specified parameters.

-- PREDICTORS for t only :

  procedure Real_Predictor
                ( t : in out Complex_Number; target : in Complex_Number;
                  h,tol : in double_float; pow : in natural32 := 1;
                  hh : out double_float );

  -- DESCRIPTION :
  --   Linear prediction for t with real increment.

  -- ON ENTRY :
  --   t          current value for the continuation parameter;
  --   target     final target value for t, kind of upper bound on t;
  --   h          step size;
  --   tol        tolerance to decide whether zero;
  --   pow        to deal with power prediction.

  -- ON RETURN :
  --   t          predicted value for t;
  --   hh         step size which is different from h if t would have
  --              had otherwise a larger real part than the target.

  procedure Complex_Predictor
                ( t : in out Complex_Number; target : in Complex_Number;
                  h,tol : in double_float; hh : out double_float;
                  distance : in double_float; trial : in natural32 );

  -- DESCRIPTION :
  --   Takes complex values, when trial /= 0 and stays at the given
  --   distance from the target.

  procedure Circular_Predictor
               ( t : in out Complex_Number; theta : in out double_float;
                 t0_min_target,target : in Complex_Number;
                 h : in double_float );

  -- DESCRIPTION :
  --   Predictor used to wind around the target value.

   procedure Geometric_Predictor
               ( t : in out Complex_Number; target : in Complex_Number;
                 h,tol : in double_float );

   -- DESCRIPTION :
   --   Geometric progression towards the target, no danger of overshooting.

-- PREDICTORS for x and t :

  procedure Secant_Single_Real_Predictor 
                ( x : in out Vector; prev_x : in Vector;
                  t : in out Complex_Number; prev_t,target : in Complex_Number;
                  h,tol : in double_float; pow : in natural32 := 1 );

  procedure Secant_Multiple_Real_Predictor 
                ( x : in out Solution_Array; prev_x : in Solution_Array;
                  t : in out Complex_Number; prev_t,target : in Complex_Number;
                  h,tol,dist_x : in double_float; pow : in natural32 := 1 );

  -- DESCRIPTION :
  --   Secant predictor for x and a linear predictor for t.

  -- ON ENTRY :
  --   x          the current approximation(s) for t;
  --   prev_x     the approximation(s) for a previous t;
  --   t          the current value of the continuation parameter;
  --   prev_t     is the previous value for t;
  --   target     is the target value for t;
  --   h          is the step size;
  --   tol        tolerance to decide when t = target;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > d, for k in 1..n;
  --   pow        power of t in the homotopy towards 1.

  -- ON RETURN :
  --   x          the predicted approximation(s);
  --   t          the predicted value of the continuation parameter.

  procedure Quadratic_Single_Real_Predictor
                ( x : in out Vector; x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t1,t0,target : in Complex_Number;
                  h,tol : in double_float; pow : in natural32 := 1 );

  -- DESCRIPTION :
  --   Quadratic predictor for x and a linear predictor for t.

  -- REQUIRED : t > t1 > t0.

  -- ON ENTRY :
  --   x          the current approximation for t;
  --   x1         the approximation for t = t1;
  --   x0         the approximation for t = t0;
  --   t          the current value of the continuation parameter;
  --   t1         is the previous value for t;
  --   t0         is the value for t, prior to t1;
  --   target     is the target value for t;
  --   h          is the step size;
  --   tol        tolerance to decide when t = target;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > d, for k in 1..n;
  --   pow        power of t in the homotopy towards 1.

  -- ON RETURN :
  --   x          the predicted approximation;
  --   t          the predicted value of the continuation parameter.

  procedure Cubic_Single_Real_Predictor
                ( x : in out Vector; x2,x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t2,t1,t0,target : in Complex_Number;
                  h,tol : in double_float; pow : in natural32 := 1 );

  -- DESCRIPTION :
  --   Cubic predictor for x and a linear predictor for t.

  -- REQUIRED : t > t2 > t1 > t0.

  -- ON ENTRY :
  --   x          the current approximation for t;
  --   x2         the approximation for t = t2;
  --   x1         the approximation for t = t1;
  --   x0         the approximation for t = t0;
  --   t          the current value of the continuation parameter;
  --   t2         is the previous value for t;
  --   t1         is the value for t, prior to t2;
  --   t0         is the value for t, prior to t1;
  --   target     is the target value for t;
  --   h          is the step size;
  --   tol        tolerance to decide when t = target;
  --   pow        power of t in the homotopy towards 1.

  -- ON RETURN :
  --   x          the predicted approximation;
  --   t          the predicted value of the continuation parameter.

  procedure Quartic_Single_Real_Predictor
                ( x : in out Vector; x3,x2,x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t3,t2,t1,t0,target : in Complex_Number;
                  h,tol : in double_float; pow : in natural32 := 1 );

  -- DESCRIPTION :
  --   Quartic predictor for x and a linear predictor for t.

  -- REQUIRED : t > t3 > t2 > t1 > t0.

  -- ON ENTRY :
  --   x          the current approximation for t;
  --   x3         the approximation for t = t3;
  --   x2         the approximation for t = t2;
  --   x1         the approximation for t = t1;
  --   x0         the approximation for t = t0;
  --   t          the current value of the continuation parameter;
  --   t3         is the previous value for t;
  --   t2         is the value for t, prior to t3;
  --   t1         is the value for t, prior to t2;
  --   t0         is the value for t, prior to t1;
  --   target     is the target value for t;
  --   h          is the step size;
  --   tol        tolerance to decide when t = target;
  --   pow        power of t in the homotopy towards 1.

  -- ON RETURN :
  --   x          the predicted approximation;
  --   t          the predicted value of the continuation parameter.

  procedure Quintic_Single_Real_Predictor
                ( x : in out Vector; x4,x3,x2,x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t4,t3,t2,t1,t0,target : in Complex_Number;
                  h,tol : in double_float; pow : in natural32 := 1 );

  -- DESCRIPTION :
  --   Quintic predictor for x and a linear predictor for t.

  -- REQUIRED : t > t4 > t3 > t2 > t1 > t0.

  -- ON ENTRY :
  --   x          the current approximation for t;
  --   x4         the approximation for t = t4;
  --   x3         the approximation for t = t3;
  --   x2         the approximation for t = t2;
  --   x1         the approximation for t = t1;
  --   x0         the approximation for t = t0;
  --   t          the current value of the continuation parameter;
  --   t4         is the previous value for t;
  --   t3         is the value for t, prior to t4;
  --   t2         is the value for t, prior to t3;
  --   t1         is the value for t, prior to t2;
  --   t0         is the value for t, prior to t1;
  --   target     is the target value for t;
  --   h          is the step size;
  --   tol        tolerance to decide when t = target;
  --   pow        power of t in the homotopy towards 1.

  -- ON RETURN :
  --   x          the predicted approximation;
  --   t          the predicted value of the continuation parameter.

  procedure Secant_Single_Complex_Predictor 
                ( x : in out Vector; prev_x : in Vector;
                  t : in out Complex_Number; prev_t,target : in Complex_Number;
                  h,tol,dist_t : in double_float; trial : in natural32 );

  procedure Secant_Multiple_Complex_Predictor 
                ( x : in out Solution_Array; prev_x : in Solution_Array;
                  t : in out Complex_Number; prev_t,target : in Complex_Number; 
                  h,tol,dist_x,dist_t : in double_float; 
                  trial : in natural32 );

  -- DESCRIPTION :
  --   Secant predictor for x and complex predictor for t.

  -- ON ENTRY :
  --   x          the current approximation(s) for t;
  --   prev_x     the approximation(s) for a previous t;
  --   t          the current value of the continuation parameter;
  --   prev_t     is the previous value for t;
  --   target     is the target value for t;
  --   h          is the step size;
  --   tol        tolerance to decide when two numbers are equal;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > d, for k in 1..n;
  --   dist_t     t must keep a distance to the target;
  --   trial      indicates the number of trials for starting out of
  --              the previous value for t.

  -- ON RETURN :
  --   x          the predicted approximation(s);
  --   t          the predicted value of the continuation parameter.

  procedure Secant_Circular_Predictor
                ( x : in out Vector; prev_x : in Vector;
                  t : in out Complex_Number; theta : in out double_float;
                  prev_t,t0_min_target,target : in Complex_Number;
                  h,tol : in double_float );

  -- DESCRIPTION :
  --   Secant predictor for x and circular predictor for t, around target.

  -- NOTE : This is the link between t and theta :
  --   t = target + t0_min_target * ( cos(theta) + i sin(theta) )

  -- ON ENTRY :
  --   x,prev_x,t,prev_t,target as before;
  --   theta      the angle for t.
  --   t0_min_target is t0-target, where t0 is the start point;
  --   h          is step size for theta !!

  -- ON RETURN :
  --   x          the predicted approximation;
  --   t          the predicted value of the continuation parameter;
  --   theta      the predicted angle.

  procedure Secant_Geometric_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out Complex_Number; prev_t,target : in Complex_Number;
                 h,tol : in double_float );

  -- DESCRIPTION :
  --   Secant predictor for x and a geometric predictor for t.

  -- ON ENTRY :
  --   x          the current approximation(s) for t;
  --   prev_x     the approximation(s) for a previous t;
  --   t          the current value of the continuation parameter;
  --   prev_t     is the previous value for t;
  --   target     is the target value for t;
  --   h          ratio between two consecutive distance to target, 0<h<1;
  --   tol        tolerance to decide when t = target;

  -- ON RETURN :
  --   x          the predicted approximation;
  --   t          the predicted value of the continuation parameter.

  generic

    with function Norm ( x : Vector ) return quad_double;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Single_Real_Predictor 
                ( x : in out Vector; t : in out Complex_Number;
                  target : in Complex_Number; h,tol : in double_float;
                  pow : in natural32 := 1 );

  generic

    with function Norm ( x : Vector ) return quad_double;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Multiple_Real_Predictor
                ( x : in out Solution_Array; t : in out Complex_Number;
                  target : in Complex_Number; h,tol,dist_x : in double_float;
                  nsys : in out natural32; pow : in natural32 := 1 );

  -- DESCRIPTION :
  --   Tangent predictor for x and a linear predictor for t.

  -- ON ENTRY :
  --   x          current approximation for the solution;
  --   t          current value of the continuation parameter;
  --   target     target value for the continuation parameter;
  --   h          steplength;
  --   tol        tolerance to decide when t = target;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > dist_x, for k in 1..n;
  --   nsys       must be initally equal to zero, used for counting;
  --   pow        power of t in the homotopy.

  -- ON RETURN :
  --   x          predicted approximation for the solution;
  --   t          new value of the continuation parameter;
  --   nsys       the number of linear systems solved.

  generic

    with function Norm ( x : Vector ) return quad_double;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Single_Complex_Predictor
                ( x : in out Vector; t : in out Complex_Number; 
                  target : in Complex_Number;
                  h,tol,dist_t : in double_float; trial : in natural32 );

  generic

    with function Norm ( x : Vector) return quad_double;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Multiple_Complex_Predictor
                ( x : in out Solution_Array; t : in out Complex_Number;
                  target : in Complex_Number;
                  h,tol,dist_x,dist_t : in double_float;
                  trial : in natural32; nsys : in out natural32 );

  -- DESCRIPTION :
  --   Tangent predictor for x and a complex predictor for t.

  -- ON ENTRY :
  --   x          current approximation for the solution;
  --   t          current value of the continuation parameter;
  --   target     target value for the continuation parameter;
  --   h          steplength;
  --   tol        tolerance to decide when two numbers are equal;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > d, for k in 1..n;
  --   dist_t     t must keep distance to the target;
  --   trial      indicates the number of trials for starting out of
  --              the previous value for t;
  --   nsys       must be initially equal to zero.

  -- ON RETURN :
  --   x          predicted approximation for the solution;
  --   t          new value of the continuation parameter;
  --   nsys       the number of linear systems solved.

  generic

    with function Norm ( x : Vector) return quad_double;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Circular_Predictor
                ( x : in out Vector; t : in out Complex_Number;
                  theta : in out double_float;
                  t0_min_target,target : in Complex_Number;
                  h,tol : in double_float );

  -- DESCRIPTION :
  --   This is a tangent predictor for x and a circular predictor for t
  --   For information on the parameters, see Secant_Circular_Predictor.

  generic

    with function Norm ( x : Vector) return quad_double;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Geometric_Predictor
               ( x : in out Vector; t : in out Complex_Number;
                 target : in Complex_Number; h,tol : in double_float );

  -- DESCRIPTION :
  --   Tangent predictor for x and a geometric predictor for t.
  --   For information on the parameters, see Secant_Geometric_Predictor.

  generic

    with function Norm ( x : Vector) return quad_double;
    with function dH ( x : Vector; t : Complex_Number ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : Complex_Number ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Hermite_Single_Real_Predictor
                ( x : in out Vector; prev_x : in Vector;
                  t : in out Complex_Number; prev_t,target : in Complex_Number;
                  v : in out Vector; prev_v : in Vector;
                  h,tol : in double_float; pow : in natural32 := 1 );

  -- DESCRIPTION :
  --   Third-order extrapolation based on previous values of the solution
  --   paths along with corresponding first derivatives.

 -- ON ENTRY :
  --   x          current approximation for the solution;
  --   prev_x     previous approximation of the solution at prev_t;
  --   t          current value of the continuation parameter;
  --   prev_t     previous value of the continuation parameter;
  --   target     target value for the continuation parameter;
  --   v          will be used at work space;
  --   prev_v     direction of the path at prev_t;
  --   h          steplength;
  --   tol        tolerance to decide when t = target;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > dist_x, for k in 1..n;
  --   pow        power of t in the homotopy.

  -- ON RETURN :
  --   x          predicted approximation for the solution;
  --   t          new value of the continuation parameter;
  --   v          direction of the path computed for value of t on entry.

end QuadDobl_Predictors;
