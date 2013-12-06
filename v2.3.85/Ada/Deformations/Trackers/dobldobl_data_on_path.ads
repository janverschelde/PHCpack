with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with Double_Double_Vectors;
with Double_Double_VecVecs;
with DoblDobl_Continuation_Data;         use DoblDobl_Continuation_Data;
with DoblDobl_Complex_Solutions;         use DoblDobl_Complex_Solutions;
with Continuation_Parameters;            use Continuation_Parameters;

package DoblDobl_Data_on_Path is

-- DESCRIPTION :
--   This package provides data management operations for path tracking
--   in complex double double arithmetic.

-- UTILITIES :

  function At_End ( t,target : Complex_Number; distance,tol : double_float )
                  return boolean;

  -- DESCRIPTION :
  --   Decides whether at end of continuation stage.

  -- ON ENTRY :
  --   t        current value of the continuation parameter;
  --   target   target value for t;
  --   distance is the distance to invoke the end game;
  --   tol      numbers less than tol are considered zero.

  function Stop ( p : Pred_Pars; t,target : Complex_Number;
                  step : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if either the step size is smaller than p.minstep, or
  --   or alternatively, in case of geometric predictor, if the distance to
  --   the target has become smaller than p.minstep.

  -- ON ENTRY :
  --   p        settings of the predictor;
  --   t        current value of the continuation parameter;
  --   target   target value for t;
  --   step     current step size.

-- MANAGEMENT OF DATA DURING PATH FOLLOWING :

  procedure Linear_Single_Initialize 
              ( s : in Solu_Info; p : in Pred_Pars;
                old_t,prev_t : out Complex_Number;
                prev_v,old_solution,prev_solution
                  : out DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Initialization for linear path following of one path.

  -- ON ENTRY :
  --   s                solution at beginning of path;
  --   p                predictor parameters.

  -- ON RETURN :
  --   old_t            back up value for continuation parameter;
  --   prev_t           previous value of continuation parameter;
  --   old_solution     back up value for solution vector;
  --   prev_solution    previous value of solution vector;

  procedure Linear_Single_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t : in out Complex_Number;
                old_solution,prev_solution,old_v,prev_v,vv
                  : in out DoblDobl_Complex_Vectors.Vector;
                step : in out double_float;
                nsuccess,trial : in out natural32;
                success : in out boolean );

  -- DESCRIPTION :
  --   Management of data after correction during linear path following.

  -- PARAMETERS :
  --   s                current solution;
  --   p                predictor parameters;
  --   c                corrector parameters;
  --   old_t            back up value for continuation parameter;
  --   prev_t           previous value of continuation parameter;
  --   old_solution     back up value for solution vector;
  --   prev_solution    previous value of solution vector;
  --   old_v            back up value for tangent direction;
  --   prev_v           previous value for tangent direction;
  --   vv               current tangent direction;
  --   step             current step size;
  --   nsuccess         number of consecutive successes;
  --   trial            number of trials after failure;
  --   success          successful correction step.

  procedure Linear_Single_Quadratic_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t,prev_t0 : in out Complex_Number;
                old_s,prev_s,prev_s0 : in out DoblDobl_Complex_Vectors.Vector;
                step : in out double_float;
                nsuccess,trial : in out natural32;
                success : in out boolean );

  -- DESCRIPTION :
  --   Management of data after correction during linear path following
  --   for a quadratic predictor.

  -- PARAMETERS :
  --   s           current solution;
  --   p           predictor parameters;
  --   c           corrector parameters;
  --   old_t       back up value for continuation parameter;
  --   prev_t      previous value of continuation parameter;
  --   prev_t0     value of continuation parameter prior to prev_t;
  --   old_s       back up value for solution vector;
  --   prev_s      previous value of solution vector for t = prev_t;
  --   prev_s0     previous value of solution vector for t = prev_t0;
  --   step        current step size;
  --   nsuccess    number of consecutive successes;
  --   trial       number of trials after failure;
  --   success     successful correction step.

  procedure Linear_Single_Cubic_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t,prev_t1,prev_t0 : in out Complex_Number;
                old_s,prev_s,prev_s1,prev_s0
                  : in out DoblDobl_Complex_Vectors.Vector;
                step : in out double_float;
                nsuccess,trial : in out natural32;
                success : in out boolean );

  -- DESCRIPTION :
  --   Management of data after correction during linear path following
  --   for a cubic predictor.

  -- PARAMETERS :
  --   s           current solution;
  --   p           predictor parameters;
  --   c           corrector parameters;
  --   old_t       back up value for continuation parameter;
  --   prev_t      previous value of continuation parameter;
  --   prev_t1     value of continuation parameter prior to prev_t;
  --   prev_t0     value of continuation parameter prior to prev_t1;
  --   old_s       back up value for solution vector;
  --   prev_s      previous value of solution vector for t = prev_t;
  --   prev_s1     previous value of solution vector for t = prev_t1;
  --   prev_s0     previous value of solution vector for t = prev_t0;
  --   step        current step size;
  --   nsuccess    number of consecutive successes;
  --   trial       number of trials after failure;
  --   success     successful correction step.

  procedure Linear_Single_Quartic_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prv_t,prv_t2,prv_t1,prv_t0 : in out Complex_Number;
                old_s,prv_s,prv_s2,prv_s1,prv_s0
                  : in out DoblDobl_Complex_Vectors.Vector;
                step : in out double_float;
                nsuccess,trial : in out natural32;
                success : in out boolean );

  -- DESCRIPTION :
  --   Management of data after correction during linear path following
  --   for a quartic predictor.

  -- PARAMETERS :
  --   s           current solution;
  --   p           predictor parameters;
  --   c           corrector parameters;
  --   old_t       back up value for continuation parameter;
  --   prv_t       previous value of continuation parameter;
  --   prv_t2      value of continuation parameter prior to prv_t;
  --   prv_t1      value of continuation parameter prior to prv_t2;
  --   prv_t0      value of continuation parameter prior to prv_t1;
  --   old_s       back up value for solution vector;
  --   prv_s       previous value of solution vector for t = prv_t;
  --   prv_s2      previous value of solution vector for t = prv_t2;
  --   prv_s1      previous value of solution vector for t = prv_t1;
  --   prv_s0      previous value of solution vector for t = prv_t0;
  --   step        current step size;
  --   nsuccess    number of consecutive successes;
  --   trial       number of trials after failure;
  --   success     successful correction step.

  procedure Linear_Multiple_Initialize 
              ( s : in Solu_Info_Array; p : in Pred_Pars;
                t,old_t,prev_t : out Complex_Number;
                sa,old_sa,prev_sa : in out Solution_Array );

  -- DECRIPTION :
  --   Initialization for linear path following for more than one path.

  procedure Linear_Multiple_Management
              ( s : in out Solu_Info_array;
                sa,old_sa,prev_sa : in out Solution_Array;
                t,old_t,prev_t : in out Complex_Number; p : in Pred_Pars; 
                step : in out double_float; pivot : in integer32; 
                nsuccess,trial : in out natural32;
                success : in boolean );

  -- DESCRIPTION :
  --   Management of data after correction during linear path following.

  -- PARAMETERS :
  --   s            current solutions with information statistics;
  --   sa           current solutions;
  --   old_sa       back up value for solutions;
  --   prev_sa      previous solutions;
  --   t            current value of continuation parameter;
  --   old_t        back up value for continuation parameter;
  --   prev_t       previous value of continuation parameter;
  --   p            predictor parameters;
  --   step         current step size;
  --   pivot        solution where correction is started;
  --   nsuccess     number of consecutive successes;
  --   trial        number of trials after failure;
  --   success      successful correction step.

  procedure Circular_Management
              ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                old_t,prev_t : in out Complex_Number;
                start_solution : in DoblDobl_Complex_Vectors.Vector;
                old_solution,prev_solution,w_sum,w_all_sum
                  : in out DoblDobl_Complex_Vectors.Vector;
                twopi,epslop,tol : in double_float;
                theta,old_theta,step : in out double_float;
                nsuccess,n_sum,n_all_sum,w_c : in out natural32;
                max_wc : in natural32; stop,success : in out boolean );

  -- DESCRIPTION :
  --   Management of circular path following.

  -- PARAMETERS :
  --   s                current solution;
  --   p                predictor parameters;
  --   c                corrector parameters;
  --   old_t            back up value for continuation parameter;
  --   prev_t           previous value of continuation parameter;
  --   start_solution   solution vector at start of continuation;
  --   old_solution     back up value for solution vector;
  --   prev_solution    previous value of solution vector;
  --   w_sum            sum of cycle;
  --   w_all_sum        sum of all cycles;
  --   twopi            two times PI;
  --   epslop           tolerance to decide whether two vectors are equal;
  --   theta            current value of theta;
  --   old_theta        back up value for theta;
  --   step             current step size;
  --   nsuccess         number of consecutive successes;
  --   n_sum            number of cycles;
  --   n_all_sum        total number of cycles;
  --   w_c              winding number;
  --   max_wc           upper bound on winding number;
  --   stop             true when winding number has been computed;
  --   success          successful correction step.

-- UPDATE OF PATH DIRECTION :

  procedure Update_Direction
              ( proj : in boolean;
                freqcnt,defer : in out natural32; 
                r,m,estm : in out integer32; cntm : in out natural32;
                thresm : in natural32; er : in out integer32;
                t,target : in Complex_Number;
                x : in DoblDobl_Complex_Vectors.Vector; 
                dt,s,logs : in out Double_Double_Vectors.Vector;
                logx,wvl0,wvl1,wvl2 : in out Double_Double_VecVecs.VecVec;
                v,errv : in out Double_Double_Vectors.Vector;
                err : in out double_double );

  -- DESCRIPTION :
  --   Computes an approximation of the direction of the path.

  -- ON ENTRY :
  --   file       to write intermediate data on, may be omitted;
  --   proj       whether solution vector lies in projective space;
  --   freqcnt    counts the frequency of calls;
  --   defer      only if freqcnt = defer, calculations are done;
  --   r          current order in extrapolation formula;
  --   m          current value for multiplicity;
  --   estm       current estimated for multiplicity;
  --   cntm       number of consecutive times estm has been guessed;
  --   thresm     threshold for augmenting m to estm;
  --   er         order of extrapolator on the errors;
  --   t          current value of continuation parameter;
  --   target     target value of continuation parameter;
  --   x          current solution vector;
  --   dt         vector with distances to target;
  --   s          s-values w.r.t. the current value m;
  --   logs       logarithms of the s-values;
  --   logx       logarithms of the solution vectors;
  --   wvl0       consecutive estimates for previous direction;
  --   wvl1       consecutive estimates for current direction;
  --   wvl2       used as work space for wvl0 and wvl1;
  --   v          current approximate direction of the path;
  --   errv       vector of errors used to estimate m;
  --   err        norm of errv.

  -- ON RETURN :
  --   All in-out variables are updated, provided freqcnt = defer.

  procedure Update_Direction
              ( file : in file_type; proj : in boolean;
                freqcnt,defer : in out natural32;
                r,m,estm : in out integer32; cntm : in out natural32;
                thresm : in natural32; er : in out integer32;
                t,target : in Complex_Number;
                x : in DoblDobl_Complex_Vectors.Vector; 
                dt,s,logs : in out Double_Double_Vectors.Vector;
                logx,wvl0,wvl1,wvl2 : in out Double_Double_VecVecs.VecVec;
                v,errv : in out Double_Double_Vectors.Vector;
                err : in out double_double );

  -- DESCRIPTION :
  --   Computes an approximation of the direction of the path and produces
  --   intermediate output to the file.

end DoblDobl_Data_on_Path;
