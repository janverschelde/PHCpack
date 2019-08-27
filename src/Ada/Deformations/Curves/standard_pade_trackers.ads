with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;  
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Solutions;
with Standard_Complex_Series_VecVecs;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Hessians;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_SysFun;
with Standard_CSeries_Jaco_Matrices;
with Standard_Pade_Approximants;
with Homotopy_Continuation_Parameters;

package Standard_Pade_Trackers is

-- DESCRIPTION :
--   Path trackers which compute Pade approximants in standard precision
--   are provided, in silent or verbose versions.

  function Maximum ( a,b : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the maximum of a and b.

  function Minimum ( a,b : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the minimum of a and b.

  procedure Set_Step
              ( t,step : in out double_float;
                maxstep,target : in double_float );

  -- DESCRIPTION :
  --   Sets the step size and the new value for t,
  --   taking into account the maximal step size maxstep
  --   and the value for the target for t.

  procedure Update_Step_Sizes
              ( minsize,maxsize : in out double_float;
                step : in double_float );

  -- DESCRIPTION :
  --   Given in step is a size of a step,
  --   the smallest and largest sizes in minsize and maxsize are updated.

  function Residual_Prediction
              ( sol : Standard_Complex_Vectors.Vector;
                t : double_float ) return double_float;

  -- DESCRIPTION :
  --   Given a predicted solution and a value for t,
  --   returns the norm of the evaluated homotopy at the solution.

  -- REQUIRED :
  --   The Standard_Homotopy must have been properly initialized.

  function Residual_Prediction
              ( abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                sol : Standard_Complex_Vectors.Vector;
                t : double_float ) return double_float;
  function Residual_Prediction
              ( file : file_type;
                abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                sol : Standard_Complex_Vectors.Vector;
                t : double_float ) return double_float;

  -- DESCRIPTION :
  --   Given a predicted solution and a value for t,
  --   returns the mixed residual of the solution.

  -- REQUIRED :
  --   The Standard_Homotopy must have been properly initialized.

  -- ON ENTRY :
  --   file     for extra output on the residual computation;
  --   abh      homotopy polynomials with absolute value coefficients;
  --   sol      current solution for the value of t;
  --   t        current value of the continuation parameter.

  procedure Predictor_Feedback
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                sol : out Standard_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep : in double_float;
                cntcut : in out natural32 );
  procedure Predictor_Feedback
              ( file : in file_type;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                sol : out Standard_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep : in double_float;
                cntcut : in out natural32 );
  procedure Predictor_Feedback
              ( file : in file_type; verbose : in boolean;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                sol : out Standard_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep : in double_float;
                cntcut : in out natural32 );

  -- DESCRIPTION :
  --   Runs the predictor feedback loop, cutting the step size in half
  --   if the residual of the predicted solution is too large.

  -- ON ENTRY :
  --   file     to write extra output and diagnostics;
  --   verbose  flag indicates if extra output will be written to file;
  --   abh      homotopy with absolute value coefficients;
  --   pv       vector of Pade approximants;
  --   t        current value of the continuation parameter;
  --   step     current value of the step size;
  --   tolpres  tolerance on the predictor residual;
  --   minstep  minimum step size;
  --   cntcut   current count of cuts in the step size.

  -- ON RETURN :
  --   sol      predicted solution as evaluated at the Pade vector;
  --   predres  the residual of the predicted solution;
  --   t        updated value of the continuation parameter;
  --   step     updated value of the step size;
  --   cntcut   updated value of the counter of step size cuts.

  procedure Predictor_Corrector
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                sol : out Standard_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep,tolcres : in double_float;
                maxit,extra : in natural32; nbrcorrs : in out natural32;
                err,rco,res : out double_float;
                cntcut,cntfail : in out natural32; fail : out boolean );
  procedure Predictor_Corrector
              ( file : in file_type; verbose : in boolean;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                sol : out Standard_Complex_Vectors.Vector;
                predres : out double_float;
                t,step : in out double_float;
                tolpres,minstep,tolcres : in double_float;
                maxit,extra : in natural32; nbrcorrs : in out natural32;
                err,rco,res : out double_float;
                cntcut,cntfail : in out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Runs the predictor-corrector feedback loop.

  -- ON ENTRY :
  --   file     to write extra output and diagnostics;
  --   verbose  flag indicates if extra output will be written to file;
  --   abh      homotopy with absolute value coefficients;
  --   pv       vector of Pade approximants;
  --   t        current value of the continuation parameter;
  --   step     current value of the step size;
  --   tolpres  tolerance on the predictor residual;
  --   minstep  minimum step size;
  --   tolcres  tolerance on the corrector residual;
  --   maxit    maximum number of corrector iterations;
  --   extra    number of extra corrector iterations;
  --   nbrcorrs is current number of corrector iterations;
  --   cntcut   current count of cuts in the step size;
  --   cntfail  current count of corrector failures.

  -- ON RETURN :
  --   sol      predicted solution as evaluated at the Pade vector;
  --   predres  the residual of the predicted solution;
  --   t        updated value of the continuation parameter;
  --   step     updated value of the step size;
  --   err      forward error;
  --   rco      estimate for the inverse condition number;
  --   res      residual of the solution;
  --   cntcut   updated value of the counter of step size cuts;
  --   cntfail  updated value of the counter of the corrector failures;
  --   fail     true if failed to meet the tolerance tolcres.

  procedure Minimum_Step_Size
              ( sstp,dstp,pstp : in double_float; minstp : out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32 );
  procedure Minimum_Step_Size
              ( file : in file_type;
                sstp,dstp,pstp : in double_float; minstp : out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32 );

  -- DESCRIPTION :
  --   Determines the smallest step and updates the corresponding counter.

  -- ON ENTRY :
  --   file     output file for writing diagnostics;
  --   sstp     step size determined by the power series;
  --   dstp     step size determined by the curvature;
  --   pstp     step size determined by the closest pole;
  --   cntsstp  counts the number of times sstp was smallest;
  --   cntdstp  counts the number of times dstp was smallest;
  --   cntpstp  counts the number of times pstp was smallest.

  -- ON RETURN :
  --   minstp   the smallest value of sstp, dstp, and pstp;
  --   cntsstp  updated count of the number of times sstp was smallest;
  --   cntdstp  updated count of the number of times dstp was smallest;
  --   cntpstp  updated count of the number of times pstp was smallest.

  procedure Step_Control
              ( jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out Standard_Pade_Approximants.Pade_Vector;
                poles : in out Standard_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 );
  procedure Step_Control
              ( file : in file_type; verbose : in boolean;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out Standard_Pade_Approximants.Pade_Vector;
                poles : in out Standard_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Determines the step size for the next step on a path.

  -- ON ENTRY :
  --   file     to write extra output and diagnostics;
  --   verbose  flag indicates if extra output will be written to file;
  --   jm       Jacobian matrix for the polynomial homotopy;
  --   hs       Hessians for all equations in the polynomial homotopy;
  --   hom      a homotopy with series coefficients;
  --   sol      current solution on a path;
  --   maxdeg   largest degree of the power series;
  --   nit      number of iterations allowed in Newton for power series;
  --   pars     values of the parameters and tolerances;
  --   pv       space allocated for a vector of Pade approximants;
  --   poles    space allocated for the poles of the Pade approximants;
  --   t        current value of the continuation parameter;
  --   step     current value of the step size;
  --   cntsstp  counts the number of times sstp was smallest;
  --   cntdstp  counts the number of times dstp was smallest;
  --   cntpstp  counts the number of times pstp was smallest;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   pv       vector of Pade approximants;
  --   poles    poles of the Pade approximants;
  --   t        updated value of the continuation parameter;
  --   step     updated value of the step size;
  --   cntsstp  updated count of the number of times sstp was smallest;
  --   cntdstp  updated count of the number of times dstp was smallest;
  --   cntpstp  updated count of the number of times pstp was smallest.

  procedure Step_Control
              ( jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in Standard_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out Standard_Pade_Approximants.Pade_Vector;
                poles : in out Standard_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 );
  procedure Step_Control
              ( file : in file_type; verbose : in boolean;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in Standard_Complex_Vectors.Vector;
                maxdeg,nit : in integer32;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                pv : in out Standard_Pade_Approximants.Pade_Vector;
                poles : in out Standard_Complex_VecVecs.VecVec;
                t,step : in out double_float;
                cntsstp,cntdstp,cntpstp : in out natural32;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Determines the step size for the next step on a path,
  --   for a coefficient-parameter homotopy.

  -- ON ENTRY :
  --   file     to write extra output and diagnostics;
  --   verbose  flag indicates if extra output will be written to file;
  --   jm       Jacobian matrix for the polynomial homotopy;
  --   hs       Hessians for all equations in the polynomial homotopy;
  --   fhm      coefficient-parameter homotopy for efficient evaluation,
  --            the series parameter is the continuation parameter;
  --   fcf      coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   sol      current solution on a path;
  --   maxdeg   largest degree of the power series;
  --   nit      number of iterations allowed in Newton for power series;
  --   pars     values of the parameters and tolerances;
  --   pv       space allocated for a vector of Pade approximants;
  --   poles    space allocated for the poles of the Pade approximants;
  --   t        current value of the continuation parameter;
  --   step     current value of the step size;
  --   cntsstp  counts the number of times sstp was smallest;
  --   cntdstp  counts the number of times dstp was smallest;
  --   cntpstp  counts the number of times pstp was smallest;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   pv       vector of Pade approximants;
  --   poles    poles of the Pade approximants;
  --   t        updated value of the continuation parameter;
  --   step     updated value of the step size;
  --   cntsstp  updated count of the number of times sstp was smallest;
  --   cntdstp  updated count of the number of times dstp was smallest;
  --   cntpstp  updated count of the number of times pstp was smallest.

  procedure Track_One_Path
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntcut,cntfail : out natural32;
                minsize,maxsize : out double_float;
                cntsstp,cntdstp,cntpstp : out natural32;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution sol using the homotopy hom.
  --   This version remains silent and does not write any output.

  -- ON ENTRY :
  --   abh      homotopy with absolute value coefficients;
  --   jm       Jacobian matrix for the polynomial homotopy;
  --   hs       Hessians for all equations in the polynomial homotopy;
  --   hom      a homotopy with series coefficients;
  --   sol      start solution in the homotopy;
  --   pars     values of the parameters and tolerances;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   sol      solution at the end of the path;
  --   nbrsteps is the total number of steps on the path;
  --   nbrcorrs is the total number of corrector iterations on the path;
  --   cntcut   is the total number of steps cut by predictor residual;
  --   cntfail  is the total number of corrector failures on the path;
  --   minsize  is the smallest step size on the path;
  --   maxsize  is the largest step size on the path;
  --   cntsstp  counts the number of times sstp was smallest;
  --   cntdstp  counts the number of times dstp was smallest;
  --   cntpstp  counts the number of times pstp was smallest.

  procedure Track_One_Path
              ( file : in file_type;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                nbrsteps,nbrcorrs,cntcut,cntfail : out natural32;
                minsize,maxsize : out double_float;
                cntsstp,cntdstp,cntpstp : out natural32;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks one path starting at the solution sol using the homotopy hom.
  --   This version is verbose and writes extra diagnostic output to file.

  -- ON ENTRY :
  --   file     for writing output during the computations;
  --   abh      homotopy with absolute value coefficients;
  --   jm       Jacobian matrix for the polynomial homotopy;
  --   hs       Hessians for all equations in the polynomial homotopy;
  --   hom      a homotopy with series coefficients;
  --   sol      start solution in the homotopy;
  --   pars     values of the parameters and tolerances;
  --   verbose  flag indicates if extra output will be written to file;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   sol      solution at the end of the path;
  --   nbrsteps is the total number of steps on the path;
  --   nbrcorrs is the total number of corrector iterations on the path;
  --   cntcut   is the total number of steps cut by predictor residual;
  --   cntfail  is the total number of corrector failes on the paths;
  --   minsize  is the smallest step size on the path;
  --   maxsize  is the largest step size on the path;
  --   cntsstp  counts the number of times sstp was smallest;
  --   cntdstp  counts the number of times dstp was smallest;
  --   cntpstp  counts the number of times pstp was smallest.

-- VERSIONS WITH COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure Track_One_Path
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat; 
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                nbrsteps,nbrcorrs,cntcut,cntfail : out natural32;
                minsize,maxsize : out double_float;
                cntsstp,cntdstp,cntpstp : out natural32;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks one path starting using a coefficient-parameter homotopy.
  --   This version remains silent and does not write any output.

  -- ON ENTRY :
  --   abh      homotopy with absolute value coefficients;
  --   jm       Jacobian matrix for the polynomial homotopy;
  --   hs       Hessians for all equations in the polynomial homotopy;
  --   fhm      coefficient-parameter homotopy for efficient evaluation,
  --            the series parameter is the continuation parameter;
  --   fcf      coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   sol      start solution in the homotopy;
  --   pars     values of the parameters and tolerances;
  --   mhom     0 for affine, 1 for 1-homogenization, m for m-homogenization;
  --   idz      index representation of the partition z, for mhom > 1;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   sol      solution at the end of the path;
  --   nbrsteps is the total number of steps on the path;
  --   nbrcorrs is the total number of corrector iterations on the path;
  --   cntcut   is the total number of steps cut by predictor residual;
  --   cntfail  is the total number of corrector failures on the path;
  --   minsize  is the smallest step size on the path;
  --   maxsize  is the largest step size on the path;
  --   cntsstp  counts the number of times sstp was smallest;
  --   cntdstp  counts the number of times dstp was smallest;
  --   cntpstp  counts the number of times pstp was smallest.

  procedure Track_One_Path
              ( file : in file_type;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
                hs : in Standard_Complex_Hessians.Link_to_Array_of_Hessians;
                fhm : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                fcf : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                sol : in out Standard_Complex_Solutions.Solution;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in natural32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                nbrsteps,nbrcorrs,cntcut,cntfail : out natural32;
                minsize,maxsize : out double_float;
                cntsstp,cntdstp,cntpstp : out natural32;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks one path starting using a coefficient-parameter homotopy.
  --   This version is verbose and writes extra diagnostic output to file.

  -- ON ENTRY :
  --   file     for writing output during the computations;
  --   abh      homotopy with absolute value coefficients;
  --   jm       Jacobian matrix for the polynomial homotopy;
  --   hs       Hessians for all equations in the polynomial homotopy;
  --   fhm      coefficient-parameter homotopy for efficient evaluation,
  --            the series parameter is the continuation parameter;
  --   fcf      coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   sol      start solution in the homotopy;
  --   pars     values of the parameters and tolerances;
  --   mhom     0 for affine, 1 for 1-homogenization, m for m-homogenization;
  --   idz      index representation of the partition z, for mhom > 1;
  --   verbose  flag indicates if extra output will be written to file;
  --   vrblvl   the verbose level.

  -- ON RETURN :
  --   sol      solution at the end of the path;
  --   nbrsteps is the total number of steps on the path;
  --   nbrcorrs is the total number of corrector iterations on the path;
  --   cntcut   is the total number of steps cut by predictor residual;
  --   cntfail  is the total number of corrector failes on the paths;
  --   minsize  is the smallest step size on the path;
  --   maxsize  is the largest step size on the path;
  --   cntsstp  counts the number of times sstp was smallest;
  --   cntdstp  counts the number of times dstp was smallest;
  --   cntpstp  counts the number of times pstp was smallest.

end Standard_Pade_Trackers;
