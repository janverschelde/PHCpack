with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_SysFun;

package Homotopy_Newton_Steps is

-- DESCRIPTION :
--   Defines steps with Newton's method on the polynomial homotopies
--   stored in the packages {Standard,DoblDobl,QuadDobl}_Homotopy.

  procedure Standard_LU_Newton_Step
              ( nq : in integer32;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float );
  procedure DoblDobl_LU_Newton_Step
              ( nq : in integer32;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float );
  procedure QuadDobl_LU_Newton_Step
              ( nq : in integer32;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float );

  -- DESCRIPTION :
  --   Does one step with Newton's method, using LU factorization to
  --   solve the linear system in the computation of the update to x,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED :
  --   The proper homotopy data is initialized properly,
  --   in Standard_Homotopy, DoblDobl_Homotopy, or QuadDobl_Homotopy,
  --   respectively for double, double double, or quad double precision.

  -- ON ENTRY :
  --   nq       number of equations in the homotopy;
  --   t        current value of the continuation parameter;
  --   x        approximation for a point on a solution path at t.

  -- ON RETURN :
  --   x        corrected point after one Newton step;
  --   err      magnitude of the correction, estimates the forward error;
  --   rco      estimate for the condition number of the Jacobian matrix;
  --   res      magnitude of the residual, estimates the backward error.

  procedure Standard_LU_Newton_Step
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float );
  procedure DoblDobl_LU_Newton_Step
              ( abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float );
  procedure QuadDobl_LU_Newton_Step
              ( abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float );

  -- DESCRIPTION :
  --   Does one step with Newton's method, using LU factorization to
  --   solve the linear system in the computation of the update to x,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED :
  --   The proper homotopy data is initialized properly,
  --   in Standard_Homotopy, DoblDobl_Homotopy, or QuadDobl_Homotopy,
  --   respectively for double, double double, or quad double precision.

  -- ON ENTRY :
  --   abh      homotopy polynomials with absolute coefficients;
  --   t        current value of the continuation parameter;
  --   x        approximation for a point on a solution path at t.

  -- ON RETURN :
  --   x        corrected point after one Newton step;
  --   err      magnitude of the correction, estimates the forward error;
  --   rco      estimate for the condition number of the Jacobian matrix;
  --   res      mixed residual.

  procedure Standard_SVD_Newton_Step
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                t : in Standard_Complex_Numbers.Complex_Number;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float );

  -- DESCRIPTION :
  --   Does one step with Newton's method, using SVD factorization to
  --   solve the linear system in the computation of the update to x,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED :
  --   The proper homotopy data is initialized properly,
  --   in Standard_Homotopy, DoblDobl_Homotopy, or QuadDobl_Homotopy,
  --   respectively for double, double double, or quad double precision.

  -- ON ENTRY :
  --   abh      homotopy polynomials with absolute coefficients;
  --   t        current value of the continuation parameter;
  --   x        approximation for a point on a solution path at t.

  -- ON RETURN :
  --   x        corrected point after one Newton step;
  --   err      magnitude of the correction, estimates the forward error;
  --   rco      estimate for the condition number of the Jacobian matrix;
  --   res      mixed residual.

  procedure Correct
              ( nq : in integer32; t,tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 );
  procedure Correct
              ( nq : in integer32; t,tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 );
  procedure Correct
              ( nq : in integer32; t,tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 );

  -- DESCRIPTION :
  --   Applies Newton's method to correct the solution, silent version,
  --   in standard double, double double, or quad double precision.

  -- REQUIRED :
  --   Depending on the double, double double, or quad double precision,
  --   the corresponding Standard_Homotopy, DoblDobl_Homotopy, or
  --   QuadDobl_Homotopy must have been properly initialized.

  -- ON ENTRY :
  --   nq       number of equations in the homotopy;
  --   t        current value of the continuation parameter;
  --   tolres   tolerance on the residual, stops when res <= tolres;
  --   maxit    maximum number of steps done with Newton's method;
  --   sol      predicted value for the solution;
  --   extra    number of extra steps done as long as err goes down.

  -- ON RETURN :
  --   nbrit    number of iterations done;
  --   sol      the corrected value for the solution;
  --   err      magnitude of the correction term;
  --   rco      estimate for the inverse condition number;
  --   res      magnitude of the residual;
  --   fail     true if tolres was not met within maxit steps,
  --            and/or Newton's method diverged,
  --            false if Newton's method converged well.

  procedure Correct
              ( file : in file_type;
                nq : in integer32; t,tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false );
  procedure Correct
              ( file : in file_type;
                nq : in integer32; t,tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false );
  procedure Correct
              ( file : in file_type;
                nq : in integer32; t,tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Applies Newton's method to correct the solution, verbose version,
  --   in standard double or double double precision.

  -- REQUIRED :
  --   Depending on the double, double double, or quad double precision,
  --   the corresponding Standard_Homotopy, DoblDobl_Homotopy, or
  --   QuadDobl_Homotopy must have been properly initialized.

  -- ON ENTRY :
  --   file     for writing extra diagnostic output, if verbose;
  --   t        current value of the continuation parameter;
  --   tolres   tolerance on the residual, stops when res <= tolres;
  --   maxit    maximum number of steps to do with Newton's method;
  --   sol      predicted value for the solution;
  --   extra    number of extra steps done as long as err goes down.
  --   verbose  to indicate that extra output is wanted.

  -- ON RETURN :
  --   sol      corrected value of the solution;
  --   nbrit    number of iterations done;
  --   err      magnitude of the correction term;
  --   rco      estimate for the inverse condition number;
  --   res      magnitude of the residual;
  --   fail     true if tolres was not met within maxit steps,
  --            and/or Newton's method diverged,
  --            false if Newton's method converged well.

  procedure Correct
              ( abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float;
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 );
  procedure Correct
              ( abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 );
  procedure Correct
              ( abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0 );

  -- DESCRIPTION :
  --   Applies Newton's method to correct the solution, silent version,
  --   in standard double, double double, or quad double precision,
  --   with the computation of mixed residuals.

  -- REQUIRED :
  --   Depending on the double, double double, or quad double precision,
  --   the corresponding Standard_Homotopy, DoblDobl_Homotopy, or
  --   QuadDobl_Homotopy must have been properly initialized.

  -- ON ENTRY :
  --   abh      homotopy polynomials with absolute coefficients;
  --   t        current value of the continuation parameter;
  --   tolres   tolerance on the residual, stops when res <= tolres;
  --   maxit    maximum number of steps done with Newton's method;
  --   extra    number of extra steps done as long as err goes down.
  --   sol      predicted value for the solution.

  -- ON RETURN :
  --   nbrit    number of iterations done;
  --   sol      the corrected value for the solution;
  --   err      magnitude of the correction term;
  --   rco      estimate for the inverse condition number;
  --   res      the mixed residual of the solution;
  --   fail     true if tolres was not met within maxit steps,
  --            and/or Newton's method diverged,
  --            false if Newton's method converged well.

  procedure Correct
              ( file : in file_type;
                abh : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false );
  procedure Correct
              ( file : in file_type;
                abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false );
  procedure Correct
              ( file : in file_type;
                abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                t,tolres : in double_float; 
                maxit : in natural32; nbrit : out natural32;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out double_float; fail : out boolean;
                extra : in natural32 := 0;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Applies Newton's method to correct the solution, verbose version,
  --   in standard double or double double precision,
  --   with the computation of mixed residuals.

  -- REQUIRED :
  --   Depending on the double, double double, or quad double precision,
  --   the corresponding Standard_Homotopy, DoblDobl_Homotopy, or
  --   QuadDobl_Homotopy must have been properly initialized.

  -- ON ENTRY :
  --   file     for writing extra diagnostic output, if verbose;
  --   abh      homotopy polynomials with absolute coefficients;
  --   t        current value of the continuation parameter;
  --   tolres   tolerance on the residual, stops when res <= tolres;
  --   maxit    maximum number of steps to do with Newton's method;
  --   sol      predicted value for the solution;
  --   extra    number of extra steps done as long as err goes down.
  --   verbose  to indicate that extra output is wanted.

  -- ON RETURN :
  --   sol      corrected value of the solution;
  --   nbrit    number of iterations done;
  --   err      magnitude of the correction term;
  --   rco      estimate for the inverse condition number;
  --   res      magnitude of the residual;
  --   fail     true if tolres was not met within maxit steps,
  --            and/or Newton's method diverged,
  --            false if Newton's method converged well.

end Homotopy_Newton_Steps;
