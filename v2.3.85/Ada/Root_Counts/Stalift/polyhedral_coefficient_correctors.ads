with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;      use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;    use Standard_Complex_Laur_Jacomats;

package Polyhedral_Coefficient_Correctors is

-- DESCRIPTION :
--   This package offers various implementations of Newton's method,
--   to serve as a corrector in a polyhedral coefficient homotopy.

  procedure Silent_Newton_Step
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y : out double_float );

  procedure Silent_Newton_Step
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y,rcond : out double_float );

  procedure Reporting_Newton_Step
              ( file : in file_type; hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y : out double_float );

  procedure Reporting_Newton_Step
              ( file : in file_type; hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                x,y : in out Standard_Complex_Vectors.Vector;
                norm_dx,norm_y,rcond : out double_float );

  -- DESCRIPTION :
  --   Applies one Newton step for a fixed value of the continuation
  --   parameter, so the coefficients ctm are fixed as well. 

  -- ON ENTRY :
  --   file     for printing norm of increment vector and residual;
  --   hq       evaluable form of coefficient Laurent system;
  --   ctm      coefficients for fixed value of continuation parameter
  --   jacmat   coefficient Jacobian matrix;
  --   mulfac   multiplication factors for Jacobian matrix;
  --   x        current approximation for the solution;
  --   y        solution evaluated at x.

  -- ON RETURN :
  --   x        updated solution vector;
  --   y        system evaluated at the new x;
  --   norm_dx  norm of the increment added to x;
  --   norm_y   norm of the vector y;
  --   rcond    estimate for the inverse of the condition number.

  procedure Silent_Apply_Newton
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                tol : in double_float; max : in natural32;
                x,y : in out Standard_Complex_Vectors.Vector;
                nit : out natural32; fail : out boolean );

  procedure Reporting_Apply_Newton
              ( file : in file_type; hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                tol : in double_float; max : in natural32;
                x,y : in out Standard_Complex_Vectors.Vector;
                nit : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Applies Newton's method until the tolerance is reached,
  --   or until the maximul number allowed steps is exhausted,
  --   in which case failure is reported.

  -- ON ENTRY :
  --   file     for intermediate output;
  --   hq       evaluable form of coefficient Laurent system;
  --   ctm      coefficients for fixed value of continuation parameter
  --   jacmat   coefficient Jacobian matrix;
  --   mulfac   multiplication factors for Jacobian matrix;
  --   tol      the iteration stops if the increment vector or
  --            the residual is smaller in norm than tol;
  --   max      maximal number of iterations allowed to reach
  --            the desired accuracy;
  --   x        current approximation for the solution;
  --   y        solution evaluated at x.

  -- ON RETURN :
  --   x        updated solution vector;
  --   y        system evaluated at the new x.
  --   nit      number of iterations performed;
  --   fail     true if the required accuracy was not reached
  --            within max Newton steps, false otherwise.

end Polyhedral_Coefficient_Correctors;
