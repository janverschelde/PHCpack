with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_JacoMats;
with QuadDobl_Complex_Solutions;

package QuadDobl_Root_Refiners is

-- DESCRIPTION :
--   Provides root refinement in complex quad double arithmetic.

  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double );
  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double );

  -- DESCRIPTION :
  --   Does one Newton step in quad double complex arithmetic.

  -- ON ENTRY :
  --   f        evaluable form of a (Laurent) polynomial system;
  --   jf       Jacobian matrix of f;
  --   x        current approximate solution.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   err      norm of the update vector;
  --   rco      estimate for the inverse condition number;
  --   res      residual, norm of the function value.

  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in quad_double; numit : in out natural32;
                max : in natural32; fail : out boolean );
  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in quad_double; numit : in out natural32;
                max : in natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Applies Newton's method to refine an approximate root x of f.
  --   Stops when one conditions is satisfied:
  --   (1) numit >= max (reached maximum number of iterations),
  --   (2) x.err < epsxa (update factor to x is less than epsxa),
  --   (3) x.res < epsfa (residual smaller than epsfa).

  -- ON ENTRY :
  --   f        evaluable form of a (Laurent) polynomial system;
  --   jf       Jacobian matrix of f;
  --   x        current approximate solution,
  --   epsxa    accuracy requirement on update factor;
  --   epsfa    accuracy requirement on residual;
  --   numit    number of iterations, must be zero on entry,
  --   max      maximum number of iterations allowed.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   numit    number of iterations spent on refining x;
  --   fail     true if spent max number of iterations
  --            and none of the accuracy requirements is met.

  procedure QuadDobl_Root_Refiner
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                s : in QuadDobl_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Refines the solution s of the system f with Jacobi matrix jf.

  procedure QuadDobl_Root_Refiner
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Applies Newton's method to the solutions s of the system p.

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa : in quad_double;
                 numit : in out natural32; max : in natural32 );

  -- DESCRIPTION :
  --   Applies Newton's method to refine roots of p in s.
  --   Stops when one conditions is satisfied:
  --   (1) numit >= max (reached maximum number of iterations),
  --   and for solutions x in s:
  --   (2) x.err < epsxa (update factor to x is less than epsxa),
  --   (3) x.res < epsfa (residual smaller than epsfa).

  -- ON ENTRY :
  --   p        a polynomial system;
  --   s        current approximate solutions;
  --   epsxa    accuracy requirement on update factor;
  --   epsfa    accuracy requirement on residual;
  --   numit    number of iterations, must be zero on entry,
  --   max      maximum number of iterations allowed.

  -- ON RETURN :
  --   s        updated approximate solutions;
  --   numit    number of iterations spent on refining x;

end QuadDobl_Root_Refiners;
