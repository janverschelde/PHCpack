with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Laur_JacoMats;
with DoblDobl_Jacobian_Circuits;
with DoblDobl_Complex_Solutions;

package DoblDobl_Root_Refiners is

-- DESCRIPTION :
--   Provides root refinement in complex double double arithmetic.

  procedure DoblDobl_Newton_Step
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double );
  procedure DoblDobl_Newton_Step
              ( f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double );

  -- DESCRIPTION :
  --   Does one Newton step in double double complex arithmetic.

  -- ON ENTRY :
  --   f        evaluable form of a (Laurent) polynomial system;
  --   jf       Jacobian matrix of f;
  --   x        current approximate solution.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   err      norm of the update vector;
  --   rco      estimate for the inverse condition number;
  --   res      residual, norm of the function value.

  procedure DoblDobl_Newton_Step
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Jacobian_Circuits.Circuit;
                x : in out DoblDobl_Complex_Vectors.Vector;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double );

  -- DESCRIPTION :
  --   Does one Newton step using a circuit to evaluate and differentiate.

  -- ON ENTRY :
  --   f        evaluable form of a polynomial system;
  --   jf       Jacobian matrix of f, defined as a circuit;
  --   x        current approximate solution;
  --   wrk      work space allocated for evaluated monomials.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   err      norm of the update vector;
  --   rco      estimate for the inverse condition number;
  --   res      residual, norm of the function value.

  procedure Silent_Newton
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean );
  procedure Reporting_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean );
  procedure Silent_Newton
              ( f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean );
  procedure Reporting_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Applies Newton's method to refine an approximate root x of f.
  --   Stops when one conditions is satisfied:
  --   (1) numit >= max (reached maximum number of iterations),
  --   (2) x.err < epsxa (update factor to x is less than epsxa),
  --   (3) x.res < epsfa (residual smaller than epsfa).

  -- ON ENTRY :
  --   file     to write intermediate diagnostics at each step;
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

  procedure Silent_Newton
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Jacobian_Circuits.Circuit;
                x : in out DoblDobl_Complex_Solutions.Solution;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean );
  procedure Reporting_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Jacobian_Circuits.Circuit;
                x : in out DoblDobl_Complex_Solutions.Solution;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Applies Newton's method to refine an approximate root x of f.
  --   Stops when one conditions is satisfied:
  --   (1) numit >= max (reached maximum number of iterations),
  --   (2) x.err < epsxa (update factor to x is less than epsxa),
  --   (3) x.res < epsfa (residual smaller than epsfa).
  --   The reporting version writes information about the solution to file.

  -- ON ENTRY :
  --   file     for writing diagnostics about the Newton updates;
  --   f        evaluable form of a (Laurent) polynomial system;
  --   jf       Jacobian matrix of f, defined as a circuit.
  --   x        current approximate solution,
  --   wrk      work space for the evaluated monomials;
  --   epsxa    accuracy requirement on update factor;
  --   epsfa    accuracy requirement on residual;
  --   numit    number of iterations, must be zero on entry,
  --   max      maximum number of iterations allowed.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   numit    number of iterations spent on refining x;
  --   fail     true if spent max number of iterations
  --            and none of the accuracy requirements is met.

  procedure DoblDobl_Root_Refiner
              ( f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                s : in DoblDobl_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Refines the solution s of the system f with Jacobian matrix jf,
  --   applying three Newton steps.

  procedure DoblDobl_Root_Refiner
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Jacobian_Circuits.Circuit;
                s : in DoblDobl_Complex_Solutions.Link_to_Solution;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Refines the solution s of the system f with Jacobian matrix jf,
  --   defined as a circuit, applying three Newton steps.

  -- ON ENTRY :
  --   f        polynomial system in evaluable form;
  --   jf       Jacobian matrix defined as a circuit;
  --   s        pointer to a solution;
  --   wrk      work space for the evaluated monomials.

  -- ON RETURN :
  --   s        content where the pointer refers to is updated;
  --   wrk      modified work space for the evaluated monomials.

  procedure DoblDobl_Root_Refiner
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                s : in out DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Applies Newton's method to the solutions s of the system p.

  procedure DoblDobl_Root_Refiner
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Applies Newton's method to the solutions s of the system p,
  --   using the circuit representation for the Jacobian matrix.

  procedure Silent_Root_Refiner
               ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa : in double_double;
                 numit : in out natural32; max : in natural32 );
  procedure Silent_Root_Refiner
               ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa : in double_double;
                 numit : in out natural32; max : in natural32 );
  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa : in double_double;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean );
  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa : in double_double;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean );

  -- DESCRIPTION :
  --   Applies Newton's method to refine roots of p in s.
  --   Stops when one conditions is satisfied:
  --   (1) numit >= max (reached maximum number of iterations),
  --   and for solutions x in s:
  --   (2) x.err < epsxa (update factor to x is less than epsxa),
  --   (3) x.res < epsfa (residual smaller than epsfa).

  -- ON ENTRY :
  --   file     for writing intermediate output and diagnostics;
  --   p        a polynomial system;
  --   s        current approximate solutions;
  --   epsxa    accuracy requirement on update factor;
  --   epsfa    accuracy requirement on residual;
  --   numit    number of iterations, must be zero on entry,
  --   max      maximum number of iterations allowed;
  --   wout     if true, then information about each Newton update
  --            is written to file.

  -- ON RETURN :
  --   s        updated approximate solutions;
  --   numit    number of iterations spent on refining x;

end DoblDobl_Root_Refiners;
