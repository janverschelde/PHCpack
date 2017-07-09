with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;
with Multprec_Complex_Poly_Systems;      use Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;       use Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;     use Multprec_Complex_Jaco_Matrices;
with Multprec_Complex_Laur_Systems;
with Multprec_Complex_Laur_SysFun;       use Multprec_Complex_Laur_SysFun;
with Multprec_Complex_Laur_JacoMats;
with Multprec_Complex_Solutions;         use Multprec_Complex_Solutions;

package Multprec_Root_Refiners is

-- DESCRIPTION:
--   The aim of this package is to provide some routines to
--   perform root refining, to validate the approximate solutions.
--   It can be used as a postprocessor to check the computed roots,
--   or, as preprocessor, to find some suitable starting values for
--   the continuation procedure.

--   The basic root refiner is Newton's method.  There are four versions
--   provided : silent and reporting, each combined with optional estimator
--   of multiplicities.  The reporting version puts intermediate results on
--   file during the iterations, while the silent version simply returns 
--   the refined solutions.  

-- ONE NEWTON STEP :

  procedure Multprec_Newton_Step
              ( f : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out Multprec_Complex_Vectors.Vector;
                err,rco,res : out Floating_Number );
  procedure Multprec_Newton_Step
              ( f : in Multprec_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in Multprec_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out Multprec_Complex_Vectors.Vector;
                err,rco,res : out Floating_Number );

  -- DESCRIPTION :
  --   Does one Newton step in multiprecision complex arithmetic.

  -- ON ENTRY :
  --   f        evaluable form of a polynomial system;
  --   jf       Jacobian matrix of f;
  --   x        current approximate solution.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   err      norm of the update vector;
  --   rco      estimate for the inverse condition number;
  --   res      residual, norm of the function value.

-- NEWTON's METHOD :

  procedure Silent_Newton 
               ( p_eval : in Eval_Poly_Sys;
                 j_eval : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in Floating_Number; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean );
  procedure Silent_Newton 
               ( p_eval : in Eval_Laur_Sys; 
                 j_eval : in Multprec_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in Floating_Number; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean );

  procedure Reporting_Newton
               ( file : in file_type; p_eval : in Eval_Poly_Sys;
                 j_eval : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 fail : out boolean );
  procedure Reporting_Newton
               ( file : in file_type; p_eval : in Eval_Laur_Sys;
                 j_eval : in Multprec_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Newton's method is applied to refine the approximation of a root.
  --   The stopping criteria are:
  --     * numit > max   (maximum number of iterations is exceeded);
  --     * zero.err < epsxa (accuracy for x is reached);
  --     * zero.res < epsfa (tolerance for residual is reached).
  --   When one of these conditions is fulfilled, the procedure stops.

  -- ON ENTRY :
  --   file      to write intermediate results on;
  --   p_eval    evaluable form of the polynomial system;
  --   j_eval    evaluable form of the Jacobian matrix;
  --   zero      starting value;
  --   epsxa     maximum absolute error on the zero;
  --   epsfa     maximum absolute value of the residue;
  --   numit     number of iterations, to be initialized with zero;
  --   max       maximum number of iterations.

  -- ON RETURN :
  --   zero      refined root;
  --   orders    estimated orders of multiplicities, if estord was on;
  --   numit     number of iterations performed;
  --   fail      is true when the desired precision is not reached.

-- FOR OVERDETERMINED SYSTEMS :

  procedure Silent_Gauss_Newton
               ( p_eval : in Eval_Poly_Sys; j_eval : in Eval_Jaco_Mat;
                 zero : in out Solution; tol : in double_float;
                 epsxa,epsfa : in Floating_Number; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean );

  procedure Reporting_Gauss_Newton
               ( file : in file_type;
                 p_eval : in Eval_Poly_Sys; j_eval : in Eval_Jaco_Mat;
                 zero : in out Solution; tol : in double_float;
                 epsxa,epsfa : in Floating_Number; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Applies Gauss-Newton to an overconstrained system.

  -- ON ENTRY :
  --   file      to write intermediate results on;
  --   p_eval    evaluable form of the polynomial system;
  --   j_eval    evaluable form of the Jacobian matrix;
  --   zero      starting value;
  --   tol       tolerance for the numerical rank;
  --   epsxa     maximum absolute error on the zero;
  --   epsfa     maximum absolute value of the residue;
  --   numit     number of iterations, to be initialized with zero;
  --   max       maximum number of iterations.

  -- ON RETURN :
  --   zero      refined root;
  --   orders    estimated orders of multiplicities, if estord was on;
  --   numit     number of iterations performed;
  --   fail      is true when the desired precision is not reached.

-- APPLICATION OF NEWTON's METHOD TO LISTS OF SOLUTIONS.
--   The silent versions simply perform the calculations.  
--   The reporting root refiners allow the output of intermediate results
--   and produce a summary of the calculations.
--   With each version, an additional output parameter can be supplied to
--   contain only those solutions that satisfy the accuracy requirements.
  
  procedure Silent_Root_Refiner
               ( p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate : in boolean );

  procedure Silent_Root_Refiner
               ( p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate : in boolean );

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate,wout : in boolean );

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate,wout : in boolean );

  -- DESCRIPTION :
  --   The list of solutions sols is refined w.r.t. the system p.
  --   The multiplicity of each solution in sols is determined as follows:
  --     m = 0 : if the solution is singular and probably non isolated
  --             or if the solution lies at infinity ( in fact no solution );
  --     m = 1 : if the solution is regular;
  --     m > 1 : a multiple solution with multiplicity m.

  -- ON ENTRY :
  --   file      file for writing diagnostics on;
  --   p         a polynomial system;
  --   sols      the start solutions;
  --   epsxa     maximum absolute error on the zero;
  --   epsfa     maximum absolute value for the residue;
  --   tolsing   tolerance on inverse condition number for singular solution;
  --   numit     the number of iterations, to be initialized with zero;
  --   max       maximum number of iterations per zero;
  --   deflate   if true, deflate singular solutions; 
  --   wout      has to be true when intermediate output is wanted.

  -- ON RETURN :
  --   sols      a list of computed solutions;
  --   refsols   only those solutions which satisfy the given accuracy;
  --   numit     the number of iterations.

-- APPLICATION OF GAUSS-NEWTON's METHOD TO LISTS OF SOLUTIONS.
--   The silent versions simply perform the calculations.  
--   The reporting root refiners allow the output of intermediate results
--   and produce a summary of the calculations.
--   With each version, an additional output parameter can be supplied to
--   contain only those solutions that satisfy the accuracy requirements.
  
  procedure Silent_Root_Sharpener
               ( p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate : in boolean );

  procedure Silent_Root_Sharpener
               ( p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate : in boolean );

  procedure Reporting_Root_Sharpener
               ( file : in file_type;
                 p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate,wout : in boolean );

  procedure Reporting_Root_Sharpener
               ( file : in file_type;
                 p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate,wout : in boolean );

  -- DESCRIPTION :
  --   The list of solutions sols is refined w.r.t. the system p.
  --   The multiplicity of each solution in sols is determined as follows:
  --     m = 0 : if the solution is singular and probably non isolated
  --             or if the solution lies at infinity ( in fact no solution );
  --     m = 1 : if the solution is regular;
  --     m > 1 : a multiple solution with multiplicity m.

  -- ON ENTRY :
  --   file      file for writing diagnostics on;
  --   p         a polynomial system;
  --   sols      the start solutions;
  --   epsxa     maximum absolute error on the zero;
  --   epsfa     maximum absolute value for the residue;
  --   tolsing   tolerance on inverse condition number for singular solution;
  --   numit     the number of iterations, to be initialized with zero;
  --   max       maximum number of iterations per zero;
  --   deflate   if true, deflate singular solutions; 
  --   wout      has to be true when intermediate output is wanted.

  -- ON RETURN :
  --   sols      a list of computed solutions;
  --   refsols   only those solutions which satisfy the given accuracy;
  --   numit     the number of iterations.

end Multprec_Root_Refiners;
