with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
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
with DoblDobl_Complex_Solutions;         use DoblDobl_Complex_Solutions;
with DoblDobl_Point_Lists;               use DoblDobl_Point_Lists;

package DoblDobl_Root_Refiners is

-- DESCRIPTION :
--   Provides root refinement in complex double double arithmetic.

-- ROOT ACCOUNTING :

  procedure Write_Info
              ( file : in file_type; zero : in Solution;
                initres : in double_double;
                i,numb,nbdef : in natural32;
                fail,infty : in boolean );

  -- DESCRIPTION :
  --   The information concerning the zero is written to file.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   zero     a solution;
  --   initres  the initial residual
  --   i        number of solution;
  --   numb     number of iterations;
  --   nbdef    number of deflations;
  --   fail     if failed or not;
  --   infty    at infinity or not.

  procedure Multiplicity
              ( h1,h2 : in DoblDobl_Complex_Vectors.Vector;
                pl : in out Point_List; ls : in Link_to_Solution;
                nb : in natural32; sols : in out Solution_List;
                fail,infty,deflate : in boolean;
                tolsing,tolclus : in double_float );

  -- DESCRIPTION :
  --   Checks whether the solution ls at position nb in sa is clustered.

  -- ON ENTRY :
  --   h1       first hash function for solution vectors;
  --   h2       second hash function for solution vectors;
  --   pl       list of hashed points;
  --   ls       current solution;
  --   nb       position of the solution ls in sa;
  --   sols     list of solutions;
  --   fail     if Newton failed to converge;
  --   infty    if at infinity;
  --   deflate  if deflation will be applied;
  --   tolsing  tolerance on singular solution;
  --   tolclus  tolerance for two solutions to be clustered.

  -- ON RETURN :
  --   pl       updated list of hashed points;
  --   ls       ls.m will reflect if multiple or clustered:
  --            if clustered, then ls.m < 0 and refers to the
  --            other solution in sols that is equal to ls.

  procedure Write_Type
              ( file : in file_type; ls : in Link_to_Solution;
                fail,infty : in boolean;
                tolsing : in double_float;
                nbfail,nbinfty : in out natural32;
                nbreal,nbcomp,nbreg,nbsing : in out natural32 );

  -- DESCRIPTION :
  --   Writes the type of the solution to file without cluster analysis.

  -- ON ENTRY :
  --   file     file opened for output;
  --   ls       the current solution;
  --   fail     if Newton failed to converge;
  --   infty    if at infinity;
  --   tolsing  tolerance for singular solutions;
  --   nbfail   current number of failures;
  --   nbinfty  current number of at infinity;
  --   nbreal   current number of real solutions;
  --   nbcomp   current number of complex solutions;
  --   nbreg    current number of regular solutions;
  --   nbsing   current number of singular solutions.

  -- ON RETURN :
  --   nbfail   updated number of failures;
  --   nbinfty  updated number of at infinity;
  --   nbreal   updated number of real solutions;
  --   nbcomp   updated number of complex solutions;
  --   nbreg    updated number of regular solutions;
  --   nbsing   updated number of singular solutions.

  procedure Write_Global_Info
              ( file : in file_type; tot,nbfail,nbinfty,
                nbreal,nbcomp,nbreg,nbsing,nbclus : in natural32 );

  -- DESCRIPTION :
  --   Writes the global information summary at the end of the solutions.

  -- ON ENTRY :
  --   file     file openen for output;
  --   tot      total number of solutions;
  --   nbfail   number of failures;
  --   nbinfty  number of solutions at infinity;
  --   nbreal   number of real solutions;
  --   nbcomp   number of complex solutions;
  --   nbreg    number of regular solutions;
  --   nbsing   number of singular solutions;
  --   nbclus   number of clustered solutions.

-- ONE NEWTON STEP :

  procedure Write_Diagnostics
              ( file : in file_type; step : natural32;
                err,rco,res : in double_double );

  -- DESCRIPTION :
  --   Writes diagnostics in one line about the Newton step to file.
  --   This procedure defines the formatting in the Reporting
  --   versions of a sequence of Newton steps.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   step     the current step number;
  --   err      forward error, magnitude of the update;
  --   rco      estimate for the inverse of the condition number;
  --   res      backward error, residual.

  procedure DoblDobl_SVD_Newton_Step
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );
  procedure DoblDobl_SVD_Newton_Step
              ( f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );
  procedure DoblDobl_SVD_Newton_Step
              ( f,abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does one Newton step in double double complex arithmetic,
  --   using the Singular Value Decomposition to compute the update to x
  --   and for the inverse condition number rco.
  --   Applies the method of Gauss-Newton, valid for f'last >= x'last.
  --   If abh is present as parameter, then res is the mixed residual.

  -- ON ENTRY :
  --   f        evaluable form of a (Laurent) polynomial system;
  --   abh      homotopy polynomials with absolute value coefficients;
  --   jf       Jacobian matrix of f;
  --   x        current approximate solution;
  --   verbose  is the verbose level.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   err      norm of the update vector;
  --   rco      estimate for the inverse condition number;
  --   res      residual, norm of the function value, if abh is absent,
  --            otherwise res is the mixed residual.

  procedure DoblDobl_LU_Newton_Step
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );
  procedure DoblDobl_LU_Newton_Step
              ( f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );
  procedure DoblDobl_LU_Newton_Step
              ( f,abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does one Newton step in double double complex arithmetic,
  --   using LU factorization to compute the update to x,
  --   and to estimate the inverse of the condition number rco.
  --   If abh is present as parameter, then res is the mixed residual.

  -- REQUIRED : f'last = x'last.

  -- ON ENTRY :
  --   f        evaluable form of a (Laurent) polynomial system;
  --   abh      homotopy polynomials with absolute value coefficients;
  --   jf       Jacobian matrix of f;
  --   x        current approximate solution;
  --   verbose  is the verbose level.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   err      norm of the update vector;
  --   rco      estimate for the inverse condition number;
  --   res      residual, norm of the function value.

  procedure DoblDobl_SVD_Newton_Step
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Jacobian_Circuits.Circuit;
                x : in out DoblDobl_Complex_Vectors.Vector;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );

  -- DESCRPTION :
  --   Does one Newton step in double double complex arithmetic,
  --   using a circuit to evaluate and differentiate,
  --   using the Singular Value Decomposition to compute the update to x
  --   and for the inverse condition number rco.
  --   Applies the method of Gauss-Newton, valid for f'last >= x'last.

  -- ON ENTRY :
  --   f        evaluable form of a polynomial system;
  --   jf       Jacobian matrix of f, defined as a circuit;
  --   x        current approximate solution;
  --   wrk      work space allocated for evaluated monomials;
  --   verbose  is the verbose level.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   err      norm of the update vector;
  --   rco      estimate for the inverse condition number;
  --   res      residual, norm of the function value.

  procedure DoblDobl_LU_Newton_Step
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Jacobian_Circuits.Circuit;
                x : in out DoblDobl_Complex_Vectors.Vector;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does one Newton step using a circuit to evaluate and differentiate,
  --   using LU factorization to compute the update to x,
  --   and to estimate the inverse of the condition number rco.

  -- REQUIRED : f'last = x'last.

  -- ON ENTRY :
  --   f        evaluable form of a polynomial system;
  --   jf       Jacobian matrix of f, defined as a circuit;
  --   x        current approximate solution;
  --   wrk      work space allocated for evaluated monomials;
  --   verbose  is the verbose level.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   err      norm of the update vector;
  --   rco      estimate for the inverse condition number;
  --   res      residual, norm of the function value.

-- WRAPPING ONE NEWTON STEP :

  procedure DoblDobl_Newton_Step
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );
  procedure DoblDobl_Newton_Step
              ( f,abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );
  procedure DoblDobl_Newton_Step
              ( f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does one Newton step in double double complex arithmetic,
  --   using LU for f'last = x'last and SVD for f'last > x'last.
  --   If the parameter abh is present, then res is the mixed residual.

  -- ON ENTRY :
  --   f        evaluable form of a (Laurent) polynomial system;
  --   abh      homotopy polynomials with absolute value coefficients,
  --   jf       Jacobian matrix of f;
  --   x        current approximate solution;
  --   verbose  is the verbose level.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   err      norm of the update vector;
  --   rco      estimate for the inverse condition number;
  --   res      residual, norm of the function value, if abh is absent,
  --            otherwise, equals the mixed residual.

  procedure DoblDobl_Newton_Step
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Jacobian_Circuits.Circuit;
                x : in out DoblDobl_Complex_Vectors.Vector;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does one Newton step in double double complex arithmetic,
  --   using a circuit to evaluate and differentiate, and
  --   using LU for f'last = x'last and SVD for f'last > x'last.

  -- ON ENTRY :
  --   f        evaluable form of a (Laurent) polynomial system;
  --   jf       Jacobian matrix of f;
  --   x        current approximate solution;
  --   wrk      work space allocated for evaluated monomials.

  -- ON RETURN :
  --   x        updated approximate solution;
  --   err      norm of the update vector;
  --   rco      estimate for the inverse condition number;
  --   res      residual, norm of the function value.

-- SEVERAL NEWTON STEPS :

  procedure Silent_Newton
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean );
  procedure Reporting_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean );
  procedure Silent_Newton
              ( f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean );
  procedure Reporting_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
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
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean );
  procedure Reporting_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Jacobian_Circuits.Circuit;
                x : in out DoblDobl_Complex_Solutions.Solution;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                epsxa,epsfa : in double_float; numit : in out natural32;
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

-- REFINING A LIST OF SOLUTIONS :

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

-- THE MAIN ROOT REFINERS :

  procedure Silent_Root_Refiner
               ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 verbose : in integer32 := 0 );
  procedure Silent_Root_Refiner
               ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s,refs : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 verbose : in integer32 := 0 );
  procedure Silent_Root_Refiner
               ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 verbose : in integer32 := 0 );
  procedure Silent_Root_Refiner
               ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 s,refs : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies Newton's method to refine roots of p in s.
  --   Stops when one conditions is satisfied:
  --   (1) numit >= max (reached maximum number of iterations),
  --   and for solutions x in s:
  --   (2) x.err < epsxa (update factor to x is less than epsxa),
  --   (3) x.res < epsfa (residual smaller than epsfa).
  --   The silent versions write no output.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   s        current approximate solutions;
  --   epsxa    accuracy requirement on update factor;
  --   epsfa    accuracy requirement on residual;
  --   tolsing  tolerance on inverse condition number of the Jacobian
  --            matrix at the root for consideration as singular solution;
  --   numit    number of iterations, must be zero on entry,
  --   max      maximum number of iterations allowed;
  --   verbose  is the verbose level.

  -- ON RETURN :
  --   s        updated approximate solutions;
  --   refs     list that does not include the path failures;
  --   numit    number of iterations spent on refining x;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean; verbose : in integer32 := 0 );
  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s,refs : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean; verbose : in integer32 := 0 );
  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean; verbose : in integer32 := 0 );
  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 s,refs : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies Newton's method to refine roots of p in s.
  --   Stops when one conditions is satisfied:
  --   (1) numit >= max (reached maximum number of iterations),
  --   and for solutions x in s:
  --   (2) x.err < epsxa (update factor to x is less than epsxa),
  --   (3) x.res < epsfa (residual smaller than epsfa).
  --   The reporting versions write output to file.

  -- ON ENTRY :
  --   file     for writing intermediate output and diagnostics;
  --   p        a polynomial system;
  --   s        current approximate solutions;
  --   epsxa    accuracy requirement on update factor;
  --   epsfa    accuracy requirement on residual;
  --   tolsing  tolerance on inverse condition number of the Jacobian
  --            matrix at the root for consideration as singular solution;
  --   numit    number of iterations, must be zero on entry,
  --   max      maximum number of iterations allowed;
  --   wout     if true, then information about each Newton update
  --            is written to file;
  --   verbose  is the verbose level.

  -- ON RETURN :
  --   s        updated approximate solutions;
  --   refs     list that does not include the path failures;
  --   numit    number of iterations spent on refining x.

-- REFINEMENT with DEFLATION : 

  procedure Silent_Root_Refiner
               ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; verbose : in integer32 := 0 );
  procedure Silent_Root_Refiner
               ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s,refs : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies Newton's method to refine roots of p in s,
  --   with deflation applied to singular solution if deflate is true.
  --   Stops when one conditions is satisfied:
  --   (1) numit >= max (reached maximum number of iterations),
  --   and for solutions x in s:
  --   (2) x.err < epsxa (update factor to x is less than epsxa),
  --   (3) x.res < epsfa (residual smaller than epsfa).
  --   The silent versions write no output.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   s        current approximate solutions;
  --   epsxa    accuracy requirement on update factor;
  --   epsfa    accuracy requirement on residual;
  --   tolsing  tolerance on inverse condition number of the Jacobian
  --            matrix at the root for consideration as singular solution;
  --   numit    number of iterations, must be zero on entry,
  --   max      maximum number of iterations allowed;
  --   deflate  to ask for deflation of the singular solutions.

  -- ON RETURN :
  --   s        updated approximate solutions;
  --   refs     list that does not include the path failures;
  --   numit    number of iterations spent on refining x;
  --   deflate  set to false if system is too large.

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; wout : in boolean;
                 verbose : in integer32 := 0 );
  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s,refs : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; wout : in boolean;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies Newton's method to refine roots of p in s,
  --   with deflation applied to singular solution if deflate is true.
  --   Stops when one conditions is satisfied:
  --   (1) numit >= max (reached maximum number of iterations),
  --   and for solutions x in s:
  --   (2) x.err < epsxa (update factor to x is less than epsxa),
  --   (3) x.res < epsfa (residual smaller than epsfa).
  --   The reporting versions write output to file.

  -- ON ENTRY :
  --   file     for writing intermediate output and diagnostics;
  --   p        a polynomial system;
  --   s        current approximate solutions;
  --   epsxa    accuracy requirement on update factor;
  --   epsfa    accuracy requirement on residual;
  --   tolsing  tolerance on inverse condition number of the Jacobian
  --            matrix at the root for consideration as singular solution;
  --   numit    number of iterations, must be zero on entry,
  --   max      maximum number of iterations allowed;
  --   deflate  to ask for deflation of the singular solutions;
  --   wout     if true, then information about each Newton update
  --            is written to file;
  --   verbose  is the verbose level.

  -- ON RETURN :
  --   s        updated approximate solutions;
  --   refs     list that does not include the path failures;
  --   numit    number of iterations spent on refining x;
  --   deflate  set to false if system is too large.

-- REFINEMENT with mixed residuals :

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 abh : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; wout : in boolean;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Refines the roots in sols for the system p, computing mixed residuals
  --   with the absolute coefficient homotopy polynomials in abh.

  -- ON ENTRY :
  --   file      file for writing diagnostics on;
  --   p         a polynomial system;
  --   abh       homotopy polynomials with absolute coefficients;
  --   sols      the start solutions;
  --   epsxa     maximum absolute error on the zero;
  --   epsfa     maximum absolute value for the residue;
  --   tolsing   tolerance on inverse condition number for singular solution;
  --   numit     the number of iterations, to be initialized with zero;
  --   max       maximum number of iterations per zero;
  --   deflate   if true, apply deflation to singular solutions;
  --   wout      has to be true when intermediate output is wanted;
  --   verbose   is the verbose level.

  -- ON RETURN :
  --   sols      a list of computed solutions;
  --   numit     the number of iterations;
  --   deflate   if set to false, then the system was too large to deflate.

end DoblDobl_Root_Refiners;
