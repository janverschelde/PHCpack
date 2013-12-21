with text_io,integer_io;                 use text_io,integer_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;     use Standard_Complex_Laur_Jacomats;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Exponent_Vectors;                   use Exponent_Vectors;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Polyhedral_Coefficient_Tracking is

-- DESCRIPTION :
--   In a polyhedral coefficient homotopy, the coefficients are random
--   complex numbers multiplied with some (possibly) high powers of the
--   continuation parameter.  To deal with the numerical instabilities
--   caused by these high powers, an exponentional transformation of
--   on the continuation parameter is applied in the path trackers
--   provided by this package.

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
                tol : in double_float; max : in natural;
                x,y : in out Standard_Complex_Vectors.Vector;
                nit : out natural; fail : out boolean );

  procedure Reporting_Apply_Newton
              ( file : in file_type; hq : in Eval_Coeff_Laur_Sys;
                ctm : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                tol : in double_float; max : in natural;
                x,y : in out Standard_Complex_Vectors.Vector;
                nit : out natural; fail : out boolean );

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

  procedure Silent_Track_One_Path
              ( hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                pow : in Standard_Floating_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                sol : in out Solution;
                nsteps : out natural; fail : out boolean  );

  procedure Reporting_Track_One_Path
              ( file : in file_type;
                hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                pow : in Standard_Floating_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                sol : in out Solution;
                nsteps : out natural; fail : out boolean  );

  -- DESCRIPTION :
  --   Tracks one path, starting at one solution "sol" of start system
  --   corresponding to a mixed cell with inner normal given in "nor".

  -- ON ENTRY :
  --   file     to write all kinds of diagnostics of the path;
  --   hq       evaluable form of coefficient Laurent system;
  --   ctm      work space to hold coefficient for fixed value
  --            of the continuation parameter;
  --   pow      powers of the continuation parameter in the homotopy
  --            defined by one mixed cell;
  --   coeffv   coefficients in the polyhedral coefficient homotopy;
  --   jacmat   coefficient Jacobian matrix;
  --   mulfac   multiplication factors for Jacobian matrix;
  --   sol      one start solution.

  -- ON RETURN :
  --   sol      solution at the end of the path, if not fail;
  --   nsteps   number of steps done in the path tracking;
  --   fail     true if path failed to converged in the allowed #steps.

  procedure Track_Paths_for_Cell
              ( file : in file_type; lq : in Laur_Sys;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                ctm : in out Standard_Complex_VecVecs.VecVec;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector;
                ind : in natural; mic : in Mixed_Cell;
                qsols,last : in out Solution_List );

  -- DESCRIPTION :
  --   Tracks the paths using the polyhedral homotopies defined by one
  --   mixed cell, appending to the list of solutions.

  -- ON ENTRY :
  --   file     for writing diagnostics;
  --   lq       random coefficient start system;
  --   ctm      work space for the coefficients of the system;

  -- ON RETURN :
  --   ctm      modified work space;
  --   qsols    updated solution list;
  --   last     pointer to the last element in qsols.

  procedure Track_Paths_for_Subdivision
              ( file : in file_type; lq : in Laur_Sys;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector;
                sub : in Mixed_Subdivision; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Tracks the paths using the polyhedral homotopies defined by the
  --   mixed cells in the subdivision.

  -- ON ENTRY :
  --   file     for intermediate output;
  --   lq       a random coefficient Laurent polynomial system;
  --   ls       lifted supports of the polynomial system lq;

  -- ON RETURN :
  --   qsols    solutions to the system lq.

  procedure Track_Paths_for_Subdivision
              ( infile,outfile : in file_type; m : in natural;
                lq : in Laur_Sys;
                vs : in Standard_Floating_VecVecs.Array_of_VecVecs;
                ls : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                hq : in Eval_Coeff_Laur_Sys;
                coeffv : in Standard_Complex_VecVecs.VecVec;
                expvec : in Exponent_Vectors_Array;
                jacmat : in Eval_Coeff_Jaco_Mat; mulfac : in Mult_Factors;
                mix : in Standard_Integer_Vectors.Vector;
                qsols : out Solution_List );

  -- DESCRIPTION :
  --   Tracks the paths using the polyhedral homotopies defined by the
  --   mixed cells on the input file.

  -- REQUIRED :
  --   The current position at the input file must be at the beginning
  --   of a mixed cell.
  
  -- ON ENTRY :
  --   infile   input file with a regular mixed-cell configuration,
  --            properly positioned at the beginning of a mixed cell;
  --   outfile  output file for intermediate results;
  --   m        total number of mixed cells on the input file;
  --   lq       a random coefficient Laurent polynomial system;
  --   vs       lifted supports of the polynomial system lq;
  --   ls       lifted supports of the polynomial system lq;
  --   hq       polyhedral coefficient homotopy;

  -- ON RETURN :
  --   qsols    solutions to the system lq.

  procedure Polyhedral_Continuation
              ( file : in file_type; n : in natural; p : in Poly_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                sub : in Mixed_Subdivision;
                q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Creates a random coefficient system with the vertices used in the
  --   mixed-cell configuration, constructs a polyhedral homotopy, and
  --   then calls the path trackers to solve the random coefficient system.
 
  -- ON ENTRY :
  --   file     output file for intermediate results and diagnostics;
  --   n        ambient dimension before the lifting;
  --   p        polynomial system for which a random coefficient
  --            start system will be constructed;
  --   mix      type of mixture read from infile counts the number
  --            of occurrences of each support in the subdivision.

  -- ON RETURN :
  --   q        a random coefficient start system to solve p;
  --   qsols    solutions of q.

  procedure Polyhedral_Continuation
              ( infile,outfile : in file_type; 
                n : in natural; p : in Poly_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Creates a random coefficient system with the vertices used in the
  --   mixed-cell configuration, constructs a polyhedral homotopy, and
  --   then calls the path trackers to solve the random coefficient system.
  --   This jumpstarting version does not require for all mixed cells
  --   to be in the internal memory at the same time.

  -- REQUIRED :
  --   The regular mixed-cell configuration on the input file is in
  --   its labeled presentation, with the lifted sets in front.

  -- ON ENTRY :
  --   infile   input file with a regular mixed-cell configuration,
  --            positioned at the start of the lifted support sets;
  --   outfile  output file for intermediate results and diagnostics;
  --   n        ambient dimension before the lifting;
  --   p        polynomial system for which a random coefficient
  --            start system will be constructed;
  --   mix      type of mixture read from infile counts the number
  --            of occurrences of each support in the subdivision.

  -- ON RETURN :
  --   q        a random coefficient start system to solve p;
  --   qsols    solutions of q.
  
end Polyhedral_Coefficient_Tracking;
