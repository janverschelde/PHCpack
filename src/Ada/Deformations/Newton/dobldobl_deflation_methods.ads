with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Standard_Natural_Vectors;
with DoblDobl_Complex_Vectors;
with Standard_Natural64_VecVecs;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Solutions;
with DoblDobl_Jacobian_Trees;
with DoblDobl_Evaluate_Deflation;

package DoblDobl_Deflation_Methods is

-- DESCRIPTION :
--   This package offers several implementations of Newton's method
--   with deflation for isolated singularities in double double arithmetic.
--   There are two main differences in the implementations:
--   (1) the easiest way is to treat the deflated systems symbolically,
--       just like any other polynomial system and apply Newton's method;
--   (2) a more efficient and algorithmic way is to exploit the block
--       structure of the Jacobian matrices in the deflation.

-- INTERACTIVE NEWTON with DEFLATION :

  procedure Display_Results_of_One_Newton_Step
              ( file : in file_type; z : in DoblDobl_Complex_Vectors.Vector;
                tol : double_float; err,rco,res : in double_double;
                rank,cnt : in natural32 );

  -- DESCRIPTION :
  --   Shows the results of one step with Newton's method.

  procedure One_Symbolic_Newton_Step
              ( ep : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                z : in out DoblDobl_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out double_double;
                s : out DoblDobl_Complex_Vectors.Vector;
                rank : out natural32 );
  procedure One_Symbolic_Newton_Step
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                z : in out DoblDobl_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out double_double;
                s : out DoblDobl_Complex_Vectors.Vector;
                rank : out natural32 );

  -- DESCRIPTION :
  --   Performs one step with Newton's method.

  -- ON ENTRY :
  --   file     for intermediate results and logging diagnostics;
  --   ep       evaluable form of a polynomial system;
  --   ejm      evaluable form of Jacobian matrix of the system;
  --   z        initial approximation for a solution of f;
  --   tol      tolerance to determine the numerical rank.
 
  -- ON RETURN :
  --   z        new approximation for a solution of f;
  --   err      magnitude of the correction added to z;
  --   rco      estimate for the inverse condition number of z;
  --   res      magnitude of the residual f(z);
  --   s        singular values of Jacobian matrix;
  --   rank     numerical rank of Jacobian matrix.

  procedure One_Algorithmic_Newton_Step
              ( ep : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in DoblDobl_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                z : out DoblDobl_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out double_double;
                s : out DoblDobl_Complex_Vectors.Vector;
                rank : out natural32 );
  procedure One_Algorithmic_Newton_Step
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in DoblDobl_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                z : out DoblDobl_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out double_double;
                s : out DoblDobl_Complex_Vectors.Vector;
                rank : out natural32 );

  -- DESCRIPTION :
  --   Performs one step with Newton's method.

  -- REQUIRED : nv'range = nq'range = 0..k,
  --   R1'range = B'range = h'range = 1..k, x'range = 0..k.

  -- ON ENTRY :
  --   file     for intermediate results and logging diagnostics;
  --   ep       evaluable form of a polynomial system;
  --   ejm      evaluable form of Jacobian matrix of the system;
  --   evt      remember table for deflation matrices;
  --   nd       remember table with Jacobian matrices;
  --   monkeys  keys used to hash the Jacobian matrices;
  --   k        current deflation stage;
  --   nv       nv(i) is number of columns in i-th deflation matrix;
  --   nq       nq(i) is number of rows in i-th deflation matrix;
  --   R1       R1(i) is the number of multipliers in i-th stage;
  --   B        B(i) is the random matrix used in the i-th stage;
  --   h        h(i) is the random vector to scale the i-th multipliers;
  --   x        x(0) contains values for the original variables,
  --            x(i) contains values for the i-th multipliers;
  --   tol      tolerance to determine the numerical rank.
 
  -- ON RETURN :
  --   x        new approximation for a solution of f;
  --   z        updated vector with new approximations,
  --            same information as in x, but in one long vector;
  --   err      magnitude of the correction added to z;
  --   rco      estimate for the inverse condition number of z;
  --   res      magnitude of the residual f(z);
  --   s        singular values of Jacobian matrix;
  --   rank     numerical rank of Jacobian matrix.

  procedure Interactive_Symbolic_Newton
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                z : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double;
                tol : in double_float; rank : out natural32 );

  -- DESCRIPTION :
  --   Calls Newton's method to find a better approximation of a root
  --   of p, starting at the vector z.

  procedure Interactive_Algorithmic_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in DoblDobl_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                tol : in double_float; rank : out natural32 );

  -- DESCRIPTION :
  --   Calls Newton's method at deflation stage k.

  -- REQUIRED : nv'range = nq'range = 0..k,
  --   R1'range = B'range = h'range = 1..k, x'range = 0..k.

  -- ON ENTRY :
  --   f        polynomial system deflated k times;
  --   evt      remember table for deflation matrices;
  --   nd       remember table with Jacobian matrices;
  --   monkeys  keys used to hash the Jacobian matrices;
  --   k        current deflation stage;
  --   nv       nv(i) is number of columns in i-th deflation matrix;
  --   nq       nq(i) is number of rows in i-th deflation matrix;
  --   R1       R1(i) is the number of multipliers in i-th stage;
  --   B        B(i) is the random matrix used in the i-th stage;
  --   h        h(i) is the random vector to scale the i-th multipliers;
  --   x        x(0) contains values for the original variables,
  --            x(i) contains values for the i-th multipliers;
  --   tol      tolerance to decide the rank of the Jacobian matrix.

  -- ON RETURN :
  --   x        refinement of the approximate solutions;
  --   err      magnitude of the last correction term;
  --   rco      estimate for the inverse of the condition number;
  --   res      residual of the solution.

  procedure Deflate
              ( p : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32 );

  -- DESCRIPTION :
  --   Replaces the system in p with the deflated system,
  --   added with m multipliers.

  procedure Add_Multipliers
              ( file : in file_type; output : in boolean;
                z : in out DoblDobl_Complex_Vectors.Link_to_Vector;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                m : in natural32; res : out double_double );

  -- DESCRIPTION :
  --   Adds values for the multipliers to the current root z.

  -- ON ENTRY :
  --   file     for intermediate output if output is true;
  --   output   flag to indicate whether to write multipliers and residual;
  --   z        current approximation for the root;
  --   f        polynomial system after deflation with m multipliers;
  --   m        number of multipliers to be added to z.

  -- ON RETURN :
  --   z        vector extended with m values for the multipliers;
  --   res      residual of the multiplier computation,
  --            only if < 0.1 will z be extended with the multipliers.

  procedure Add_Multipliers_for_Corank_One
              ( file : in file_type; output : in boolean;
                z : in out DoblDobl_Complex_Vectors.Link_to_Vector;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                n,m : in natural32; res : out double_double );

  -- DESCRIPTION :
  --   Adds values for the multipliers to the current root z,
  --   in the case the Jacobian matrix has corank one.

  -- ON ENTRY :
  --   file     for intermediate output if output is true;
  --   output   flag to indicate whether to write multipliers and residual; 
  --   z        current approximation for the root;
  --   f        polynomial system after deflation with m multipliers;
  --   n        number of new equations added in the deflation;
  --   m        number of multipliers to be added to z.

  -- ON RETURN :
  --   z        vector extended with m values for the multipliers;
  --   res      residual of the multiplier computation,
  --            only if < 0.1 will z be extended with the multipliers.

  procedure Add_Multipliers
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                res : out double_double );
  procedure Add_Multipliers
              ( file : in file_type; output : in boolean;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                res : out double_double );

  -- DESCRIPTION :
  --   Uses least squares to compute an initial value of the multipliers
  --   in the k-th deflation step.

  -- ON ENTRY :
  --   file     for intermediate results and logging diagnostics;
  --   output   flag to indicate whether to write multipliers and residual;
  --   f        evaluable form of a polynomial system;
  --   jm       evaluable form of Jacobian matrix of the system;
  --   nd       remember table with Jacobian matrices;
  --   monkeys  keys used to hash the Jacobian matrices;
  --   k        current deflation stage;
  --   nv       nv(i) is number of columns in i-th deflation matrix;
  --   nq       nq(i) is number of rows in i-th deflation matrix;
  --   R1       R1(i) is the number of multipliers in i-th stage;
  --   B        B(i) is the random matrix used in the i-th stage;
  --   h        h(i) is the random vector to scale the i-th multipliers;
  --   x        x(0) contains values for the original variables,
  --            x(i) contains values for the i-th multipliers.

  -- ON RETURN :
  --   x        x(k) has initial values for the multipliers;
  --   res      residual of the multiplier computation,
  --            only if < 0.1 will z be extended with the multipliers.

  procedure Add_Deflation_Data
              ( k : in integer32; m : in natural32;
                nv,nq,R1 : in out Standard_Natural_Vectors.Vector;
                B : in out DoblDobl_Complex_VecMats.VecMat;
                h : in out DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Updates the data structures for the next level of deflation,
  --   using m multipliers.

  -- REQUIRED : k > 0 and nv, nq, R1, B, and h have already been used
  --   in the previous deflation stages.

  -- ON ENTRY :
  --   k        next level of deflation, all data on entry is for i < k;
  --   m        number of multipliers to be used in next deflation level;
  --   nv       nv(i) is number of columns in i-th deflation matrix;
  --   nq       nq(i) is number of rows in i-th deflation matrix;
  --   R1       R1(i) is the number of multipliers in i-th stage;
  --   B        B(i) is the random matrix used in the i-th stage;
  --   h        h(i) is the random vector to scale the i-th multipliers;

  -- ON RETURN :
  --   The data nv, nq, R1, B, and h are updated for i = k.

  procedure Write_Results
              ( file : in file_type; i : natural32;
                p,dp : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                z : in DoblDobl_Complex_Vectors.Vector;
                err,rco,res : in double_double );

  -- DESCRIPTION :
  --   Writes the deflated system p and the root z to file.

  procedure Interactive_Symbolic_Deflation
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector;
                tol : in double_float );
  procedure Interactive_Symbolic_Deflation
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float );

  -- DESCRIPTION :
  --   Calls Newton's method on the list of solutions of p.
  --   If the list is empty, then the user is prompted for a solution.

  -- ON ENTRY :
  --   file     for intermediate results and diagnostics;
  --   p        polynomial system;
  --   sols     list of initial approximations for the solutions of p;
  --   size     size of the numbers if multi-precision arithmetic is used;
  --   tol      tolerance to determine the numerical rank of a matrix.

  procedure Interactive_Algorithmic_Deflation
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector;
                tol : in double_float );
  procedure Interactive_Algorithmic_Deflation
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                tol : in double_float );

  -- DESCRIPTION :
  --   This version of the deflation algorithm takes into account the
  --   block structure of the Jacobian matrices and is more efficient
  --   than the symbolic version.

-- NEWTON with DEFLATION and CLUSTERING :

  procedure Deflate_Solution
              ( file : in file_type; m : in integer32; output : in boolean;
                df : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Adds m values of multipliers to ls, using the deflated system in df.
  --   If output, then diagnostics and results are written to file.

  procedure Apply_Newton_Step
              ( file : in file_type; output : in boolean; step : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                tolrnk : in double_float; rank : out natural32 );

  -- DESCRIPTION :
  --   Does one Newton step to refine a solution in ls of the system f.

  -- ON ENTRY :
  --   file     for intermediate results and diagnostics if output on;
  --   output   flag to indicate if intermediate output is wanted;
  --   step     step counter, only used for output purposes;
  --   f        evaluable form of polynomial system;
  --   jf       evaluable form of Jacobian matrix of system;
  --   ls       initial approximation for a solution of f;
  --   tolrnk   tolerance to decide the numerical rank.

  -- ON RETURN :
  --   ls       new approximation for a solution of f;
  --   rank     numerical rank of Jacobian matrix at ls.

  procedure Apply_Newton_Step
              ( file : in file_type; output : in boolean; step : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                tolrnk : in double_float; rank : out natural32 );

  -- DESCRIPTION :
  --   Does one algorithmic Newton step on a k-deflated system.

  -- REQUIRED : k > 0 and nv'range = nq'range = 0..k,
  --   R1'range = B'range = h'range = 1..k, x'range = 0..k.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   output   flag to indicate if intermediate output is wanted;
  --   step     step counter, only used for output purposes;
  --   f        evaluable form of a polynomial system;
  --   jf       evaluable form of Jacobian matrix of the system;
  --   nd       remember table with Jacobian matrices;
  --   monkeys  keys used to hash the Jacobian matrices;
  --   k        current deflation stage;
  --   nv       nv(i) is number of columns in i-th deflation matrix;
  --   nq       nq(i) is number of rows in i-th deflation matrix;
  --   R1       R1(i) is the number of multipliers in i-th stage;
  --   B        B(i) is the random matrix used in the i-th stage;
  --   h        h(i) is the random vector to scale the i-th multipliers;
  --   x        x(0) contains values for the original variables,
  --            x(i) contains values for the i-th multipliers.
  --   tolrnk   tolerance to decide the numerical rank.

  -- ON RETURN :
  --   x        new approximation for the solution of f;
  --   err      size of the Newton correction term;
  --   rco      estimate for the inverse condition number;
  --   res      magnitude of the residual at the solution;
  --   rank     numerical rank of the Jacobian matrix at x.

  procedure Apply_Newton
              ( nit : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                tolrnk : in double_float; rank : out natural32 );
  procedure Apply_Newton
              ( file : in file_type; output : in boolean; nit : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out DoblDobl_Complex_Solutions.Link_to_Solution;
                tolrnk : in double_float; rank : out natural32 );

  -- DESCRIPTION :
  --   Does as many Newton steps as the value in nit.

  -- ON ENTRY :
  --   file     for intermediate results and diagnostics if output on;
  --   output   flag to indicate if intermediate output is wanted;
  --   nit      number of Newton steps that will be performed; 
  --   f        evaluable form of polynomial system;
  --   jf       evaluable form of Jacobian matrix of system;
  --   ls       initial approximation for a solution of f;
  --   tolrnk   tolerance to decide the numerical rank.

  -- ON RETURN :
  --   ls       new approximation for a solution of f;
  --   rank     numerical rank of Jacobian matrix at ls.

  procedure Apply_Newton
              ( nit : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                tolrnk : in double_float; rank : out natural32 );
  procedure Apply_Newton
              ( file : in file_type; output : in boolean; nit : in natural32;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in DoblDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in DoblDobl_Complex_VecMats.VecMat;
                h : in DoblDobl_Complex_VecVecs.VecVec;
                x : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double;
                tolrnk : in double_float; rank : out natural32 );

  -- DESCRIPTION :
  --   Does as many algorithmic Newton steps on a k-deflated system
  --   as the value of nit.

  -- REQUIRED : k > 0 and nv'range = nq'range = 0..k,
  --   R1'range = B'range = h'range = 1..k, x'range = 0..k.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   output   flag to indicate if intermediate output is wanted;
  --   nit      number of Newton steps that will be performed;
  --   f        evaluable form of a polynomial system;
  --   jf       evaluable form of Jacobian matrix of the system;
  --   nd       remember table with Jacobian matrices;
  --   monkeys  keys used to hash the Jacobian matrices;
  --   k        current deflation stage;
  --   nv       nv(i) is number of columns in i-th deflation matrix;
  --   nq       nq(i) is number of rows in i-th deflation matrix;
  --   R1       R1(i) is the number of multipliers in i-th stage;
  --   B        B(i) is the random matrix used in the i-th stage;
  --   h        h(i) is the random vector to scale the i-th multipliers;
  --   x        x(0) contains values for the original variables,
  --            x(i) contains values for the i-th multipliers.
  --   tolrnk   tolerance to decide the numerical rank.

  -- ON RETURN :
  --   x        new approximation for the solution of f;
  --   err      size of the Newton correction term;
  --   rco      estimate for the inverse condition number;
  --   res      magnitude of the residual at the solution;
  --   rank     numerical rank of the Jacobian matrix at x.
 
  procedure Symbolic_Deflation_and_Clustering
              ( file : in file_type; name : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                output : in boolean; maxitr,maxdef : in natural32;
                tolerr,tolres,tolrnk : in double_float );

  -- DESCRIPTION :
  --   Applies deflation to the given list of solutions and groups
  --   clusters of solutions to determine the multiplicity.
  --   The "Symbolic_" does it symbolically and is less efficient.

  -- ON ENTRY :
  --   file     must be opened for output for diagnostics;
  --   name     file name for the deflated systems and solutions;
  --   p        a polynomial system;
  --   sols     initial approximations for the solution of p;
  --   output   flag to indicate whether extra output will go to file;
  --   maxitr   maximal number of Newton iterations allowed on a root;
  --   maxdef   maximal number of deflations;
  --   tolerr   tolerance on the error sol.err to decide success/failure;
  --   tolres   tolerance on the residual to decide success or failure;
  --   tolrnk   tolerance to decide whether a matrix is singular.

  -- ON RETURN :
  --   As many files as the different types of deflation,
  --   all beginning with the given file name in name.
  --   The suffix of the file reflects the deflation path.

  procedure Algorithmic_Deflation_and_Clustering
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                maxitr,maxdef : in natural32;
                tolerr,tolres,tolrnk : in double_float );

  -- DESCRIPTION :
  --   Silent version of the other version with file and filename below.

  procedure Algorithmic_Deflation_and_Clustering
              ( file : in file_type; name : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                output : in boolean; maxitr,maxdef : in natural32;
                tolerr,tolres,tolrnk : in double_float );

  -- DESCRIPTION :
  --   Applies deflation to the given list of solutions and groups
  --   clusters of solutions to determine the multiplicity.
  --   The "Algorithmic_" is more efficient.

  -- ON ENTRY :
  --   file     must be opened for output for diagnostics;
  --   name     file name for the deflated systems and solutions;
  --   p        a polynomial system;
  --   sols     initial approximations for the solution of p;
  --   output   flag to indicate whether extra output will go to file;
  --   maxitr   maximal number of Newton iterations allowed on a root;
  --   maxdef   maximal number of deflations;
  --   tolerr   tolerance on the error sol.err to decide success/failure;
  --   tolres   tolerance on the residual to decide success or failure;
  --   tolrnk   tolerance to decide whether a matrix is singular.

  -- ON RETURN :
  --   sols     more accurate solutions, grouped with multiplicities.

end DoblDobl_Deflation_Methods;
