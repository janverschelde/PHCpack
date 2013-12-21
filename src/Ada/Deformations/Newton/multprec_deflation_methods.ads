with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural64_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_VecMats;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with Multprec_Complex_Solutions;
with Multprec_Jacobian_Trees;
with Multprec_Deflation_Trees;
with Multprec_Evaluate_Deflation;

package Multprec_Deflation_Methods is

-- DESCRIPTION :
--   This package offers several implementations of Newton's method
--   with deflation for isolated singularities in standard arithmetic.
--   There are two main differences in the implementations:
--   (1) the easiest way is to treat the deflated systems symbolically,
--       just like any other polynomial system and apply Newton's method;
--   (2) a more efficient and algorithmic way is to exploit the block
--       structure of the Jacobian matrices in the deflation.

-- INTERACTIVE NEWTON with DEFLATION :

  procedure Display_Results_of_One_Newton_Step
              ( file : in file_type;
                z : in Multprec_Complex_Vectors.Vector; tol : in double_float;
                err,rco,res : in Floating_Number; rank,cnt : in natural32 );

  -- DESCRIPTION :
  --   Shows the results of one step with Newton's method.

  procedure One_Symbolic_Newton_Step
              ( ep : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                z : in out Multprec_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out Floating_Number;
                s : out Multprec_Complex_Vectors.Vector;
                rank : out natural32 );
  procedure One_Symbolic_Newton_Step
              ( file : in file_type;
                ep : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                z : in out Multprec_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out Floating_Number;
                s : out Multprec_Complex_Vectors.Vector;
                rank : out natural32 );

  -- DESCRIPTION :
  --   Performs one step with Newton's method.

  -- ON ENTRY :
  --   file     for intermediate results and logging diagnostics;
  --   ep       evaluable form of a polynomial system;
  --   ejm      evaluable form of Jacobi matrix of the system;
  --   z        initial approximation for a solution of f;
  --   tol      tolerance to determine the numerical rank.
 
  -- ON RETURN :
  --   z        new approximation for a solution of f;
  --   err      magnitude of the correction added to z;
  --   rco      estimate for the inverse condition number of z;
  --   res      magnitude of the residual f(z);
  --   s        singular values of Jacob matrix;
  --   rank     numerical rank of the Jacobi matrix.

  procedure One_Algorithmic_Newton_Step
              ( ep : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in Multprec_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in Multprec_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in Multprec_Complex_VecMats.VecMat;
                h : in Multprec_Complex_VecVecs.VecVec;
                x : in out Multprec_Complex_VecVecs.VecVec;
                z : in out Multprec_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out Floating_Number;
                s : out Multprec_Complex_Vectors.Vector;
                rank : out natural32 );
  procedure One_Algorithmic_Newton_Step
              ( file : in file_type;
                ep : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                ejm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in Multprec_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in Multprec_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in Multprec_Complex_VecMats.VecMat;
                h : in Multprec_Complex_VecVecs.VecVec;
                x : in out Multprec_Complex_VecVecs.VecVec;
                z : in out Multprec_Complex_Vectors.Vector;
                tol : in double_float; err,rco,res : out Floating_Number;
                s : out Multprec_Complex_Vectors.Vector;
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
  --   z        same information as in x but in one long vector;
  --   tol      tolerance to determine the numerical rank.
 
  -- ON RETURN :
  --   x        new approximation for a solution of f;
  --   z        updated vector with new approximations;
  --   err      magnitude of the correction added to z;
  --   rco      estimate for the inverse condition number of z;
  --   res      magnitude of the residual f(z);
  --   s        singular values of Jacobian matrix;
  --   rank     numerical rank of Jacobian matrix.

  procedure Interactive_Symbolic_Newton
              ( file : in file_type;
                p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                z : in out Multprec_Complex_Vectors.Vector;
                err,rco,res : out Floating_Number;
                tol : in double_float; rank : out natural32 );

  -- DESCRIPTION :
  --   Calls Newton's method to find a better approximation of a root
  --   of p, starting at the vector z.

  procedure Interactive_Algorithmic_Newton
              ( file : in file_type;
                f : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                evt : in Multprec_Evaluate_Deflation.Link_to_Eval_Tree;
                nd : in Multprec_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in Multprec_Complex_VecMats.VecMat;
                h : in Multprec_Complex_VecVecs.VecVec;
                x : in out Multprec_Complex_VecVecs.VecVec;
                err,rco,res : out Floating_Number;
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
              ( p : in out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                m,size : in natural32 );

  -- DESCRIPTION :
  --   Replaces the system in p with the deflated system,
  --   added with m multipliers.

  procedure Add_Multipliers
              ( z : in out Multprec_Complex_Vectors.Link_to_Vector;
                f : in Multprec_Complex_Poly_Systems.Poly_Sys;
                m : in natural32 );

  -- DESCRIPTION :
  --   Adds values for the multipliers to the current root z.

  -- ON ENTRY :
  --   z        current approximation for the root;
  --   f        polynomial system after deflation with m multipliers;
  --   m        number of multipliers to be added to z.

  -- ON RETURN :
  --   z        vector extended with m values for the multipliers.

  procedure Add_Multipliers
              ( file : in file_type;
                f : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                nd : in Multprec_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R1 : in Standard_Natural_Vectors.Vector;
                B : in Multprec_Complex_VecMats.VecMat;
                h : in Multprec_Complex_VecVecs.VecVec;
                x : in out Multprec_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Uses least squares to compute an initial value of the multipliers
  --   in the k-th deflation step.

  -- ON ENTRY :
  --   file     for intermediate results and logging diagnostics;
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
  --   x        x(k) has initial values for the multipliers.

  procedure Write_Results
              ( file : in file_type; i : natural32;
                p,dp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                z : in Multprec_Complex_Vectors.Vector;
                err,rco,res : in Floating_Number );

  -- DESCRIPTION :
  --   Writes the deflated system p and the root z to file.

  procedure Interactive_Symbolic_Deflation
              ( file : in file_type;
                p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                size : in natural32; tol : in double_float );

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
                 p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 size : in natural32;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 tol : in double_float );

  -- DESCRIPTION :
  --   This version of the deflation algorithm takes into account the
  --   block structure of the Jacobian matrices and is more efficient
  --   than the symbolic version.

-- NEWTON with DEFLATION and CLUSTERING :

  function Equal ( n : natural32; tol : double_float;
                   s1,s2 : Multprec_Complex_Vectors.Vector )
                 return boolean;

  -- DESCRIPTION :
  --   Compares the first n components of the two vectors.
  --   Returns true if none of the first n components differs
  --   in absolute value more than the given tolerance.

  procedure Set_Multiplicity
              ( sols : in out Multprec_Complex_Solutions.Solution_List;
                s : in Multprec_Complex_Solutions.Solution;
                tol : in double_float; n,m : in natural32 );

  -- DESCRIPTION :
  --   Every solution in sols close to s within the given tolerance tol
  --   will be given the multiplicity m.

  function Number_of_Occurrences 
              ( sols : Multprec_Complex_Solutions.Solution_List;
                s : Multprec_Complex_Solutions.Solution;
                tol : in double_float; n : in natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of times the solution s occurs in the list.

  function Is_In ( sols : Multprec_Complex_Solutions.Solution_List;
                   v : Multprec_Complex_Vectors.Vector;
                   tol : double_float; n : natural32 ) return boolean;

  -- DESCRIPTION :
  --   Return true if the first n components of the solution vector v
  --   belong to the list sols within the given tolerance tol.

  procedure Remove_Duplicates
              ( sols : in out Multprec_Complex_Solutions.Solution_List;
                tol : in double_float; n : in natural32 );

  -- DESCRIPTION :
  --   Removes all duplicates from the list sols.

  procedure Compute_Multiplicities
              ( sols : in out Multprec_Complex_Solutions.Solution_List;
                tol : in double_float; n : in natural32 );

  -- DESCRIPTION :
  --   Sets the multiplicity for every solution in the list,
  --   grouping the solutions according to their clusters,
  --   using the tolerance tol as cluster radius.

  procedure Compute_Multiplicities
              ( nd : in out Multprec_Deflation_Trees.Node;
                tol : in double_float; n : in natural32 );

  -- DESCRIPTION :
  --   Computes the multiplicities of the solution lists in the tree.
  --   Duplicate entries in the solution lists are removed.

  procedure Deflate_Solution
              ( file : in file_type; m : in natural32; output : in boolean;
                df : in Multprec_Complex_Poly_Systems.Poly_Sys;
                ls : in out Multprec_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Adds m values of multipliers to ls, using the deflated system in df.
  --   If output, then diagnostics and results are written to file.

  procedure Apply_Newton
              ( file : in file_type; output : in boolean; cnt : in natural32;
                f : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out Multprec_Complex_Solutions.Link_to_Solution;
                tolrnk : in double_float; rank : out natural32 );

  -- DESCRIPTION :
  --   Applies Newton's method to refine the solution ls of the system f.

  -- ON ENTRY :
  --   file     for intermediate results and diagnostics if output on;
  --   output   flag to indicate if intermediate output is wanted;
  --   cnt      counter for number of Newton iterations;
  --   f        evaluable form of polynomial system;
  --   jf       evaluable form of Jacobi matrix of system;
  --   ls       initial approximation for a solution of f;
  --   tolrnk   tolerance to decide the numerical rank.

  -- ON RETURN :
  --   ls       new approximation for a solution of f;
  --   rank     numerical rank of Jacobi matrix at ls.

  procedure Symbolic_Deflation_and_Clustering
              ( file : in file_type; outfilename : in string;
                p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in Multprec_Complex_Solutions.Solution_List;
                output : in boolean; maxitr,maxdef,size : in natural32;
                tolerr,tolres,tolrnk : in double_float );
 
  -- DESCRIPTION :
  --   Applies deflation to the given list of solutions and computes
  --   the multiplicity of each solution by grouping clusters.

end Multprec_Deflation_Methods;
