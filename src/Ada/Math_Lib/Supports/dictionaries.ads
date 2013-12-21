with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Integer_Vectors;

package Dictionaries is

-- DESCRIPTION :
--   This package manages dictionaries, used to solve
--   the linear programming problem, in standard primal form
--
--      max <c,x>
--          <a,x> <= 0,   m linear inequalities
--   
--   where x = (1,x1,x2,..,xn);
-- 
--   and for the dual formulation of the problem
--
--      min <c,x>
--          <a,x> >= 0,   m linear inequalities
--
--   where x = (1,x1,x2,..,xn).

--   A dictionary is the linear system that corresponds to a feasible solution.
--   The right-hand side vector yields the primal solution, whereas the
--   vector of the objective function determines the dual solution.
--   A dictionary is reprented by a matrix of dimension (0..m,0..n) and
--   by two vectors which indicate the unknowns that are in/out the basis,
--   of respective dimensions (1..m) and (1..n). 

-- USAGE : be consequent when using primal and dual models.

-- IMPORTANT ASSUMPTIONS :
--   The initializers do not check on the validity on the input tableau.
--   In particular:
--    * primal : all right-hand sides to be positive;
--    * dual   : all coefficients of the objective function must be positive.
--   Otherwise, the claims on unboundedness and infeasibility may be false.

-- INITIALIZERS :

  procedure Init_Basis
                ( in_bas,out_bas : in out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION : Initializes the basis.

  -- ON ENTRY :
  --   in_bas     unknowns in the basis, vector with range 1..m,
  --   out_bas    unknowns out the basis, vector with range 1..n.

  function Init_Primal_Dictionary
                ( c : Standard_Floating_Vectors.Vector; a : Matrix )
                return Matrix;

  function Init_Dual_Dictionary
                ( c : Standard_Floating_Vectors.Vector; a : Matrix )
                return Matrix;

  -- DESCRIPTION : Initializes the matrix for the dictionary.

  -- MATRIX DIMENSIONS : c(0..m), a(1..m,0..n).

  -- ON ENTRY :
  --   c          coefficients of the target function;
  --   a          coefficients of the constraints, stored row-wise,
  --              with right-hand side at column zero.

  -- ON RETURN :
  --   Matrix of dimension (0..m,0..n), with first row contains -c, the
  --   rest is filled with a, or with -a in the dual case.

  procedure Primal_Init
                ( c : in Standard_Floating_Vectors.Vector;
                  a : in Matrix; dic : out Matrix;
                  in_bas,out_bas : in out Standard_Integer_Vectors.Vector );

  procedure Dual_Init
                ( c : in Standard_Floating_Vectors.Vector;
                  a : in Matrix; dic : out Matrix;
                  in_bas,out_bas : in out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   This procedure initializes the dictionary by consecutive application
  --   of the producedure Init_Basis and Init_Primal/Dual_Dictionary.

  -- ON ENTRY :
  --   c          the coefficients of the target function: c(0..n);
  --   a          the coefficients of the constraints: a(1..m,0..n).

  -- ON RETURN :
  --   dic        dictionary a matrix of ranges 0..m,0..n;
  --   in_bas     unknowns in the basis, vector with range 1..m;
  --   out_bas    unknowns out the basis, vector with range 1..n.

-- MODIFIERS :

  procedure Primal_Update
                ( dic : in out Matrix;
                  in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                  eps : in double_float; unbounded : out boolean );

  procedure Dual_Update 
                ( dic : in out Matrix;
                  in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                  eps : in double_float; feasible: out boolean );

  -- DESCRIPTION :
  --   This procedure performs one step of the simplex algorithm.

  -- ON ENTRY :
  --   dic        current matrix for the dictionary;
  --   in_bas     unknowns in the basis;
  --   out_bas    unknowns not in the basis;
  --   eps        used to determine wether a number is zero.

  -- ON RETURN :
  --   dic        modified matrix for the dictionary;
  --   in_bas     updated list of basis unknowns;
  --   out_bas    updated list of non-basis unknowns;
  --   unbounded  true when the solution is unbounded, false otherwise;
  --   feasible   false when infeasibility is detected, true otherwise.

-- SELECTORS :

  function Primal_Optimal ( dic : Matrix; eps : double_float ) return boolean;
  function Dual_Optimal   ( dic : Matrix; eps : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the dictionary offers an optimal solution.
  --   The parameter eps is used determine whether a number equals zero or not.

  function Optimum ( dic : Matrix ) return double_float;

  -- DESCRIPTION :
  --   Returns the current value of the objective function.

  function Primal_Solution
                  ( dic : Matrix;
                    in_bas,out_bas : Standard_Integer_Vectors.Vector )
                  return Standard_Floating_Vectors.Vector;

  function Dual_Solution
                  ( dic : Matrix;
                    in_bas,out_bas : Standard_Integer_Vectors.Vector )
                  return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the current values of the n unknowns in the primal or the
  --   dual formulation of the linear program.  If the program is dual, 
  --   then the signs of the dual solution should be reversed.

end Dictionaries;
