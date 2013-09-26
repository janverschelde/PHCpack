with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;

package Linear_Minimization is

-- DESCRIPTION :
--   Minimization of a linear cost function subject to linear inequalities.

-- THE MODEL :
--   minimize cost*x
--   subject to: cff(i)*x >= rhs(i)

-- GLOBAL MEANING OF THE DATA :
--   n            number of variables;
--   m            number of inequalities;
--   cost         vector of range 1..#variables
--   cff          matrix of dimension (1..#variables,1..#constraints)
--   rhs          right-hand side vector of range  1..#constraints
--   x            vector of unknowns, of range 1..#variables
--   active       indices to active constraints, range 1..#variables
--   passive      complement of active, range 1..#constraints-#variables
--   basis        square matrix, rows are coefficients of active constraints
--   binv         inverse of the basis, columns are directions

-- EVALUATORS :

  function Eval ( n,k : integer32; 
                  cff : Standard_Floating_Matrices.Matrix;
                  x : Standard_Floating_Vectors.Vector )
                return double_float;

  -- DESCRIPTION :
  --   Returns the dot product of cff(1..n,k) with x(1..n).

  function Eval ( n : integer32; 
                  cost,x : Standard_Floating_Vectors.Vector )
                return double_float;

  -- DESCRIPTION :
  --   Returns the dot product of cost(1..n) with x(1..n).

  procedure Eval ( n,m : integer32;
                   cff : in Standard_Floating_Matrices.Matrix;
                   rhs,cost,x : in Standard_Floating_Vectors.Vector;
                   res : out Standard_Floating_Vectors.Vector;
                   val : out double_float );

  -- DESCRIPTION :
  --   Evaluates the vector x in the optimization model.

  -- REQUIRED :
  --   cff'range(1), cost'range, and x'range contain 1..n,
  --   cff'range(2), rhs'range, and res'range contain 1..m.

  -- ON ENTRY :
  --   n          number of variables <= number of rows in cff;
  --   m          number of constraints <= number of columns in cff;
  --   cff        matrix with coefficients for the constraints;
  --   rhs        right-hand side vector;
  --   cost       cost vector, to be minimized;
  --   x          test vector.

  -- ON RETURN :
  --   res        residuals, cff*x - rhs, should be positive for a solution;
  --   val        cost*x, the value of the test vector.

-- SET UP ROUTINES :

  procedure Feasibility_Model
                ( n,m : in integer32;
                  cff : in Standard_Floating_Matrices.Matrix;
                  rhs : in Standard_Floating_Vectors.Vector;
                  cols : in out Standard_Integer_Vectors.Vector;
                  init : in out Standard_Floating_Vectors.Vector;
                  extcff : out Standard_Floating_Matrices.Matrix;
                  extrhs,cost : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns the model to decide on the feasibility of a system of
  --   linear inequalities.  We add one additional slack variable s,
  --   and minimize s, subject to cff*x + s >= rhs, s >= 0.
  --   If the optimal solution has s = 0, then the system is feasible.

  -- REQUIRED :
  --   extcff'range(1) = 1..n+1 = cost'range = init'range = cols'range,
  --   extcff'range(2) = 1..m+1 = extrhs'range,

  -- ON ENTRY :
  --   n          number of variables;
  --   m          number of constraints;
  --   cff        n*m-matrix with the coefficient of the constraints;
  --   rhs        right-hand side vector;
  --   cols       selected columns from the coefficient matrix;
  --   init       initial solution, vertex at the selected rows;
  --              init(n+1) is the lower bound for the slack variable.

  -- ON RETURN :
  --   cols       cols(n+1) contains the index newest inequality;
  --   init       init(n+1) is the value for the slack variable,
  --              if init(n+1) = 0, then the problem is already feasible.
  --   extcff     coefficient matrix extended with a column of ones;
  --   extrhs     right-hand side vector extended with zero;
  --   cost       the cost vector consists of zeros, except the last one.

  procedure Initialize 
                ( n,m : in integer32;
                  cff : in Standard_Floating_Matrices.Matrix;
                  rhs,x : in Standard_Floating_Vectors.Vector;
                  active : out Standard_Integer_Vectors.Vector;
                  extcff,basinv : out Standard_Floating_Matrices.Matrix;
                  extrhs : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   For a given feasible solution x, satifying cff*x >= rhs,
  --   the coefficient matrix and right-hand side vector are extended so
  --   that the indentity matrix is a basis.

  -- ON ENTRY :
  --   n          dimension of the problem, equals #variables;
  --   m          number of constraints;
  --   cff        n*m-matrix with coefficients of the constraints;
  --   cost       cost vector of range 1..n;
  --   rhs        right-hand side vector of range 1..m;
  --   x          initial feasible solution satisfies cff*x >= rhs.

  -- ON RETURN :
  --   active     active constrains, equals 1..n;
  --   extcff     cff, with n*n idendity matrix inserted as first n rows;
  --   extrhs     right-hand side vector augmented with x in first n rows;
  --   basinv     n*n identity matrix.

-- FEASIBILITY TEST :

  procedure Feasible
                ( n,m : in integer32; 
                  cff : in Standard_Floating_Matrices.Matrix;
                  rhs : in Standard_Floating_Vectors.Vector;
                  tol : in double_float;
                  sol : out Standard_Floating_Vectors.Vector;
                  fail : out boolean );

  -- DESCRIPTION :
  --   Decides on the feasibility of the system of linear inequalities.

  -- ON ENTRY :
  --   n          number of variables;
  --   m          number of inequalities;
  --   cff        n*m-matrix with the coefficients of the inequalities;
  --   rhs        vector of range 1..m with right-hand sides;
  --   tol        tolerance to decide whether a number is zero or not.

  -- ON RETURN :
  --   sol        a solution to the system, if the problem is feasible;
  --   fail       if true, then the system is not feasible.

  procedure Feasible
                ( n,m : in integer32; 
                  cff : in Standard_Floating_Matrices.Matrix;
                  rhs : in Standard_Floating_Vectors.Vector;
                  tol : in double_float;
                  sol : out Standard_Floating_Vectors.Vector;
                  binv : out Standard_Floating_Matrices.Matrix;
                  active : out Standard_Integer_Vectors.Vector;
                  fail : out boolean );

  -- DESCRIPTION :
  --   Decides on the feasibility of the system of linear inequalities.
  --   This version allows to proceed with a minimization in case the
  --   problem is deemed feasible.

  -- ON ENTRY :
  --   n          number of variables;
  --   m          number of inequalities;
  --   cff        n*m-matrix with the coefficients of the inequalities;
  --   rhs        vector of range 1..m with right-hand sides;
  --   tol        tolerance to decide whether a number is zero or not.

  -- ON RETURN :
  --   sol        a solution to the system, if the problem is feasible;
  --   binv       basis inverse to current solution, if feasible;
  --   active     current active constraints, if feasible;
  --   fail       if true, then the system is not feasible.

-- MINIMIZERS :

  procedure Minimize
                ( n,m : in integer32;
                  cff : in Standard_Floating_Matrices.Matrix;
                  cost,rhs : in Standard_Floating_Vectors.Vector;
                  x : in out Standard_Floating_Vectors.Vector;
                  tol : in double_float; infty : out boolean );

  -- DESCRIPTION :
  --   Minimizes a linear cost function subject to linear constraints.

  -- ON ENTRY :
  --   n          dimension of the problem, equals #variables;
  --   m          number of constraints;
  --   cff        n*m-matrix with coefficients of the constraints;
  --   cost       cost vector of range 1..n;
  --   rhs        right-hand side vector of range 1..m;
  --   x          initial feasible solution of range 1..n;
  --   tol        number with absolute value less than tol is considered zero.

  -- ON RETURN :
  --   x          optimized solution;
  --   infty      true if the optimization problem is unbounded.

  procedure Minimize
                ( n,m : in integer32;
                  cff : in Standard_Floating_Matrices.Matrix;
                  cost,rhs : in Standard_Floating_Vectors.Vector;
                  x : in out Standard_Floating_Vectors.Vector;
                  binv : in out Standard_Floating_Matrices.Matrix;
                  active : in out Standard_Integer_Vectors.Vector;
                  tol : in double_float; infty : out boolean );

  -- DESCRIPTION :
  --   Minimizes a linear cost function subject to linear constraints.
  --   This version takes as additional input the inverse of a basis
  --   with the corresponding indices to the active constraints.

  -- ON ENTRY :
  --   n          dimension of the problem, equals #variables;
  --   m          number of constraints;
  --   cff        n*m-matrix with coefficients of the constraints;
  --   cost       cost vector of range 1..n;
  --   rhs        right-hand side vector of range 1..m;
  --   x          initial feasible solution of range 1..n;
  --   binv       inverse of the basis, which is an n*n-matrix;
  --   active     indices to the active constraints, of range 1..n;
  --   tol        number with absolute value less than tol is considered zero.

  -- ON RETURN :
  --   x          optimized solution;
  --   binv       updated inverse of the basis;
  --   active     currently active constraints;
  --   infty      true if the optimization problem is unbounded.

  generic

    with procedure Report ( sol : in Standard_Floating_Vectors.Vector;
                            active : in Standard_Integer_Vectors.Vector;
                            continue : out boolean );

  procedure Reporting_Minimize
                ( n,m : in integer32;
                  cff : in Standard_Floating_Matrices.Matrix;
                  cost,rhs : in Standard_Floating_Vectors.Vector;
                  x : in out Standard_Floating_Vectors.Vector;
                  tol : in double_float; infty : out boolean );

  -- DESCRIPTION :
  --   After each update of the solution, the procedure Report is called.
  --   The iteration stops when continue is set to false.
  --   This procedure is useful for testing purposes.
  --   The variables have the same meaning as in the Minimize above.

-- ENUMERATOR :

  generic

    with procedure Report ( binv : in Standard_Floating_Matrices.Matrix;
                            active : in Standard_Integer_Vectors.Vector;
                            sol : in Standard_Floating_Vectors.Vector;
                            continue : out boolean );

    -- DESCRIPTION :
    --   This procedure is called each time a new basis is found.

    -- ON ENTRY :
    --   binv     newly updated inverse of the basis;
    --   active   columns of coefficient matrix used;
    --   sol      current solution.

    -- ON RETURN :
    --   continue if true, the iteration continues, otherwise it stops.

  procedure Enumerate_Feasible_Vertices
                ( n,m : in integer32;
                  cff : in Standard_Floating_Matrices.Matrix;
                  cost,rhs : in Standard_Floating_Vectors.Vector;
                  tol : in double_float; fail,infty : out boolean );

  -- DESCRIPTION :
  --   Enumerates all feasible vertices for the minimization problem.
  --   Each time a new vertex is found, Report is called.
  --   This procedure applies the reverse search principle. 

  -- REQUIRED : binv'range(2) = binv'range(1) = 1..n;
  --   mat'range(1) contains 1..n,
  --   mat'range(2) and mat_cols'range contain 1..m.

  -- ON ENTRY :
  --   cff        column-oriented coefficient matrix;
  --   cost       cost vector that needs to be minimized;
  --   rhs        right-hand side vector;
  --   tol        tolerance for the absolute value of nonzero numbers.

  -- ON RETURN :
  --   fail       true if the optimization problem is unfeasible;
  --   infty      true if the optimization problem is unbounded.

end Linear_Minimization;
