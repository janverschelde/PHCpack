with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;

package Floating_Linear_Inequality_Solvers is

-- DESCRIPTION :
--   This package provides an incremental solving procedure for systems of
--   linear inequalities with floating point entries.

-- FORMAT OF INEQUALITIES :
--   A system of linear inequalities is represented by a matrix of
--   dimensions (1..n+1,1..m).  Each column contains one inequality,
--   for j in m'range(2), the jth inequality is given by
--     m(1,j)*x(1) + m(2,j)*x(2) + .. + m(m'last(1)-1,j)*x(x'last)
--                                                      >= m(m'last(1),j),
--   where x is any vector of range m'first(1)..m'last(1)-1.

-- OPERATIONS :
--   We distinguish three types :
--     selectors : extracting feasibility information;
--     modifiers : modifying the inequalities;
--     solvers   : solving means either to compute a solution
--                                   or to provide an inconsistency proof.

-- SELECTORS :

  function Evaluate ( m : Matrix; i : integer32; x : Vector )
                    return double_float;

  -- REQUIRED : x'first = m'first(1) and x'last = m'last(1)-1.

  -- DESCRIPTION :
  --   Returns the evaluation of the vector x at the ith inequality, i.e.:
  --     m(1,j)*x(1) + m(2,i)*x(2) + .. + m(m'last(1)-1,i)*x(x'last).

  function Satisfies ( m : Matrix; i : integer32; x : Vector;
                       tol : double_float ) return boolean;
  function Satisfies ( m : Matrix; x : Vector; tol : double_float )
                     return boolean;
  function Satisfies ( m : Matrix; fst,lst : integer32;
                       x : Vector; tol : double_float ) return boolean;

  -- REQUIRED : x'first = m'first(1) and x'last = m'last(1)-1.

  -- DESCRIPTION :
  --   Returns true if the given vector x is a feasible solution to the
  --   system of linear inequalities, defined in the columns of the matrix m.
  --   When the index i is provided, then only for the ith inequality the
  --   feasibility of the solution is checked.
  --   If fst and lst are provided then only the columns between fst and lst
  --   will be investigated.

  function First_Violated ( m : Matrix; x : Vector; tol : double_float )
                          return integer32;
  function First_Violated ( m : Matrix; fst,lst : integer32;
                            x : Vector; tol : double_float ) return integer32;

  -- DESCRIPTION :
  --   Returns the index of the first column that contains an inequality
  --   not satisfied by x.  If Satisfies(m,x), then m'last(2)+1 is returned.
  --   If fst and lst are provided then only the columns between fst and lst
  --   will be investigated, if Satisfies(m,fst,lst,x), then lst+1 is returned.

-- INCONSISTENCY CHECKS :

  function Inconsistent ( m : Matrix; i : integer32; tol : double_float )
                        return boolean;

  -- DESCRIPTION :
  --   Returns true if the ith inequality represents an inconsistent
  --   inequality like 0 >= c > 0, otherwise false is returned.

  function Inconsistent ( m : Matrix; cols : Standard_Integer_Vectors.Vector;
                          x : Vector; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Checks the inconsistency proof of the system of inequalities.
  --   The sum of x(i)*m(*,cols(i)), for i in x'range=cols'range, yields
  --   an inconsistent inequality, i.e.: the jth component of the sum
  --   is in absolute value smaller than tol and the last component is 
  --   a positive number larger than tol.

  function Inconsistent ( m : Matrix; i : integer32;
                          cols : Standard_Integer_Vectors.Vector; x : Vector;
                          tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   The inconsistency is here due to the addition of the ith inequality
  --   to the linear inequality system defined by m.
  --   The corresponding coefficient in the positive combination of the
  --   inequalities equals 1.0 and is not contained in the vector x.
  --   Note that this is the format of the proof as delivered by the solvers.

-- CONSTRUCTORS :

  procedure Scale ( m : in out Matrix; i : in integer32;
                    tol : in double_float );
  procedure Scale ( m : in out Matrix; tol : in double_float );

  -- DESCRIPTION :
  --   Scales the ith inequality or the whole system, by dividing to
  --   the square root of the inner product of the normal vector.

  procedure Center ( m : in out Matrix; x : in Vector );
  function  Center ( m : Matrix; x : Vector ) return Matrix;

  -- DESCRIPTION :
  --   Centers the inequalities in the system, so that the origin lies at x.

  procedure Intersect2D ( m : in Matrix; i,j : in integer32;
                          tol : in double_float;
                          x : out Vector; sing : out boolean );

  -- REQUIRED : x'length = 2 = m'last(1)-1.

  -- DESCRIPTION :
  --   Computes the intersection of the ith with the jth hyperplane.
  --   If sing, then both hyperplanes are (close to being) parallel,
  --   otherwise, the solution is equal to the vector x.

-- SOLVERS FOR THE INEQUALITY SYSTEMS :

  procedure Solve ( m : in Matrix; i : in integer32; tol : in double_float;
                    x : in out Vector; fail : out boolean;
                    cols : out Standard_Integer_Vectors.Vector );

  -- REQUIRED : i = First_Violated(m,x) <= m'last(2), x'range = cols'range.

  -- DESCRIPTION :
  --   Solves the inequality system, starting at the ith inequality.

  -- ON ENTRY :
  --   m         a matrix representation of an inequality system;
  --   i         first inequality that is not satisfied;
  --   tol       tolerance on the rounding errors;
  --   x         start solution.

  -- ON RETURN :
  --   x         solution to the first i inequalities in m,
  --             when not fail, then Satisfies(m(*,1..i),x);
  --             when fail, then x contains the coefficients for the
  --             inconsistency proof;
  --   fail      if false, then x satisfies the first i inequalities,
  --             otherwise, the system is inconsistent;
  --   cols      if fail, then cols indicates the columns of m to construct 
  --             the inconsistency proof, otherwise cols has no meaning.

  procedure Solve ( m : in Matrix; tol : in double_float; x : in out Vector;
                    k : out integer32;
                    cols : out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the inequality system by succesively applying the Solve above.
  --   When the parameter k > m'last(2), then the vector x is a solution,
  --   otherwise, k is the first inequality that makes the system inconsistent.

  generic

    with procedure Report ( x : in Vector; i : in integer32;
                            cont : out boolean );

    -- DESCRIPTION :
    --   x = current solution of the first i inequalities;
    --   when cont is set to false, the computations will be stopped,
    --   otherwise, the solver will continue.

  procedure Iterated_Solve ( m : in Matrix; tol : in double_float;
                             x : in out Vector; k : out integer32;
                             cols : out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the inequality system by a succesive application of the
  --   incremental solver.  After each new value of the solution, the procedure
  --   Report is called.  This enables to trace the flow of the computations
  --   and to stop then, whenever necessary.

end Floating_Linear_Inequality_Solvers;
