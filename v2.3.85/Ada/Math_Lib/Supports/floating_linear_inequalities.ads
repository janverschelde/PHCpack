with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;

package Floating_Linear_Inequalities is

-- DESCRIPTION :
--   This procedure contains some routines for verifying the Farkas lemma.

  procedure Complementary_Slackness
                 ( tableau : in out Matrix; lastcol : in integer32;
                   rhs : in out Standard_Floating_Vectors.Vector; 
                   tol : in double_float;
                   solution : out Standard_Floating_Vectors.Vector;
                   columns : out Standard_Integer_Vectors.Vector;
                   feasible : out boolean );

  procedure Complementary_Slackness
                 ( tableau : in out Matrix;
                   rhs : in out Standard_Floating_Vectors.Vector;
                   tol : in double_float;
                   solution : out Standard_Floating_Vectors.Vector;
                   columns : out Standard_Integer_Vectors.Vector;
                   feasible : out boolean );

  -- DESCRIPTION :
  --   Solves the complementary slackness problem: determines
  --   whether there exists a positive combination of the columns
  --   such that the right hand side is satisfied.

  -- REQUIRED :
  --   rhs'range = solution'range = columns'range = tableau'range(1)

  -- ON ENTRY :
  --   tableau     inequalities as columns;
  --   lastcol     indicates the last significant column in the tableau,
  --                if not given, then lastcol = tableau'last(2);
  --   tol         tolerance to decide whether a number equals zero.
  --   rhs         right hand side vector;

  -- ON RETURN :
  --   tableau     modified tableau of inequalities;
  --   rhs         modified right hand side;
  --   solution    the computed solution;
  --   columns     indicates which columns has been used;
  --   feasible    if true then the solution is feasible.

end Floating_Linear_Inequalities;
