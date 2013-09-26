with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_integer_Vectors;           use Standard_integer_Vectors;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;

package Integer_Linear_Inequalities is

-- DESCRIPTION :
--   This package provides some routines for implementing the Farkas lemma.

  procedure Integer_Complementary_Slackness
                 ( tableau : in out matrix; feasible : out boolean );
  procedure Integer_Complementary_Slackness
                 ( tableau : in out matrix; lastcol : in integer32; 
                   feasible : out boolean );

  -- DESCRIPTION :
  --   Solves the complementary slackness problem: determines
  --   whether there exists a positive combination of the columns
  --   such that the right-hand side is satisfied.

  -- ON ENTRY :
  --   tableau     inequalities as columns, last column is right hand side;
  --   lastcol     last column of interest in the tableau,
  --               when not provided, tableau'last(2)-1 is assumed.

  -- ON RETURN :
  --   tableau     modified tableau of inequalities;
  --   feasible    if true then the solution is feasible.

end Integer_Linear_Inequalities;
