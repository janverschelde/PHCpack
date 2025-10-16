with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Test_Weighted_Assignment is

-- DESCRIPTION :
--   Runs tests on solving the weighted assignment problem.

  procedure Run_Test_Example;

  -- DESCRIPTION :
  --   The example matrix on page 252 of the Papadimitriou-Steiglitz book
  --     7 2 1 9 4
  --     9 6 9 5 5
  --     3 8 3 1 8
  --     7 9 4 2 2
  --     8 4 7 4 8
  --   is used as input for the weighted assignment problem.

  procedure Run_Random_Example ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates a random matrix of dimension dim
  --   and then solves the weighted matching problem.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension, generates examples,
  --   and then runs tests.

end Test_Weighted_Assignment;
