with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Test_DEMiCs_Translated is

-- DESCRIPTION :
--   Runs tests on the mixed volume computation with the translated DEMiCs
--   on some benchmarks.

  procedure Test_Labels ( dim : in integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tests the computation of the labels to the mixed cells of
  --   a cyclic n-roots problem where n = dim.

  procedure Test_Cyclic_Roots ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs tests on the cyclic n-roots benchmark systems.

  procedure Test_Reformulated_Cyclic ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   The reformulated cyclic n-root systems is a Laurent system.
  --   This test is thus on using Laurent systems as input.

  procedure Test_Eigenvalue_Problem
              ( dim : in integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Generates a sparse eigenvalue problem of dimension dim.

  procedure Test_User_Lifting ( vrblvl : in integer32 := 0 );
 
  -- DESCRIPTION :
  --   Tests user defined lifting on a small example.

  procedure Test_Stable_Mixed_Volume ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tests the stable mixed volume computation.

  procedure Test_Call_DEMiCs ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tests the interface code to call DEMiCs.

  procedure Interactive_Test_Call_DEMICS ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for a polynomial system and then calls DEMiCs.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests and reports the results.

end Test_DEMiCs_Translated;
