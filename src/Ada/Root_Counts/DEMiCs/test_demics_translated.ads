with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Test_DEMiCs_Translated is

-- DESCRIPTION :
--   Runs tests on the mixed volume computation with the translated DEMiCs
--   on some benchmarks.

  procedure Test_Cyclic_Roots ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs tests on the cyclic n-roots benchmark systems.

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests and reports the results.

end Test_DEMiCs_Translated;
