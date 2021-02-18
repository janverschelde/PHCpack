package Test_Solution_Filters is

-- DESCRIPTION :
--   Interactive tests on solution filters,
--   in double, double double, and quad double precision.

  procedure Standard_Filter;

  -- DESCRIPTION :
  --   Filters solutions read in standard double precision.

  procedure DoblDobl_Filter;

  -- DESCRIPTION :
  --   Filters solutions read in double double precision.

  procedure QuadDobl_Filter;

  -- DESCRIPTION :
  --   Filters solutions read in quad double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user for the precision
  --   and then calls the appropriate filter.

end Test_Solution_Filters;
