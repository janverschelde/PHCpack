package Test_Systems_Pools is

-- DESCRIPTION :
--   Tests some basic operations on the systems pools.

  procedure Standard_Test;

  -- DESCRIPTION :
  --   Runs a basic test to store systems of polynomials with
  --   complex coefficients in double precision.

  procedure DoblDobl_Test;

  -- DESCRIPTION :
  --   Runs a basic test to store systems of polynomials with
  --   complex coefficients in double double precision.

  procedure QuadDobl_Test;

  -- DESCRIPTION :
  --   Runs a basic test to store systems of polynomials with
  --   complex coefficients in quad double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user for the working precision
  --   and then launches the corresponding test.

end Test_Systems_Pools;
