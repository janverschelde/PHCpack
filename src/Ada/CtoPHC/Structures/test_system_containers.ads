with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;

package Test_System_Containers is

-- DESCRIPTION :
--   Procedures in this package enable to test the operations
--   of the system containers.

  procedure Standard_Test_Retrievals;

  -- DESCRIPTION :
  --   Tests the selectors in the container for polynomial systems
  --   with coefficients in double precision.

  procedure DoblDobl_Test_Retrievals;

  -- DESCRIPTION :
  --   Tests the selectors in the container for polynomial systems
  --   with coefficients in double double precision.

  procedure QuadDobl_Test_Retrievals;

  -- DESCRIPTION :
  --   Tests the selectors in the container for polynomial systems
  --   with coefficients in quad double precision.

  procedure Multprec_Test_Retrievals;

  -- DESCRIPTION :
  --   Tests the selectors in the container for polynomial systems
  --   with coefficients in multiprecision.

  procedure Standard_Test_Laurent_Retrievals;

  -- DESCRIPTION :
  --   Test on retrieving data from the Laurent system container,
  --   for coefficients in double precision.

  procedure DoblDobl_Test_Laurent_Retrievals;

  -- DESCRIPTION :
  --   Test on retrieving data from the Laurent system container,
  --   for coefficients in double double precision.

  procedure QuadDobl_Test_Laurent_Retrievals;

  -- DESCRIPTION :
  --   Test on retrieving data from the Laurent system container,
  --   for coefficients in double double precision.

  procedure Multprec_Test_Laurent_Retrievals;

  -- DESCRIPTION :
  --   Test on retrieving data from the Laurent system container,
  --   for multiprecision complex coefficients.

  procedure Standard_Test_Additions
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Tests the constructors of the double polynomial container.

  procedure DoblDobl_Test_Additions
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Tests the constructors of the double double polynomial container.

  procedure QuadDobl_Test_Additions
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Tests the constructors of the quad double polynomial container.

  procedure Multprec_Test_Additions
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Tests the constructors of the multiprecision polynomial container.

  procedure Standard_Test_Laurent_Additions
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Tests the constructors of the double Laurent container.

  procedure DoblDobl_Test_Laurent_Additions
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Tests the constructors of the double double Laurent container.

  procedure QuadDobl_Test_Laurent_Additions
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Tests the constructors of the quad double Laurent container.

  procedure Multprec_Test_Laurent_Additions
              ( p : in Multprec_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Tests the constructors of the multiprecision Laurent container.

  procedure Main;

  -- DESCRIPTION :
  --   Displays the testing menu, prompts for a selection,
  --   and then runs the selected test.

end Test_System_Containers;
