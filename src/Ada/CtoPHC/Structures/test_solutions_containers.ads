with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;

package Test_Solutions_Containers is

-- DESCRIPTION :
--   Defines test procedures on the solutions containers.

  procedure Standard_Test_Retrievals;

  -- DESRIPTION :
  --   Test the retrieval of the solutions in the container,
  --   in double precision.

  procedure DoblDobl_Test_Retrievals;

  -- DESRIPTION :
  --   Test the retrieval of the solutions in the container,
  --   in double double precision.

  procedure QuadDobl_Test_Retrievals;

  -- DESRIPTION :
  --   Test the retrieval of the solutions in the container,
  --   in quad double precision.

  procedure Multprec_Test_Retrievals;

  -- DESRIPTION :
  --   Test the retrieval of the solutions in the container,
  --   in multiprecision.

  procedure Standard_Test_Additions
              ( sols : in Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Tests the constructors of the double solutions container.

  procedure DoblDobl_Test_Additions
              ( sols : in DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Tests the constructors of the double double solutions container.

  procedure QuadDobl_Test_Additions
              ( sols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Tests the constructors of the quad double solutions container.

  procedure Multprec_Test_Additions
              ( sols : in Multprec_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Tests the constructors of the container for multiprecision solutions.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a selection from the testing menu
  --   and then runs the selected test.

end Test_Solutions_Containers;
