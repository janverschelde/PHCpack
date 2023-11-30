with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;

package Test_Solution_Strings is

-- DESCRIPTION :
--   Tests writing of solutions into strings and their parsing
--   into solutions of various precision levels.

  procedure Standard_Test_Write
              ( ls : in Standard_Complex_Solutions.Link_to_Solution );
  procedure DoblDobl_Test_Write
              ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution );
  procedure QuadDobl_Test_Write
              ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution );
  procedure Multprec_Test_Write
              ( ls : in Multprec_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Tests the writing and parsing of ls in double, double double,
  --   quad double precision, and arbitrary multiprecision.

  procedure Standard_Test_Write_Intro
              ( ls : in Standard_Complex_Solutions.Link_to_Solution );
  procedure DoblDobl_Test_Write_Intro
              ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution );
  procedure QuadDobl_Test_Write_Intro
              ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution );
  procedure Multprec_Test_Write_Intro
              ( ls : in Multprec_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Tests the writing and parsing of the introduction of ls,
  --   in double, double double, quad double precision,
  --   and in arbitrary multiprecision.

  procedure Standard_Test_Write_Vector
              ( ls : in Standard_Complex_Solutions.Link_to_Solution );
  procedure DoblDobl_Test_Write_Vector
              ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution );
  procedure QuadDobl_Test_Write_Vector
              ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution );
  procedure Multprec_Test_Write_Vector
              ( ls : in Multprec_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Tests the writing and parsing of the vector of ls,
  --   in double, double double, quad double precision,
  --   and in arbitrary multiprecision.

  procedure Standard_Test_Write_Diagnostics
              ( ls : in Standard_Complex_Solutions.Link_to_Solution );
  procedure DoblDobl_Test_Write_Diagnostics
              ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution );
  procedure QuadDobl_Test_Write_Diagnostics
              ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution );
  procedure Multprec_Test_Write_Diagnostics
              ( ls : in Multprec_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Tests the writing and parsing of the diagnostics of ls,
  --   in double, double double, quad double precision,
  --   and in arbitrary multiprecision.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts a menu and then runs the selected test.

end Test_Solution_Strings;
