with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;

package Test_Integer_Circuits is

-- DESCRIPTION :
--   Tests the computation of circuits of integer point configurations.

  function Random_Point_Configuration
             ( n,d : in integer32 ) return Matrix;

  -- DESCRIPTION :
  --   Generates a random point configuration of n points of dimension d,
  --   after prompting for bounds on the numbers.

  function Prompt_for_Columns ( d : in integer32 ) return Vector;

  -- DESCRIPTION :
  --   Prompts for a number less than or equal to d
  --   and then reads in as many indices as that given number.

  procedure Check_Circuit ( A : in Matrix );

  -- DESCRIPTION :
  --   Prompts for a selection of rows and then uses this
  --   selection as the basis for a circuit calculation.

  procedure Enumerate_Circuits ( A : in Matrix );

  -- DESCRIPTION :
  --   Prompts for a selection of rows and then uses this
  --   selection as the basis for a circuit enumeration.

  procedure Enumerate_Bases ( A : in Matrix; d : in integer32 );

  -- DESCRIPTION :
  --   Generates all selections of d columns of A which define a basis.

  procedure Enumerate_the_Circuits ( A : in Matrix; d : in integer32 );

  -- DESCRIPTION :
  --   Enumerates all circuits of the points defined by A.

  procedure Matrix_of_Circuits ( A : in Matrix; d : in integer32 );

  -- DESCRIPTION :
  --   Computes a matrix B of circuits and displays A*B.

  procedure Circuits ( n,d : in integer32 );

  -- DESCRIPTION :
  --   Generates a point configuration of n points of dimension d.
  --   Displays a menu and prompts for a test.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the number of points and the dimension,
  --   generates a random point configuration and displays
  --   the testing menu.

end Test_Integer_Circuits;
