with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Test_Vector_Strings is

-- DESCRIPTION :
--   Writing vectors to strings and parsing strings from vectors.

  procedure Standard_Random_Test ( n : in integer32 );

  -- DESCRIPTION :
  --   Test on a random vector of dimension n in double precision.
  --   Writes the complex vector to a string
  --   and then parses the string for a vector.

  procedure DoblDobl_Random_Test ( n : in integer32 );

  -- DESCRIPTION :
  --   Test on a random vector of dimension n in double double precision.
  --   Writes the complex vector to a string
  --   and then parses the string for a vector.

  procedure QuadDobl_Random_Test ( n : in integer32 );

  -- DESCRIPTION :
  --   Test on a random vector of dimension n in quad double precision.
  --   Writes the complex vector to a string
  --   and then parses the string for a vector.

  procedure Multprec_Random_Test ( n : in integer32; size : in natural32 );

  -- DESCRIPTION :
  --   Test on a random vector of dimension n,
  --   with multiprecision numbers of the given size.
  --   Writes the vector to a string and then parses
  --   the string for the vector.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Vector_Strings;
