with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;

package Test_RanCff_Systems is

-- DESCRIPTION :
--   Test on the generation of polynomials and systems with given
--   supports and with random coefficients.

  procedure Write_to_File ( p : in Standard_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Gives the user the opportunity to write the system q to file,
  --   prompting for a file name.

  procedure Write_to_File ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Gives the user the opportunity to write the system q to file,
  --   prompting for a file name.

  procedure Write_to_File ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Gives the user the opportunity to write the system q to file,
  --   prompting for a file name.

  procedure Standard_Creation_of_Random_Polynomials;

  -- DESCRIPTION :
  --   Prompts for a list of exponents and then generates
  --   a system with random standard double complex coefficients
  --   with the given list of exponents.

  procedure DoblDobl_Creation_of_Random_Polynomials;

  -- DESCRIPTION :
  --   Prompts for a list of exponents and then generates
  --   a system with random double double complex coefficients
  --   with the given list of exponents.

  procedure QuadDobl_Creation_of_Random_Polynomials;

  -- DESCRIPTION :
  --   Prompts for a list of exponents and then generates
  --   a system with random quad double complex coefficients
  --   with the given list of exponents.

  procedure Read_Supports
              ( dim : in natural32;
                s : out Arrays_of_Integer_Vector_Lists.Array_of_Lists );

  -- DESCRIPTION :
  --   Prompts for as many lists of exponents as the range of s.

  procedure Standard_Random_Polynomial_System ( n,r : integer32 );

  -- DESCRIPTION :
  --   Generates n polynomials with r different supports
  --   with coefficients in double precision.

  procedure DoblDobl_Random_Polynomial_System ( n,r : integer32 );

  -- DESCRIPTION :
  --   Generates n polynomials with r different supports
  --   with coefficients in double double precision.

  procedure QuadDobl_Random_Polynomial_System ( n,r : integer32 );

  -- DESCRIPTION :
  --   Generates n polynomials with r different supports
  --   with coefficients in quad double precision.

  procedure Random_Polynomial_System ( precision : in character );

  -- DESCRIPTION :
  --   Prompts the user for the dimension and different supports
  --   and creates a random coefficient system with standard double,
  --   double double, or quad double coefficient, depending whether
  --   precision is '0', '1', or '2'.

  procedure Standard_Random_Coefficient_System;

  -- DESCRIPTION :
  --   Reads a polynomial system and constructs a random coefficient
  --   system in double precision.

  procedure DoblDobl_Random_Coefficient_System;

  -- DESCRIPTION :
  --   Reads a polynomial system and constructs a random coefficient
  --   system in double precision.

  procedure QuadDobl_Random_Coefficient_System;

  -- DESCRIPTION :
  --   Reads a polynomial system and constructs a random coefficient
  --   system in double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for choices.

end Test_RanCff_Systems;
