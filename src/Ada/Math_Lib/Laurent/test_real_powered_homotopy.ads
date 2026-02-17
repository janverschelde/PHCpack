with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Test_Real_Powered_Homotopy is

-- DESCRIPTION :
--   Tests the input/output procedures on a homotopy of Laurent polynomials,
--   where to every coefficient corresponds a real powered series.

  procedure Test_Random_Polynomial ( nbr,nvr,size : in integer32 );

  -- DESCRIPTION :
  --   Tests output/input of a random polynomial of nbr monomials 
  --   in nvr variables, where each series has the given size.

  procedure Test_String_Polynomial ( nbr,nvr,size : in integer32 );

  -- DESCRIPTION :
  --   Tests the string representations of a Laurent polynomial
  --   with nbr terms in nvr variables, with as coefficients
  --   real powered series of the given size.

  procedure Test_Random_System ( dim,nbr,size : in integer32 );

  -- DESCRIPTION :
  --   Generates a random homotopy system of dimension dim,
  --   with nbr monomials in each equation and where the coefficients
  --   are real powered series of the given size.
  --   Tests the output and the input.

  procedure Test_Regularity ( dim,nbr,size : in integer32 );

  -- DESCRIPTION :
  --   Generates a random homotopy system of dimension dim,
  --   with nbr monomials in each equation and where the coefficients
  --   are real powered series of the given size.
  --   Tests the regularity of the constant coefficient matrix.

  procedure Main;

  -- DESCRIPTION :
  --   Display a menu with all tests
  --   and then runs the selected test.

end Test_Real_Powered_Homotopy;
