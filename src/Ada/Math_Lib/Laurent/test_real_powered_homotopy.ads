with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Test_Real_Powered_Homotopy is

-- DESCRIPTION :
--   Tests the input/output procedures on a homotopy of Laurent polynomials,
--   where to every coefficient corresponds a real powered series.

  procedure Test_Random_Polynomial ( nbr,nvr,size : in integer32 );

  -- DESCRIPTION :
  --   Tests a random polynomial of nbr monomials in nvr variables,
  --   where each series has the given size.

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

  procedure Main;

  -- DESCRIPTION :
  --   Runs all tests.

end Test_Real_Powered_Homotopy;
