with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Test_Ordered_Evaluations is

-- DESCRIPTION :
--   Develops tests on power series of first and higher order terms, 
--   their evaluation in Laurent polynomials.

  procedure Test_First_Order ( dim,nbr,ord : in integer32 );

  -- DESCRIPTION :
  --   Generates a Laurent polynomial of nbr monomials in dim many variables
  --   and a power series of order ord.
  --   Then tests the first order evaluation.

  procedure Test_Second_Order ( dim,nbr,ord : in integer32 );

  -- DESCRIPTION :
  --   Generates a Laurent polynomial of nbr monomials in dim many variables
  --   and a power series of order ord.
  --   Then tests the second order evaluation.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and then launches a test.

end Test_Ordered_Evaluations;
