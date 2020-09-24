with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;

package Test_Standard_Laur_Homotopy is

-- DESCRIPTION :
--   Tests a Laurent polynomial homotopy in double precision.

  procedure Interactive_Test ( p,q : in Laur_Sys );

  -- DESCRIPTION :
  --   Prompts for complex numbers of the vectors where
  --   to evaluate the homotopy defined by target p and start q.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a target and start system
  --   and then runs an interactive test.

end Test_Standard_Laur_Homotopy;
