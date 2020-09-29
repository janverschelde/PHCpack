package Newton_Fabry_on_Homotopy is

-- DESCRIPTION :
--   Computes the Newton-Fabry convergence radius
--   for artificial or natural-parameter homotopies.

  function Prompt_for_Precision return character;

  -- DESCRIPTION :
  --   Returns '1', '2', '3', '4', '5', '6', or '7'
  --   respectively corresponding to double, double double,
  --   triple double, quad double, penta double, octo double,
  --   or deca double precision as selected by the user.

  procedure Run_Newton_Fabry ( precision : in character );

  -- DESCRITPION :
  --   Given '1', 2', '3, '4', '5', '6', or '7' as precision,
  --   computes the convergence radius of a solution series
  --   in the corresponding precision.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the precision
  --   and then launches the computation in that precision.

end Newton_Fabry_on_Homotopy;
