with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

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

  procedure Run_Newton_Fabry
              ( nbtasks : in natural32; precision : in character;
                vrblvl : in integer32 := 0 );

  -- DESCRITPION :
  --   Given the number of tasks in nbtasks and
  --   given '1', 2', '3, '4', '5', '6', or '7' as precision,
  --   computes the convergence radius of a solution series
  --   in the corresponding precision.
  --   The value of the verbose level is in vrblvl.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the precision
  --   and then launches the computation in that precision.

  procedure Test;

  -- DESCRIPTION :
  --   Interactive procedure to generate a test case.

end Newton_Fabry_on_Homotopy;
