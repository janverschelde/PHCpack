with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Homotopy_Continuation_Parameters;   use Homotopy_Continuation_Parameters;

package Homotopy_Continuation_Parameters_io is

-- DESCRIPTION :
--   Provides procedures for output of the values of the parameter and
--   for interactive tuning of the homotopy continuation parameters.

  procedure put ( pars : in Parameters );
  procedure put ( file : in file_type; pars : in Parameters );

  -- DESCRIPTION :
  --   Returns default values for the parameters.

  procedure Prompt_for_Selection ( nbr : out natural32 );

  -- DESCRIPTION :
  --   Prompts the user for a number between 0 and 10,
  --   and allows for a retries if the given number is larger than 10.
  --   The number is returned in nbr.

  procedure Prompt_for_Parameter
              ( pars : in out Parameters; nbr : in natural32 );

  -- DESCRIPTION :
  --   Given a number in the range 1..10, prompts for a new value
  --   and assigns the corresponding value in pars.

  procedure Tune ( pars : in out Parameters );

  -- DESCRIPTION :
  --   Runs a loop to allow the user to interactive change 
  --   the values of the parameters.
  --   A good initial choice for the parameters are provided by
  --   Homotopy_Continuation_Parameters.Default_Values.

end Homotopy_Continuation_Parameters_io;
