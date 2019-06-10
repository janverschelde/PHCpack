with String_Splitters;                   use String_Splitters;

package Option_Handlers is

-- DESCRIPTION :
--   Defines the handlers for the options of phc.

  procedure General_Help ( opt : in character );

  -- DESCRIPTION :
  --   Writes general help about an option to screen.

  procedure Enumeration_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Invokes the numerical Schubert calculus.

  procedure Decomposition_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -c of phc.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extract from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Factorization_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -f of phc.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extract from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

end Option_Handlers;
