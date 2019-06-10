with String_Splitters;                   use String_Splitters;

package Option_Handlers is

-- DESCRIPTION :
--   Defines the handlers for the options of phc.

  procedure General_Help ( opt : in character );

  -- DESCRIPTION :
  --   Writes general help about an option to screen.

  procedure Good_Format_Handler
              ( opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -g.

  -- ON ENTRY :
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure EqnByEqn_Solver_Handler
              ( opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -a.

  -- ON ENTRY :
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Component_Solver_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -B.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Scaling_Handler
              ( opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -s.

  -- ON ENTRY :
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Reduction_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -d.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Root_Count_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -r.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Mixed_Volume_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -m.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Continuation_Handler
              ( args : in Array_of_Strings; opts : in string;
                file1,file2,file3 : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -p.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   file1    name of the input file with the target system;
  --   file2    name of the output file, if not combined with -b,
  --            otherwise, as -p -b, file2 contains the start system
  --            with the start solutions;
  --   file3    name of the output file, when running as -p -b.

  procedure Jumpstart_Handler
              ( opts : in string; file1,file2,file3 : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -q.

  -- ON ENTRY :
  --   opts     options extracted from the arguments;
  --   file1    name of the input file with the target system;
  --   file2    name of the input file with the start system;
  --   file3    name of the output file.

  procedure Algorithmic_Differentiation_Handler
              ( opts : in string; file1,file2,file3 : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -j.

  -- ON ENTRY :
  --   opts     options extracted from the arguments;
  --   file1    name of the input file with the target system;
  --   file2    name of the input file with the start system;
  --   file3    name of the output file.

  procedure Enumeration_Handler
              ( opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -e.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Decomposition_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -c.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Factorization_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -f.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Series_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -u.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Verification_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -v.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Witness_Set_for_Hypersurface_Handler
              ( args : in Array_of_Strings; opts : in string;
                polyfile,logfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -l.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   polyfile is the name of the input file;
  --   logfile  name of the output file.

  procedure Witness_Set_Intersection_Handler
              ( opts,witset_one,witset_two,logfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -w.

  -- ON ENTRY :
  --   opts     options extracted from the arguments;
  --   witset_one is the name of the first witness set;
  --   witset_two is the name of the second witness set;
  --   logfile  name of the output file.

  procedure Witness_Set_Sampler_Handler
              ( opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -y.

  -- ON ENTRY :
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Maple_Format_Handler
              ( opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -z.

  -- ON ENTRY :
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Python_Format_Handler
              ( opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -x.

  -- ON ENTRY :
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

-- THE MAIN HANDLER :

  procedure Handle ( args : in Array_of_Strings; 
                     opts,a1,a2,a3 : in string );

  -- DESCRIPTION :
  --   Handles the options in opts of the command line arguments in args,
  --   with arguments in a1, a2, and a3.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments,
  --            each option is represented by one character;
  --   a1       first argument that is not an option;
  --   a2       second argument that is not an option;
  --   a3       third argument that is not an option.

  -- REQUIRED :
  --   The options in opts are sorted, actions precede options.

end Option_Handlers;
