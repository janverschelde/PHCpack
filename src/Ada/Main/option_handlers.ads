with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with String_Splitters;                   use String_Splitters;

package Option_Handlers is

-- DESCRIPTION :
--   Defines the handlers for the options of phc.
--   The option handlers define the main program,
--   in a top down design point of view of the software.

  procedure General_Help ( opt : in character );

  -- DESCRIPTION :
  --   Writes general help about an option to screen.

  procedure General_Help_Handler ( opts : in string );

  -- DESCRIPTION :
  --   Wraps the procedure General_Help.

  -- ON ENTRY :
  --   opts     options extracted from the command line arguments.

  procedure Help_Version_License
              ( args : in Array_of_Strings; name : in string );

  -- DESCRIPTION :
  --   Deals with the options --help, --version, and --license.
  --   This procedure is called when the first option is '-',
  --   the specifics of the handler are determined by the contents
  --   of the first arguments in args.
  --   The version number is written to the file with the given name
  --   if the name is not the empty string, otherwise,
  --   the version string is written to screen.

  procedure Full_Mode_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for running phc in full mode.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Good_Format_Handler
              ( opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -g.

  -- ON ENTRY :
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Find_and_Set_Seed
              ( args : in Array_of_Strings; opts : in string );

  -- DESCRIPTION :
  --   In case the '-0' option is on, the value of the seed
  --   will be extracted from the command line, and if nonzero
  --   it will be used to set the seed.
  --   If there is no value at the command line after '-0',
  --   then a fixed constant value is used as seed.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments.

  procedure EqnByEqn_Solver_Handler
              ( opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -a,
  --   the equation-by-equation solver.

  -- ON ENTRY :
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure BlackBox_Solver_Handler
              ( args : in Array_of_Strings; opts : in string;
                file1,file2,file3 : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -b,
  --   the blackbox solver for isolated solutions.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   file1    name of the input file with the target system;
  --   file2    name of the output file, if not combined with -b,
  --            otherwise, as -p -b, file2 contains the start system
  --            with the start solutions;
  --   file3    name of the output file, when running as -p -b.

  procedure Component_Solver_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -B,
  --   the blackbox numerical irreducible decomposition.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Scaling_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -s,
  --   for equation and coefficient scaling.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Reduction_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -d,
  --   for total degree reduction.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Root_Count_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -r,
  --   for root counts and start system construction.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Mixed_Volume_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -m,
  --   for mixed volumes and polyhedral homotopies.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Symbols_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -o.

  -- ON ENTRY :
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Continuation_Handler
              ( args : in Array_of_Strings; opts : in string;
                file1,file2,file3 : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -p,
  --   for continuation methods.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   file1    name of the input file with the target system;
  --   file2    name of the output file, if not combined with -b,
  --            otherwise, as -p -b, file2 contains the start system
  --            with the start solutions;
  --   file3    name of the output file, when running as -p -b.

  procedure Jumpstart_Handler
              ( args : in Array_of_Strings; opts : in string;
                file1,file2,file3 : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -q,
  --   for path tracking with jumpstarting.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   file1    name of the input file with the target system;
  --   file2    name of the input file with the start system;
  --   file3    name of the output file.

  procedure Algorithmic_Differentiation_Handler
              ( args : in Array_of_Strings; opts : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -j,
  --   for path trackers with algorithmic differentation.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments.

  procedure Enumeration_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -e,
  --   for the numerical Schubert calculus.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Feedback_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -k,
  --   to compute feedback laws with Pieri homotopies.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Decomposition_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -c,
  --   to run cascades of homotopies and diagonal homotopies.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Factorization_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -f,
  --   to filter and factor positive dimensional solution sets.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Series_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -u,
  --   to compute power series of solution curves.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Verification_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -v,
  --   to verify solutions.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Witness_Set_for_Hypersurface_Handler
              ( args : in Array_of_Strings; opts : in string;
                polyfile,logfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -l,
  --   to compute a witness set of a hypersurface.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the command line arguments;
  --   polyfile is the name of the input file;
  --   logfile  name of the output file.

  procedure Witness_Set_Intersection_Handler
              ( args : in Array_of_Strings;
                opts,witset_one,witset_two,logfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -w,
  --   to intersect two witness sets.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the comamnd line arguments;
  --   witset_one is the name of the first witness set;
  --   witset_two is the name of the second witness set;
  --   logfile  name of the output file.

  procedure Witness_Set_Sampler_Handler
              ( opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -y,
  --   to sample witness sets.

  -- ON ENTRY :
  --   opts     options extracted from the comamnd line arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Maple_Format_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -z,
  --   to interchange solutions from and into Maple formats.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

  procedure Python_Format_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string );

  -- DESCRIPTION :
  --   Defines the action for the option -x,
  --   to interchange solutions from and into Python dictionaries.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments;
  --   infile   name of the input file;
  --   outfile  name of the output file.

-- THE MAIN HANDLERS :

  procedure Handle ( args : in Array_of_Strings; 
                     opts,a1,a2,a3 : in string;
                     verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Handles the options in opts of the command line arguments in args,
  --   with arguments in a1, a2, and a3.

  -- ON ENTRY :
  --   args     command line arguments;
  --   opts     options extracted from the arguments,
  --            each option is represented by one character;
  --   a1       first argument that is not an option;
  --   a2       second argument that is not an option;
  --   a3       third argument that is not an option;
  --   verbose  the verbose level, for tracking procedure calls.

  -- REQUIRED :
  --   The options in opts are sorted, actions precede options.

  procedure Handle_no_Options ( infile,outfile : in string );

  -- DESCRIPTION :
  --   This handler is for when there are no command line arguments.

  -- ON ENTRY :
  --   infile   name of the input file;
  --   outfile  name of the output file.

end Option_Handlers;
