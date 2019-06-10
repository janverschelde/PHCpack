with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with String_Splitters;                   use String_Splitters;

package Actions_and_Options is

-- DESCRIPTION :
--   Defines the actions and options of the phc executable.

  actions : constant string := "abBcdefghjklmopqrsuvwxyz";
  -- a : solve => equation-by-equation solver
  -- b : batch or black box processing
  -- B : black box numerical irreducible decomposition
  -- c : comp => numerical irreducible decomposition
  -- d : redu => reduction w.r.t. the total degree
  -- e : enum => numerical Schubert calculus
  -- f : fac  => factor pure dimensional solution set into irreducibles
  -- g : good => check if the input is a system in the valid format
  -- h : help => writes information about the general help
  -- j : adep => algorithmic differentiation path trackers
  -- k : feed => dynamic output feedback to control linear systems
  -- l : hyp  => witness set for hypersurface cutting with random line
  -- m : mvc  => mixed-volume computation
  -- o : symbOls -> get symbols using as variables in a system
  -- p : poco => polynomial continuation
  -- q : track => tracking solution paths
  -- r : roco => root counting methods
  -- s : scal => scaling of a polynomial system
  -- u : series => Newton's method for power series solutions
  -- v : vali => validation of solutions
  -- w : wit  => intersection of witness sets using diagonal homotopies
  -- x : dix  => Python dictionary output for solutions
  -- y : sam  => sampling points from an algebraic set
  -- z : zip  => strip output for solutions into Maple format

  options : constant string := "-0htV";
  -- - : displays version, general help, or how to cite
  -- 0 : zero seed for repeatable runs
  -- h : help => write information about a certain action
  -- t : task => use multitasking
  -- V : verbose level, -Vddd sets the level to ddd, default is zero

  function Position ( s : string; c : character ) return integer32;

  -- DESCRIPTION :
  --   If the the string s contains the character c, then its position in
  --   the string s will be returned.  Otherwise s'first-1 will be returned.
  --   This function is useful when string equals actions or options,
  --   to determine whether an option is valid or not.

  function Number_of_Tasks ( args : Array_of_Strings ) return natural32;

  -- DESCRIPTION :
  --   Given in args are the command line arguments.
  --   Returns the number of tasks of the argument -t on the command line.
  --   If there is no argument -t, then the number of cores is returned.

  function Verbose_Level ( args : Array_of_Strings ) return integer32;

  -- DESCRIPTION :
  --   Given in args are the command line arguments.
  --   Returns 0 if there is no option -V, or no number following the V,
  --   otherwise returns the verbose level, the number followed the V.

  function Find_Seed ( args : Array_of_Strings ) return natural32;

  -- DESCRIPTION :
  --   Given in args are the command line arguments.
  --   Reads the digits after the '-0' and returns the seed.
  --   If there is nothing after the '-0' then 0 is returned.

  function Scan_Precision
             ( args : Array_of_Strings; opt : character ) return natural32;

  -- DESCRIPTION :
  --   Given in args are the command line arguments.
  --   Returns the precision of the option defined by the charactor opt.
  --   1 : the -opt is followed by a space (or nothing);
  --   2 : double double precision, as we have -opt2 at the command line;
  --   4 : quad double precision is given as -opt4 at the command line.

  function Scan_Options ( args : Array_of_Strings ) return string;

  -- DESCRIPTION :
  --   Returns a string with the second character in args following '-'
  --   in the strings of args.

  function Sort_Options ( opts : string ) return string;

  -- DESCRIPTION :
  --   The string on return contains the options in opts in sorted order,
  --   sorted so that characters in the constant string actions defined above
  --   always precede the characters in the constant string options.

  function Get_Argument
             ( args : Array_of_Strings; k : integer32 ) return string;

  -- DESCRIPTION :
  --   Reads the kth argument from the command line.
  --   An argument is a string not proceeded by a `-' character.
  --   The empty string is returned when there is no argument.

end Actions_and_Options;
