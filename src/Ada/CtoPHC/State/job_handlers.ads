with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;

package Job_Handlers is

-- DESCRIPTION :
--   Defines the functions to handle the jobs for the C gateway.
--   The functions return 0 if all went well, or else the job code.
--   If the verbose level vrbvlv is positive, then the name of the
--   function is written to screen, to track bugs faster.

  function Version_String ( a : C_intarrs.Pointer;
                            b : C_intarrs.Pointer;
                            vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the PHCpack version string in b,
  --   with in a the number of characters in the string.

  function Get_Seed ( a : C_intarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Gets the seed used in the random number generators
  --   and returns its value in the parameter a on return.
  --   Having the seed could be useful for reproducible runs.

  function Set_Seed ( a : C_intarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Takes the number stored in a
  --   and uses its value to initialize the seed
  --   for the random number generator.

  function Standard_Polynomial_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32;

  -- DESCRIPTION :
  --   Solves the polynomial system in the system container
  --   with the blackbox solver in standard double precision
  --   and puts the solutions into the solution container.
  --   The two parameters on entry are as follows:
  --     a[0] : 1 or 2 for to be silent or not,
  --     a[1] : the number of tasks.
  --   On return in a[0] is the root count and if not silent,
  --   then a[1] contains the number of characters in b,
  --   where b is the root counter output string.

  function Standard_Laurent_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32;

  -- DESCRIPTION :
  --   Solves the Laurent system in the system container
  --   with the blackbox solver in standard double precision
  --   and puts the solutions into the solution container.
  --   The two parameters on entry are as follows:
  --     a[0] : 1 or 2 for to be silent or not,
  --     a[1] : the number of tasks.
  --   On return in a[0] is the root count and if not silent,
  --   then a[1] contains the number of characters in b,
  --   where b is the root counter output string;

  function DoblDobl_Polynomial_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32;

  -- DESCRIPTION :
  --   Solves the polynomial system in the system container
  --   with the blackbox solver in double double precision 
  --   and puts the solutions into the solution container,
  --   The two parameters on entry are as follows:
  --     a[0] : 1 or 2 for to be silent or not,
  --     a[1] : the number of tasks.
  --   On return in a[0] is the root count and if not silent,
  --   then a[1] contains the number of characters in b,
  --   where b is the root counter output string.

  function DoblDobl_Laurent_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32;

  -- DESCRIPTION :
  --   Solves the Laurent system in the system container
  --   with the blackbox solver in double double precision 
  --   and puts the solutions into the solution container.
  --   The two parameters on entry are as follows:
  --     a[0] : 1 or 2 for to be silent or not,
  --     a[1] : the number of tasks.
  --   On return in a[0] is the root count and if not silent,
  --   then a[1] contains the number of characters in b,
  --   where b is the root counter output string.

  function QuadDobl_Polynomial_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32;

  -- DESCRIPTION :
  --   Solves the polynomial system in the system container
  --   with the blackbox solver in quad double precision 
  --   and puts the solutions into the solution container,
  --   The two parameters on entry are as follows:
  --     a[0] : 1 or 2 for to be silent or not,
  --     a[1] : the number of tasks.
  --   On return in a[0] is the root count and if not silent,
  --   then a[1] contains the number of characters in b,
  --   where b is the root counter output string.

  function QuadDobl_Laurent_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32;

  -- DESCRIPTION :
  --   Solves the Laurent system in the system container
  --   with the blackbox solver in quad double precision 
  --   and puts the solutions into the solution container.
  --   The two parameters on entry are as follows:
  --     a[0] : 1 or 2 for to be silent or not,
  --     a[1] : the number of tasks.
  --   On return in a[0] is the root count and if not silent,
  --   then a[1] contains the number of characters in b,
  --   where b is the root counter output string.

end Job_Handlers;
