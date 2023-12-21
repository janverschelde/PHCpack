with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

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

  function Get_Core_Count ( a : C_intarrs.Pointer;
                            vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Gets the number of available cores
  --   and returns its value in the parameter a on return.
  --   This is useful for running a computer at full speed.

  function Standard_Polynomial_Solver
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32;

  -- DESCRIPTION :
  --   Solves the polynomial system in the system container
  --   with the blackbox solver in standard double precision
  --   and puts the solutions into the solution container.
  --   The start system and start solutions are stored
  --   in the data management package.
  --   The two parameters on entry are as follows:
  --     a[0] : 1 or 2 for to be silent or not,
  --     a[1] : the number of tasks.
  --     a[2] : 1 if the focus is on mixed volumes and polyhedral homotopies,
  --            0 for all bounds on the number of solutions.
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
  --   The start system and start solutions are stored
  --   in the data management package.
  --   The two parameters on entry are as follows:
  --     a[0] : 1 or 2 for to be silent or not,
  --     a[1] : the number of tasks.
  --     a[2] : 1 if the focus is on mixed volumes and polyhedral homotopies,
  --            0 for all bounds on the number of solutions.
  --   The focus on mixed volumes and polyhedral homotopies is automatic
  --   if the system is a genuine Laurent system with negative exponents.
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
  --   The start system and start solutions are stored
  --   in the data management package.
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
  --   The start system and start solutions are stored
  --   in the data management package.
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
  --   The start system and start solutions are stored
  --   in the data management package.
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
  --   The start system and start solutions are stored
  --   in the data management package.
  --   The two parameters on entry are as follows:
  --     a[0] : 1 or 2 for to be silent or not,
  --     a[1] : the number of tasks.
  --   On return in a[0] is the root count and if not silent,
  --   then a[1] contains the number of characters in b,
  --   where b is the root counter output string.

  function Get_Gamma_Constant
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves the gamma constant used in the solver.

  -- ON ENTRY :
  --   a        a[0] equals 1, 2, or 4, for double, double double,
  --            or quad double precision.

  -- ON RETURN :
  --   c        c[0] is the real part of the complex gamma constant;
  --            c[1] is the imaginary part of the complex gamma constant.

  function Set_Gamma_Constant
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Retrieves the gamma constant used in the solver.

  -- ON ENTRY :
  --   a        a[0] equals 1, 2, or 4, for double, double double,
  --            or quad double precision;
  --   c        c[0] is the real part of the complex gamma constant;
  --            c[1] is the imaginary part of the complex gamma constant.

  function Mixed_Volume
             ( a : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32;

  -- DESCRIPTION :
  --   Computes the mixed volume of the system in the polynomial
  --   system container, with double precision coefficients.
  --   If there is no system in the polynomial systems container,
  --   then the Laurent system in the Laurent systems container
  --   is used as input, also with double precision coefficients.
  --   On return in a[0] is the mixed volume.
  --   The mixed cells are in the cells container.

  function Stable_Mixed_Volume
             ( a,b : C_intarrs.Pointer; vrblvl : integer32 := 0 )
             return integer32;

  -- DESCRIPTION :
  --   Computes the mixed volume and the stable mixed volume
  --   for the system in the systems container for polynomial
  --   systems with coefficients in double precision.
  --   On return in a[0] is the mixed volume
  --   and in b[0] is the stable mixed volume.
  --   The mixed cells are in the cells container.

  function Standard_Condition_Report
             ( a,b : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   For the system and solutions in the containers for double precision,
  --   computes a condition report.

  -- REQUIRED :
  --   a holds space for 6 integers;
  --   b has space for at least 48 integers;
  --   c has space for three doubles.

  -- ON ENTRY :
  --   a[0]    the maximum number of Newton iterations per solution;
  --   a[1]    verbose flag: 1 if verbose, 0 if not;
  --   a[2]    number of characters in the string for the file,
  --           for the optional output to file (0 if no file);
  --   b       contains as many characters as the value of a[2]
  --           to the define the name of the output file;
  --   c[0]    tolerance on the residual;
  --   c[1]    tolerance on the forward error;
  --   c[2]    tolerance on the inverse condition number for singularity.

  -- ON RETURN :
  --   a[0]    the number of failures;
  --   a[1]    the number of real solutions;
  --   a[2]    the number of nonreal solutions;
  --   a[3]    the number of regular solutions;
  --   a[4]    the number of singular solutions;
  --   a[5]    the number of clustered solutions;
  --   b[0..15] is the frequency table with forward errors;
  --   b[16..31] is the frequency table with condition numbers;
  --   b[32..47] is the frequency table with residuals.

end Job_Handlers;
