with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Pade_Approximants;

package Test_QuadDobl_Pade_Approximants is

-- DESCRIPTION :
--   Tests rational approximants in quad double precision.

  function quaddobl_log_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the first dim+1 coefficients of the series of log(1+x)
  --   as a vector of range 0..dim as a vector of complex numbers in
  --   quad double precision.

  function quaddobl_invfactorial
             ( n : integer32 )
             return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns 1/n! where n! is the factorial,
  --   stored as a complex number in quad double precision.

  function quaddobl_exp_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of exp(x) at x = 0, as a vector
  --   of complex numbers in quad double precision.

  function quaddobl_sin_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of sin(x) at x = 0, as a vector
  --   of complex numbers in quad double precision.

  function quaddobl_cos_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of cos(x) at x = 0, as a vector
  --   of complex numbers in quad double precision.

  procedure QuadDobl_Line_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Tests the construction of the Pade approximant in quad
  --   double precision for the line x = t.

  procedure QuadDobl_log_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Tests the construction in quad double arithmetic
  --   on the natural logarithm of 1 + x.

  procedure QuadDobl_sin_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Tests the construction in quad double arithmetic
  --   on the series of sin(x) at x = 0.

  procedure QuadDobl_exp_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Tests the construction in quad double arithmetic
  --   on the series of exp(x) at x = 0.

  procedure QuadDobl_cos_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Tests the construction in quad double arithmetic
  --   on the series of cos(x) at x = 0.

  procedure QuadDobl_Pade_Approximation
              ( nbequ,nbsteps : in integer32;
                srv : in QuadDobl_Complex_Series_Vectors.Vector;
                pv : in QuadDobl_Pade_Approximants.Pade_Vector );

  -- DESCRIPTION :
  --   The Pade approximant pv and the series srv are evaluated in 
  --   as many points as the value of nbsteps.

  procedure QuadDobl_Test_Homotopy;

  -- DESCRIPTION :
  --   Stores the test homotopy, starting at x^2 - 1, and with
  --   target system 3*x^2 - 3/2 in the QuadDobl_Homotopy data.
  --   The test homotopy (1-t)*(x^2 - 1) + t*(3*x^2 - 3/2) = 0
  --   expands into x^2 - 1 - t*x^2 + t + t*3*x^2 - 3/2*t = 0
  --   which leads to (1+2*t)*x^2 = 1 + 1/2*t and thus defines
  --   the function x(t) = ((1 + 1/2*t)/(1 + 2*t))^(1/2).

  procedure QuadDobl_Test_Start_Solutions
              ( sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Returns in sols the two start solutions +1 and -1
  --   for the test homotopy, in quad double precision.

  procedure QuadDobl_Test_Case
              ( sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Defines the homotopy for an example test case,
  --   in quad double precision.

  procedure QuadDobl_Pade_Homotopy
              ( numdeg,dendeg,nbeq,nbsteps : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Computes Pade approximations for the solution paths
  --   defined by an artificial parameter homotopy,
  --   in quad double precision.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   nbeq     number of equations;
  --   nbsteps  number of steps;
  --   sols     start solutions for an artificial-parameter homotopy.

  procedure QuadDobl_Homotopy_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Prompts for a target, start system, and start solutions.
  --   Applies Newton's method for a series development of the first
  --   start solution, in quad double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_QuadDobl_Pade_Approximants;
