with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions;
with Standard_Complex_Series_Vectors;
with Standard_Pade_Approximants;

package Test_Standard_Pade_Approximants is

-- DESCRIPTION :
--   Tests rational approximants in double precision.

  function standard_log_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the first dim+1 coefficients of the series of log(1+x)
  --   as a vector of range 0..dim as a vector of complex numbers in
  --   standard double precision.


  function standard_invfactorial
             ( n : integer32 )
             return Standard_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns 1/n! where n! is the factorial,
  --   stored as a complex number in standard double precision.

  function standard_exp_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of exp(x) at x = 0, as a vector
  --   of complex numbers in standard double precision.

  function standard_sin_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of sin(x) at x = 0, as a vector
  --   of complex numbers in standard double precision.

  function standard_cos_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..dim with the coefficients
  --   of the series expansion of cos(x) at x = 0, as a vector
  --   of complex numbers in standard double precision.

  procedure Standard_Line_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Tests the construction of the Pade approximant in standard
  --   double precision for the line x = t.

  procedure Standard_log_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the natural logarithm of 1 + x.

  procedure Standard_sin_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the series of sin(x) at x = 0.

  procedure Standard_exp_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the series of exp(x) at x = 0.

  procedure Standard_cos_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the series of cos(x) at x = 0.

  procedure Standard_Pade_Approximation
              ( nbequ,nbsteps : in integer32;
                srv : in Standard_Complex_Series_Vectors.Vector;
                pv : in Standard_Pade_Approximants.Pade_Vector );

  -- DESCRIPTION :
  --   The Pade approximant pv and the series srv are evaluated in 
  --   as many points as the value of nbsteps.

  procedure Standard_Test_Homotopy;

  -- DESCRIPTION :
  --   Stores the test homotopy, starting at x^2 - 1, and with
  --   target system 3*x^2 - 3/2 in the Standard_Homotopy data.
  --   The test homotopy (1-t)*(x^2 - 1) + t*(3*x^2 - 3/2) = 0
  --   expands into x^2 - 1 - t*x^2 + t + t*3*x^2 - 3/2*t = 0
  --   which leads to (1+2*t)*x^2 = 1 + 1/2*t and thus defines
  --   the function x(t) = ((1 + 1/2*t)/(1 + 2*t))^(1/2).

  procedure Standard_Test_Start_Solutions
              ( sols : out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Returns in sols the two start solutions +1 and -1
  --   for the test homotopy, in standard double precision.

  procedure Standard_Test_Case
              ( sols : out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Defines the homotopy for an example test case,
  --   in standard double precision.

  procedure Standard_Pade_Homotopy
              ( numdeg,dendeg,nbeq,nbsteps : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Computes Pade approximations for the solution paths
  --   defined by an artificial parameter homotopy,
  --   in standard double precision.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   nbeq     number of equations;
  --   nbsteps  number of steps;
  --   sols     start solutions for an artificial-parameter homotopy.

  procedure Standard_Homotopy_Test ( numdeg,dendeg : in integer32 );

  -- DESCRIPTION :
  --   Prompts for a target, start system, and start solutions.
  --   Applies Newton's method for a series development of the first
  --   start solution, in standard double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Standard_Pade_Approximants;
