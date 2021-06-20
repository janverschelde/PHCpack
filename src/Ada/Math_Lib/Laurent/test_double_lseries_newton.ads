with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Double_Lseries_Polynomials;        use Double_Lseries_Polynomials;

package Test_Double_Lseries_Newton is

-- DESCRIPTION :
--   Tests the development of Newton's method on Laurent series.

  procedure Add_Parameter ( tv : in Table_Vector );

  -- DESCRIPTION :
  --   To the first constant term in tv, adds t.
  --   Writes a warning if there is no constant term.

  -- REQUIRED :
  --   The degree in the table vector must be at least one.

  procedure Test_Regular_Newton
              ( p : in Laur_Sys; sol : in Standard_Complex_Vectors.Vector;
                deg : in integer32; verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs Newton's method on the system p,
  --   starting at the leading constants of a regular solution,
  --   on Laurent series where the highest power of t equals deg. 

  procedure Test_Singular_Newton
              ( p,sol : in Laur_Sys; tdx,deg : in integer32;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Runs Newton's method on the system p,
  --   starting at the series defined in sol,
  --   on Laurent series where the highest power of t equals deg. 

  -- REQUIRED : tdx /= 0 and 
  --   p has one more variable than the number of equations.

  -- ON ENTRY :
  --   p       a system with one parameter t;
  --   sol     as many univariate polynomials in t as p'length;
  --   tdx     the index of t in p;
  --   deg     precision of the series;
  --   verbose is the verbose flag.

  procedure Test_Isolated_Start;

  -- DESCRIPTION :
  --   Prompts for a system with start solutions and a degree of t.

  procedure Test_Series_Start;

  -- DESCRIPTION :
  --   Prompts for a system, a series, and a degree of t.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts to select a test.  Runs the test.

end Test_Double_Lseries_Newton;
