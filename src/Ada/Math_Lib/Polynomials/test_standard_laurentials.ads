with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;

package Test_Standard_Laurentials is

-- DESCRIPTION :
--   Tests complex Laurent polynomials in double precision.

  procedure Test_Standard_Laurent_Eval
              ( p : in Standard_Complex_Laurentials.Poly;
                e : in Standard_Complex_Laur_Functions.Eval_Poly;
                x : in Standard_Complex_Vectors.Vector;
                output_of_results : in boolean; bug : out boolean );

  -- DESCRIPTION :
  --   Evaluates the polynomial twice and compares the results.

  procedure Interactive_Standard_Laurent_Eval;

  -- DESCRIPTION :
  --   Tests the evaluation of a polynomial in several variables
  --   and with standard complex coefficients.

  procedure Random_Standard_Laurent_Eval;

  -- DESCRIPTION :
  --   Tests the evaluation of a Laurent polynomial in several variables
  --   and with standard complex coefficients.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Standard_Laurentials;
