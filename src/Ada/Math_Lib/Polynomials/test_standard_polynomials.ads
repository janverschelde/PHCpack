with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;

package Test_Standard_Polynomials is

-- DESCRIPTION :
--   Tests complex polynomials in double precision.

  procedure Test_Standard_io;

  -- DESCRIPTION :
  --   Tests the input/output of a polynomial in several variables
  --   and with complex coefficients.

  procedure Test_Vector_io;

  -- DESCRIPTION :
  --   Tests the input and output of vectors of polynomials.

  procedure Test_Matrix_io;

  -- DESCRIPTION :
  --   Tests the input and output of matrices of polynomials.

  procedure Test_Standard_Eval
              ( p : in Standard_Complex_Polynomials.Poly;
                e : in Standard_Complex_Poly_Functions.Eval_Poly;
                x : in Standard_Complex_Vectors.Vector;
                output_of_results : in boolean; bug : out boolean );

  -- DESCRIPTION :
  --   Evaluates the polynomial twice and compares the results.

  procedure Interactive_Standard_Eval;

  -- DESCRIPTION :
  --   Tests the evaluation of a polynomial in several variables
  --   and with standard complex coefficients.

  procedure Random_Standard_Eval;

  -- DESCRIPTION :
  --   Tests the evaluation of a polynomial in several variables
  --   and with standard complex coefficients.

  procedure Test_Standard_Diff;

  -- DESCRIPTION :
  --   Test on the differentiation of standard complex polynomials.

  procedure Test_System_io;

  -- DESCRIPTION :
  --   Tests input and output of polynomial systems.

  procedure Test_Eval_Standard_System;

  -- DESCRIPTION :
  --   Prompts for a polynomial system and tests the evaluation.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Standard_Polynomials;
