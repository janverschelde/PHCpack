package Test_Multprec_Polynomials is

-- DESCRIPTION :
--   Tests multiprecision complex polynomials.

  procedure Test_Multprec_io;

  -- DESCRIPTION :
  --   Tests the input/output of a polynomial in several variables
  --   and with complex coefficients.

  procedure Test_Multprec_Eval;

  -- DESCRIPTION :
  --   Tests the evaluation of a polynomial in several variables
  --   and with multiprecision complex coefficients.

  procedure Test_Multprec_Diff;

  -- DESCRIPTION :
  --   Test on the differentiation of multiprecision complex polynomials.

  procedure Test_Eval_Multprec_System;

  -- DESCRIPTION :
  --   Tests the evaluation of a multiprecision polynomial system.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Multprec_Polynomials;
