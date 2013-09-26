with text_io;                            use text_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;
with Multprec_Complex_Poly_SysFun;       use Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Solutions;         use Multprec_Complex_Solutions;

package Multprec_Residual_Evaluations is

-- DESCRIPTION :
--   This package provides routines to evaluate solutions of polynomial
--   systems using multi-precision arithmetic.

  function Residual ( p_eval : Eval_Poly_Sys; zero : Vector )
                    return Floating_Number;

  -- DESCRIPTION :
  --   Returns maximal component in absolute value of the evaluation
  --   of the zero in the system.

  procedure Residual ( file : in file_type;
                       p_eval : in Eval_Poly_Sys; sol : in Solution );

  -- DESCRIPTION :
  --   Writes the residual of the solution on the file.

  procedure Residuals ( file : in file_type;
                        p_eval : in Eval_Poly_Sys; sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes the residuals of the solutions on the file.

end Multprec_Residual_Evaluations;
