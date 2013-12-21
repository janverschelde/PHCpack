with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;

package Hybrid_Durand_Kerner is

-- DESCRIPTION :
--   This package offers a basic version of the method of Durand-Kerner,
--   (aka the method of Weierstrass) to find all complex roots of a polynomial
--   in one variable with complex coefficients.  The input polynomial
--   is assumed to be generic enough for this basic version to work well.
--   The implementation is in hybrid standard/multi-precision arithmetic:
--   to find initial approximations, standard arithmetic is used first,
--   and then refined with multi-precision numbers.

  generic 

    with procedure Write ( step : in natural32; z,res : in Vector );

    -- DESCRIPTION :
    --   This routine can write intermediate results after each iteration,
    --   such as the step number, the approximations z and the residuals res.

  procedure Reporting_Durand_Kerner
                 ( p : in Vector; z,res : in out Vector;
                   maxsteps : in natural32; eps : in Floating_Number;
                   nb : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   This routine computes all roots of a given polynomial
  --   in one unknown, applying the method of Durand-Kerner.

  -- ON ENTRY :
  --   p           the polynomial defined by
  --                 p[k] + p[k+1]*x + p[k+2]*x^2 + .. + p[k+n]*x^n,
  --               with k = p'first;
  --   z           initial approximations for the roots;
  --   res         the residuals of the roots;
  --   maxsteps    is the maximum number of steps that are allowed;
  --   eps         the required accuracy.

  -- ON RETURN :
  --   z           the computed roots;
  --   res         the residuals of the roots;
  --   nb          the number of steps;
  --   fail        if true, then failed to reach accuracy within maxsteps.

  procedure Silent_Durand_Kerner
                 ( p : in Vector; z,res : in out Vector;
                   maxsteps : in natural32; eps : in Floating_Number;
                   nb : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   Applies the method of Durand-Kerner, without any intermediate output.
  --   All parameters have the same meaning as the reporting version.

end Hybrid_Durand_Kerner;
