with text_io;                            use text_io;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Natural_VecVecs;           use Standard_Natural_VecVecs;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;

package Monodromy_Polynomial_Breakup is

-- DESCRIPTION :
--   The routines in this package specialize the algorithms for the
--   numerical irreducible decomposition to the problem of factoring
--   complex multivariate polynomials.

-- REQUIRED :
--   The generic points lie on factors of multiplicity one.

  procedure Monodromy_Breakup
                ( file : in file_type;
                  p : in Poly; ep : in Eval_Poly; b,v,gp : in Vector;
                  output : in boolean; deco : in out Link_to_VecVec );
  procedure Monodromy_Breakup
                ( p : in Poly; ep : in Eval_Poly; b,v,gp : in Vector;
                  deco : in out Link_to_VecVec );

  -- DESCRIPTION :
  --   With monodromy loops we connect generic points on the same factor.

  -- ON ENTRY :
  --   file       for diagnostics and intermediate output;
  --   p          multivariate polynomial, number of vars is b'length;
  --   ep         evaluation form of the polynomial p;
  --   b          base point of general line x(t) = b + t*v;
  --   v          direction of the general line;
  --   gp         values for the generic points p(b + gp(i)*v) = 0.
  --   output     if true, then predictor-corrector will give output,
  --              otherwise, the continuation will be silent.
  --   deco       initially, when nothing is known yet about the breakup,
  --              the decomposition is a set of singletons.

  -- ON RETURN :
  --   deco       points to vector of t-values for generic points on the 
  --              same factor; deco'length = number of factors.

end Monodromy_Polynomial_Breakup;
