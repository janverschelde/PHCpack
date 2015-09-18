with text_io;                            use text_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Functions;

package Monodromy_Polynomial_Breakup is

-- DESCRIPTION :
--   The routines in this package specialize the algorithms for the
--   numerical irreducible decomposition to the problem of factoring
--   complex multivariate polynomials, with or without intermediate output,
--   in standard double, double double, or quad double precision.

-- REQUIRED :
--   The generic points lie on factors of multiplicity one.

  procedure Monodromy_Breakup
                ( file : in file_type;
                  p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in Standard_Complex_Vectors.Vector;
                  output : in boolean;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Monodromy_Breakup
                ( file : in file_type;
                  p : in DoblDobl_Complex_Polynomials.Poly;
                  ep : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in DoblDobl_Complex_Vectors.Vector;
                  output : in boolean;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Monodromy_Breakup
                ( file : in file_type;
                  p : in QuadDobl_Complex_Polynomials.Poly;
                  ep : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in QuadDobl_Complex_Vectors.Vector;
                  output : in boolean;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Monodromy_Breakup
                ( p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in Standard_Complex_Vectors.Vector;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Monodromy_Breakup
                ( p : in DoblDobl_Complex_Polynomials.Poly;
                  ep : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in DoblDobl_Complex_Vectors.Vector;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Monodromy_Breakup
                ( p : in QuadDobl_Complex_Polynomials.Poly;
                  ep : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in QuadDobl_Complex_Vectors.Vector;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   With monodromy loops we connect generic points on the same factor.

  -- ON ENTRY :
  --   file       for diagnostics and intermediate output,
  --              if omitted, then the procedure remains silent;
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
