with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;

package Test_Standard_Coeff_Homotopy is

-- DESCRIPTION :
--   Tests the evaluation of (1-t)*f + t*g
--   where monomials of f and g are shared
--   as in a coefficient-parameter homotopy, in double precision.

  function Eval ( p,q : Standard_Complex_Polynomials.Poly;
                  x : Standard_Complex_Vectors.Vector;
                  t : double_float )
                return Standard_Complex_Numbers.Complex_Number; 
  function Eval ( p,q : Standard_Complex_Poly_Systems.Poly_Sys;
                  x : Standard_Complex_Vectors.Vector;
                  t : double_float )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the value of (1-t)*p + t*q, evaluated at x.

  function Diff ( p,q : Standard_Complex_Poly_Systems.Poly_Sys;
                  x : Standard_Complex_Vectors.Vector;
                  t : double_float )
                return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the Jacobian matrix of (1-t)*p + t*q, evaluated at x.

  procedure Eval ( n : in natural32; p,q,h : in Poly;
                   cp,cq,ch : in Standard_Complex_Vectors.Vector;
                   ip,iq : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Evaluates (1-t)*p + t*q at some random t in [0,1]
  --   and at a random vector.

  -- ON ENTRY :
  --   n        number of variables;
  --   p        first polynomial in several variables;
  --   q        second polynomial in several variables;
  --   h        labeled structure of the sum of p and q;
  --   cp       coefficient vector of p;
  --   cq       coefficient vector of q;
  --   ch       coefficient vector of h;
  --   ip       index of the labels of coefficients of p in h:
  --            ip(k) locates in ch the k-th coefficient of p,
  --            in particular: ch(ip(k)) := cp(k) sets h to p;
  --   iq       index of the labels of coefficients of q in h:
  --            ip(k) locates in ch the k-th coefficient of q,
  --            in particular: ch(iq(k)) := cq(k) sets h to q.

  procedure Write_Elements ( A,B : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the elements of the matrices A and B to screen,
  --   next to each other for easy comparison.

  procedure Eval ( n : in natural32;
                   p,q,h : in Standard_Complex_Poly_Systems.Poly_Sys;
                   cp,cq,ch : in Standard_Complex_VecVecs.VecVec;
                   ip,iq : in Standard_Integer_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Evaluates (1-t)*p + t*q at some random t in [0,1]
  --   and at a random vector.

  -- ON ENTRY :
  --   n        number of variables;
  --   p        first polynomial in several variables;
  --   q        second polynomial in several variables;
  --   h        labeled structure of the sum of p and q;
  --   cp       coefficient vector of p;
  --   cq       coefficient vector of q;
  --   ch       coefficient vector of h;
  --   ip       index of the labels of coefficients of p in h:
  --            ip(k) locates in ch the k-th coefficient of p,
  --            in particular: ch(ip(k)) := cp(k) sets h to p;
  --   iq       index of the labels of coefficients of q in h:
  --            ip(k) locates in ch the k-th coefficient of q,
  --            in particular: ch(iq(k)) := cq(k) sets h to q.

  procedure Encapsulated_Eval
              ( n : in natural32;
                p,q : in Standard_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Evaluates (1-t)*p + t*q at some random t in [0,1]
  --   and at a random vector.

  -- ON ENTRY :
  --   n        number of variables;
  --   p        first polynomial in several variables;
  --   q        second polynomial in several variables.

  procedure Standard_Compared_Encapsulated_Eval ( n : in natural32 );

  -- DESCRIPTION :
  --   Evaluates (1-t)*p + t*q at some random t in [0,1]
  --   and at a random vector in standard complex arithmetic.

  -- ON ENTRY :
  --   n        number of variables.

  procedure Test_Evaluation
              ( n : in natural32;
                p,q : in Standard_Complex_Polynomials.Poly );

  -- DESCRIPTION :
  --   Tests the evaluation of (1-t)*p + t*q for two polynomials
  --   p and q in n variables.

  procedure Test_System_Evaluation 
              ( n : in natural32;
                p,q : in Standard_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Tests the evaluation of (1-t)*p + t*q for two polynomial
  --   systems p and q in n variables.

  procedure Interactive_Test;

  -- DESCRIPTION :
  --   Prompts for the number of variables
  --   and then for two polynomials in that number of variables,
  --   before running the evaluation test.

  procedure Random_Test;

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated polynomials.

  procedure Standard_Random_Systems
              ( n : in integer32;
                p,q : out Standard_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Generates two n-dimensional systems p and q,
  --   with standard complex coefficients.

  procedure Standard_Random_Coefficient_Systems
              ( n : in integer32;
                p,q : out Standard_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Generates two n-dimensional systems p and q,
  --   with the same supports.

  procedure Random_System_Test;

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems.

  procedure Standard_Random_Encapsulation_Test;

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems,
  --   using the encapsulated interface.

  procedure Standard_Compared_Encapsulation_Test;

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems,
  --   comparing with the Standard_Homotopy package.

  procedure Standard_Homotopy_Performance ( n,m : natural32 );

  -- DESCRIPTION :
  --   Generates m random values for x and t
  --   for evaluation in the homotopy with timings taken.

  procedure Standard_Coefficient_Homotopy_Performance ( n,m : natural32 );

  -- DESCRIPTION :
  --   Generates m random values for x and t
  --   for evaluation in the coefficient homotopy with timings taken.

  procedure Standard_Performance_Test;

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems,
  --   comparing with the Standard_Homotopy package.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Standard_Coeff_Homotopy;
