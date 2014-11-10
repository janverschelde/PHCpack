with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;
with Multprec_Complex_Vectors;
with Multprec_Complex_Polynomials;

package Random_Conditioned_Evaluations is

-- DESCRIPTION :
--   The numerical conditioning of the problem of evaluating a polynomial 
--   is determined by four factors:
--   (1) the size of the coefficients of the polynomial;
--   (2) the size of the coordinates of the point where to evaluate;
--   (3) the largest degree of the monomials in the polynomial;
--   (4) the distance of the point to a root of the polynomial.
--   Based on these factors, we can predict the condition number:
--
--   (coefficient size)*(size of coordinates of point)**(largest degree)
--   -------------------------------------------------------------------.
--              distance of the point to the closest root
--   
--   Once the condition number becomes respectively larger than 1.0E+16, 
--   1.0E+32, and 1.0E+64, then the standard double, double double, 
--   and quad double precision will no longer suffice.

  procedure Random_Conditioned_Evaluation_Problem
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close : in double_float;
                f : out Standard_Complex_Polynomials.Poly;
                z : out Standard_Complex_Vectors.Vector );
  procedure Random_Conditioned_Evaluation_Problem
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close : in double_float;
                f : out DoblDobl_Complex_Polynomials.Poly;
                z : out DoblDobl_Complex_Vectors.Vector );
  procedure Random_Conditioned_Evaluation_Problem
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close : in double_float;
                f : out QuadDobl_Complex_Polynomials.Poly;
                z : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates a polynomial evaluation problem with its numerical
  --   conditioning controlled by three parameters: the size of the
  --   coefficients of the polynomial, the size of the coordinates
  --   of the point where to evaluate, and the closeness of the point
  --   to a root of the polynomial.  Several levels of precision are
  --   supported: standard double, double double, and quad double.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root.

  -- ON RETURN :
  --   f        polynomial in n variables of degrees at most d,
  --            with coefficients of the given cffsz;
  --   z        coordinates of a point of the given pntsz.

  procedure Random_Conditioned_Evaluation_Problem
              ( n,d,m,c,sz : in natural32;
                cffsz,pntsz,close : in double_float;
                f : out Multprec_Complex_Polynomials.Poly;
                z : out Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates a polynomial evaluation problem with its numerical
  --   conditioning controlled by three parameters: the size of the
  --   coefficients of the polynomial, the size of the coordinates
  --   of the point where to evaluate, and the closeness of the point
  --   to a root of the polynomial.  
  --   The computations are done in arbitrary multiprecision arithmetic.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   sz       size of the numbers determines the working precision;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root.

  -- ON RETURN :
  --   f        polynomial in n variables of degrees at most d,
  --            with coefficients of the given cffsz;
  --   z        coordinates of a point of the given pntsz.

  procedure Fix_Gradient
              ( f : in out Standard_Complex_Polynomials.Poly;
                g,z : in Standard_Complex_Vectors.Vector );
  procedure Fix_Gradient
              ( f : in out DoblDobl_Complex_Polynomials.Poly;
                g,z : in DoblDobl_Complex_Vectors.Vector );
  procedure Fix_Gradient
              ( f : in out QuadDobl_Complex_Polynomials.Poly;
                g,z : in QuadDobl_Complex_Vectors.Vector );
  procedure Fix_Gradient
              ( f : in out Multprec_Complex_Polynomials.Poly;
                g,z : in Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Adjusts the linear terms of the polynomial in f, so that the
  --   gradient of f at z takes the corresponding values in g,
  --   in particular: Eval(Diff(f,k),z) = g(k), for k in z'range.
  --   On return, the k-th partial derivative of f will evaluate
  --   to g(k) at the point z.

  -- REQUIRED : z'range = g'range = 1..Number_of_Unknowns(f).

  procedure Random_Conditioned_Gradient_Evaluation
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close : in double_float;
                g : in Standard_Complex_Vectors.Vector;
                f : out Standard_Complex_Polynomials.Poly;
                z : out Standard_Complex_Vectors.Vector );
  procedure Random_Conditioned_Gradient_Evaluation
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close : in double_float;
                g : in DoblDobl_Complex_Vectors.Vector;
                f : out DoblDobl_Complex_Polynomials.Poly;
                z : out DoblDobl_Complex_Vectors.Vector );
  procedure Random_Conditioned_Gradient_Evaluation
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close : in double_float;
                g : in QuadDobl_Complex_Vectors.Vector;
                f : out QuadDobl_Complex_Polynomials.Poly;
                z : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates a conditioned evaluation problem with given gradient.
  --   Generating polynomials with prescribed gradient, we can set the
  --   condition number of the Jacobian matrix of a polynomial system.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   g        n values for the gradient of f at z.

  -- ON RETURN :
  --   f        polynomial in n variables of degrees at most d,
  --            with coefficients of the given cffsz,
  --            and with values for its k-th partial derivative at z
  --            equal to the value of g(k);
  --   z        coordinates of a point of the given pntsz.

  procedure Random_Conditioned_Gradient_Evaluation
              ( n,d,m,c,sz : in natural32;
                cffsz,pntsz,close : in double_float;
                g : in Multprec_Complex_Vectors.Vector;
                f : out Multprec_Complex_Polynomials.Poly;
                z : out Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates a polynomial evaluation problem with its numerical
  --   conditioning controlled by three parameters: the size of the
  --   coefficients of the polynomial, the size of the coordinates
  --   of the point where to evaluate, and the closeness of the point
  --   to a root of the polynomial.  
  --   The computations are done in arbitrary multiprecision arithmetic.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   sz       size of the numbers determines the working precision;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   g        n values for the gradient of f at z.

  -- ON RETURN :
  --   f        polynomial in n variables of degrees at most d,
  --            with coefficients of the given cffsz;
  --   z        coordinates of a point of the given pntsz.

end Random_Conditioned_Evaluations;
