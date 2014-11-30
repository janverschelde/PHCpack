with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;

package Random_Conditioned_Root_Problems is

-- DESCRIPTION :
--   The numerical conditioning of a root problem is determined by the
--   condition number of Jacobian matrix at the root and the numerical
--   condition number of the polynomial evaluation problem.
--   The procedures in this package offer to generate random root problems
--   with prescribed condition numbers.

  procedure Standard_Test
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                z : in out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps on p,
  --   in standard double precision, starting at z
  --   as the initial approximation for a solution of p.

  procedure DoblDobl_Test
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                z : in out DoblDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in double double precision, starting at z 
  --   as an initial approximation for a solution of p.

  procedure QuadDobl_Test
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                z : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in quad double precision, starting at z
  --   as an initial approximation for a solution.

  procedure Multprec_Test
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                z : in out Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Runs a sequence of Newton steps using p,
  --   in arbitrary multiprecision with numbers of the given size,
  --   using z as initial approximation for a solution of p.

  procedure Standard_Conditioned_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float );

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with standard complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+16.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

  procedure DoblDobl_Conditioned_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float );

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with double double complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+32.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

  procedure QuadDobl_Conditioned_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float );

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with quad double complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+64.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

  procedure Multprec_Conditioned_Test
              ( n,d,m,c,sz : in natural32;
                cffsz,pntsz,close,condjm : in double_float );

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with multiprecision complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed, provided the value sz for the precision
  --   is sufficiently high.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   sz       size of the numbers in the working precision;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

  procedure Multprec_Conditioned_Root_Problem
              ( n,d,m,c,prcn : in natural32;
                cffsz,pntsz,close,condjm : in double_float;
                f : out Link_to_Array_of_Strings; z : out Link_to_String );

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with multiprecision complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.

  -- REQUIRED :
  --   The working precision prcn is sufficiently high to compute the
  --   prescribed conditioning of the root problem.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   prcn     number of decimal places in the working precision;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

  -- ON RETURN :
  --   f        string representation of a polynomial system;
  --   z        string representation of an initial approximation.

  procedure Random_Conditioned_Parameters
              ( n,d,m : out natural32;
                condjm,cffsize,pntsize,close,condfz : out double_float );

  -- DESCRIPTION :
  --   Prompts the user for the parameters that determined the numerical
  --   conditioning of a root problem.

  -- ON RETURN :
  --   n        number of variables and equations in the polynomial system;
  --   d        largest degree of a monomial;
  --   m        number of monomials per equation (0 for dense);
  --   condjm   condition number of the Jacobian matrix;
  --   cffsize  size of the coefficients of the polynomials;
  --   pntsize  size of the point where to start Newton's method;
  --   close    closeness to a solution of the system;
  --   condfz   condition of the evaluation problem.

  procedure Random_Conditioned_Root_Problem ( preclvl : in character );

  -- DESCRIPTION :
  --   Prompts the user for the dimensions of the numerical problem
  --   to set the condition of the Jacobian matrix and
  --   the condition of the polynomial evaluation problem.
  --   The level of precision is indicated by the character preclvl.

  function Maximum ( a,b : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the maximum of a and b.
 
  procedure Random_Conditioned_Root_Problem
              ( f : out Link_to_Array_of_Strings;
                z : out Link_to_String );

  -- DESCRIPTION :
  --   Prompts the user for the parameters of a conditioned root problem
  --   and after generating the problem with sufficiently high precision,
  --   returns the problem in the strings f and z.

  -- ON RETURN :
  --   f        string representation of a polynomial system;
  --   z        string representation of initial approximation for a root.

end Random_Conditioned_Root_Problems;
