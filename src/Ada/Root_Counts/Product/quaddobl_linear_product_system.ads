with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Natural_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Polynomials;       use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;      use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package QuadDobl_Linear_Product_System is

-- DESCRIPTION :
--   This package enables the construction and the solution of a random
--   linear-product system.  Every polynomial in a linear-product system
--   is a product of linear polynomials.

-- CONSTRUCTORS :

  procedure Init ( n : in natural32 );

  -- DESCRIPTION :
  --   The internal data for this package are initialized with n,
  --   the number of equations of the system.

  -- REQUIRED :
  --   This operation must be executed before all other operations.

  procedure Add_Hyperplane
              ( i : in natural32; h : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Add the hyperplane h to the i-the equation of the system.

  -- REQUIRED :
  --   i is in range 1..n and h'range is 0..n,
  --   where n = QuadDobl_Linear_Product_System.Dimension.

  -- ON ENTRY :
  --   i        index of the equation to add the hyperplane to;
  --   h        coefficients of a hyperplane, stored as follows:
  --            h(0) + h(1)*x(1) + h(2)*x(2) + .. + h(n)*x(n).

-- SELECTORS :

  function Dimension return natural32;

  -- DESCRIPTION :
  --   Returns the number of equations in the product system.

  function Number_of_Hyperplanes ( i : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of hyperplanes for the i-th equation.

  function Get_Hyperplane
             ( i,j : in natural32 ) return QuadDobl_Complex_Vectors.Vector;
  function Get_Hyperplane
             ( i,j : in natural32 )
             return QuadDobl_Complex_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   returns the j-th hyperplane for the i-th equation

  generic
    with procedure process ( h : in QuadDobl_Complex_Vectors.Vector;
                             continue : out boolean );
  procedure Enumerate_Hyperplanes ( i : in natural32 );

  -- DESCRIPTION :
  --   Enumerates all hyperplanes in the product for the i-th equation,
  --   each time calling process.
  --   The enumeration stopped when continue is set to false.

-- SOLVERS :

  procedure Linear_System 
              ( s : in Standard_Natural_Vectors.Vector;
                fail : out boolean;
                A : out QuadDobl_Complex_Matrices.Matrix;
                b : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in A and b the matrix and vector to define the linear
  --   system A*x = b, selecting equations according to s.
  --   Failure is reported when the selection is invalid.

  -- ON ENTRY :
  --   s        must have range 1..QuadDobl_Linear_Product_System.Dimension,
  --            s(i) defines the hyperplane for the i-th equation.

  -- ON RETURN :
  --   fail     if the selection in s is invalid;
  --   A        coefficient matrix of the linear system;
  --   b        righthand side vector of the linear system.

  procedure Solve ( s : in Standard_Natural_Vectors.Vector;
                    tol : in double_float; rcond : out quad_double;
                    fail : out boolean;
                    v : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the linear system, where the i-th hyperplane is the
  --   s(i)-th hyperplane of the linear-product system.  Reports
  --   failure if the condition number is less than the tolerance.

  -- ON ENTRY :
  --   s        must have range 1..QuadDobl_Linear_Product_System.Dimension,
  --            s(i) defines the hyperplane for the i-th equation;
  --   tol      tolerance on the inverse of the condition number rcond,
  --            if rcond is less than tol, fail is set to true;

  -- ON RETURN :
  --   rcond    estimate for the inverse of the condition number
  --            of the selected coefficient matrix;
  --   fail     true if the selection in s is invalid, or if rcond < tol,
  --            false otherwise;
  --   v        solution vector to the selected linear system,
  --            only if fail is false.

  procedure Solve ( sols : in out Solution_List; nl : out natural32 );

  -- DESCRIPTION :
  --   The random product system is solved and the solutions are
  --   put in the list sols.
  --   nl is the number of matrices that are factored.
  -- NOTE :
  --   All possible linear systems are factorized using Gaussian
  --   elimination, together with the estimation of the condition
  --   of the matrices.
  --   Systems with a bad condition are not solved.

  procedure Solve ( sols : in out Solution_List; 
                    nl : out natural32; l : in List );

  -- DESCRIPTION :
  --   Cf. Solve, but only those linear systems are factorized,
  --   for which the linear equations are as indicated in the list
  --   of positions l.

-- ENUMERATORS :

  generic

    with procedure Process ( s : in Standard_Natural_Vectors.Vector;
                             cont : out boolean );

  procedure Enumerate_Solutions
               ( tol : in double_float; rc : out natural32 );

  -- DESCRIPTION :
  --   Via incremental row reduction, enumerates all solutions.
  --   After each good combination, the procedure "Process" is called.
  --   The enumeration terminates when Process return cont as false.

  -- ON ENTRY :
  --   tol       tolerance to decide whether a number is zero or not.

  -- ON RETURN :
  --   rc        the root count, total number of solutions.

  function Count_All_Solutions ( tol : double_float ) return natural32;

  -- DESCRIPTION :
  --   Calls the enumerator to count the number of roots.

  procedure Get_First
               ( tol : in double_float;
                 s : out Standard_Natural_Vectors.Vector;
                 fail : out boolean );
  procedure Get_First
               ( file : in file_type; tol : in double_float;
                 s : out Standard_Natural_Vectors.Vector;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Allocates the data needed to invert the enumeration.
  --   Returns the position of the first solution.

  -- ON ENTRY :
  --   file      for writing diagnostics and intermediate results;
  --   tol       used to decide whether a number is zero or not.

  -- ON RETURN :
  --   s         position vector to the first solution if not fail;
  --   fail      true if the the system has no solution.

  procedure Get_Next
               ( tol : in double_float;
                 d : in Standard_Natural_Vectors.Vector;
                 s : in out Standard_Natural_Vectors.Vector;
                 fail : out boolean );
  procedure Get_Next
               ( file : in file_type; tol : in double_float;
                 d : in Standard_Natural_Vectors.Vector;
                 s : in out Standard_Natural_Vectors.Vector;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Returns the next solution following the index vector s.

  -- REQUIRED : Get_First has been executed successfully.

  -- ON ENTRY :
  --   file      for writing diagnostics and intermediate results;
  --   tol       used to decide whether a number is zero or not;
  --   d         degrees of all the equations in the system;
  --   s         position vector for the current solution.

  -- ON RETURN :
  --   s         position vector to the next solution if not fail;
  --   fail      true if the the system has no more solutions.

  procedure Get_Next 
               ( tol : in double_float;
                 s : out Standard_Natural_Vectors.Vector;
                 fail : out boolean );

  -- DESCRIPTION :
  --   This version of Get_Next has memory.

  -- REQUIRED : Get_First has been executed successfully.

  -- ON ENTRY :
  --   tol       used to decide whether a number is zero or not.

  -- ON RETURN :
  --   s         position vector to the next solution if not fail;
  --   fail      true if the the system has no more solutions.

  procedure Get_Clear;

  -- DESCRIPTION :
  --   Clears the data used to invert the enumeration.

-- EXPAND LINEAR-PRODUCT SYSTEM :

  function Polynomial ( h : QuadDobl_Complex_Vectors.Vector ) return Poly;

  -- DESCRIPTION :
  --   Returns the linear polynomial defined by the coefficients in h.

  function Polynomial_System return Poly_Sys;

  -- DESCRIPTION :
  --   A polynomial system is constructed by multiplying all
  --   the hyperplanes from the equations of the random product system.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   Deallocation of all memory space occupied by the product system.

end QuadDobl_Linear_Product_System;
