with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;

package Start_Flag_Homotopies is

-- DESCRIPTION :
--   The Littlewood-Richardson homotopies start a the solution of
--   a linear system.  The procedures in this package extract the
--   linear system out of the formulation for a Schubert problem,
--   which is given as a polynomial system.  The solving happens
--   in standard double, double double, or quad double precision.

  function Inconsistent
             ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return boolean;
  function Inconsistent
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return boolean;
  function Inconsistent
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the system contains at least one nonzero
  --   polynomial of degree zero.

  function Linear_Equations
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Linear_Equations
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Linear_Equations
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns only those polynomials in p that are linear.
  --   The polynomials on return are shared pointers, not deep copies.

  procedure Coefficients
             ( p : in Standard_Complex_Poly_Systems.Poly_Sys; 
               A : out Standard_Complex_Matrices.Matrix;
               b : out Standard_Complex_Vectors.Vector );
  procedure Coefficients
             ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys; 
               A : out DoblDobl_Complex_Matrices.Matrix;
               b : out DoblDobl_Complex_Vectors.Vector );
  procedure Coefficients
             ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys; 
               A : out QuadDobl_Complex_Matrices.Matrix;
               b : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Given a linear system encoded in the polynomials in p,
  --   returns the coefficient matrix in A and the righthandside vector
  --   of the linear system in b.

  -- REQUIRED : A'range(1) = p'range = b'range
  --   and A'last(2) equals the number of variables in p.

  procedure Solve ( A : in out Standard_Complex_Matrices.Matrix;
                    b : in Standard_Complex_Vectors.Vector;
                    x : out Standard_Complex_Vectors.Vector;
                    res : out double_float );
  procedure Solve ( A : in out DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    x : out DoblDobl_Complex_Vectors.Vector;
                    res : out double_double );
  procedure Solve ( A : in out QuadDobl_Complex_Matrices.Matrix;
                    b : in QuadDobl_Complex_Vectors.Vector;
                    x : out QuadDobl_Complex_Vectors.Vector;
                    res : out quad_double );

  -- DESCRIPTION :
  --   Applies singular value decomposition to A to solve A*x = b,
  --   in standard double, double double, or quad double precision.
  --   Returns in x the solution and in res its residual.

  procedure First_Solution
             ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
               fail : out boolean;
               x : out Standard_Complex_Vectors.Vector;
               res : out double_float );
  procedure First_Solution
             ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               fail : out boolean;
               x : out DoblDobl_Complex_Vectors.Vector;
               res : out double_double );
  procedure First_Solution
             ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               fail : out boolean;
               x : out QuadDobl_Complex_Vectors.Vector;
               res : out quad_double );

  -- DESCRIPTION :
  --   Computes the first solution for flag conditions with an
  --   initial localization pattern that turn into a linear system,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   f       flag conditions represented as a polynomial system.

  -- ON RETURN :
  --   fail    true if no solution could be found, false otherwise;
  --   x       start solution for h if not fail;
  --   res     residual of the start solution.

  procedure Start_Solution
             ( h : in Standard_Complex_Poly_Systems.Poly_Sys;
               fail : out boolean;
               x : out Standard_Complex_Vectors.Vector;
               res : out double_float );
  procedure Start_Solution
             ( h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               fail : out boolean;
               x : out DoblDobl_Complex_Vectors.Vector;
               res : out double_double );
  procedure Start_Solution
             ( h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               fail : out boolean;
               x : out QuadDobl_Complex_Vectors.Vector;
               res : out quad_double );

  -- DESCRIPTION :
  --   Returns a start solution to the homotopy or reports failure,
  --   using standard double, double double, or quad double arithmetic.

  -- REQUIRED : x'range = 1..n, n+1 = number of variables in h.

  -- ON ENTRY :
  --   h       a moving flag homotopy, t is the last variable.

  -- ON RETURN :
  --   fail    true if no solution could be found, false otherwise;
  --   x       start solution for h if not fail;
  --   res     residual of the start solution.

end Start_Flag_Homotopies;
