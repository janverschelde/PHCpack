with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Jaco_Matrices;

package QuadDobl_Newton_Series is

-- DESCRIPTION :
--   This package contains the application of Newton's method to compute
--   truncated power series approximations for a space curve.
--   The working precision is quad double precision.

-- ONE NEWTON STEP :

  procedure LU_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure LU_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure LU_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure LU_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );

  -- DESCRIPTION :
  --   Performs one step with Newton's method on the square system p,
  --   starting at the series approximation x, 
  --   calculating with power series up to the given degree,
  --   using LU factorization to solve the linear system.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, LU_Newton_Step is silent;
  --   p        a polynomial system with series coefficients;
  --   jp       Jacobi matrix of the system p;
  --   degree   the degree at which to solve the linear system;
  --   x        current approximation for the series solution.

  -- ON RETURN :
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

  procedure QR_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure QR_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure QR_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure QR_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );

  -- DESCRIPTION :
  --   Performs one step with Newton's method on the system p,
  --   starting at the series approximation x, 
  --   calculating with power series up to the given degree,
  --   using QR decomposition to solve the linear system.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, LU_Newton_Step is silent;
  --   p        a polynomial system with series coefficients;
  --   jp       Jacobi matrix of the system p;
  --   degree   the degree at which to solve the linear system;
  --   x        current approximation for the series solution.

  -- ON RETURN :
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

-- MANY NEWTON STEPS :

  procedure LU_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure LU_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );

  -- DESCRIPTION :
  --   Does a number of Newton steps on the square system p,
  --   starting at x, doubling the degree after each step,
  --   with LU factorization on the Jacobian matrix,
  --   terminating if info /= 0 or if nbrit is reached.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, LU_Newton_Step is silent;
  --   p        a polynomial system with series coefficients;
  --   jp       Jacobi matrix of the system p;
  --   degree   the degree at start of the computations;
  --   maxdeg   the maximal degree of the series in the steps;
  --   nbrit    total number of Newton steps;
  --   x        current approximation for the series solution.

  -- ON RETURN :
  --   degree   last degree of the computation;
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

  procedure QR_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure QR_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );
  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32 );

  -- DESCRIPTION :
  --   Does a number of Newton steps on the square system p,
  --   starting at x, doubling the degree after each step,
  --   with QR decomposition on the Jacobian matrix,
  --   terminating if info /= 0 or if nbrit is reached.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, QR_Newton_Step is silent;
  --   p        a polynomial system with series coefficients;
  --   jp       Jacobi matrix of the system p;
  --   degree   the degree at start of the computations;
  --   maxdeg   the maximal degree of the series in the steps;
  --   nbrit    total number of Newton steps;
  --   x        current approximation for the series solution.

  -- ON RETURN :
  --   degree   last degree of the computation;
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

end QuadDobl_Newton_Series;
