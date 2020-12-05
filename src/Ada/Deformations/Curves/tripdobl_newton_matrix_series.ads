with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Triple_Double_Numbers;             use Triple_Double_Numbers;
with TripDobl_Complex_Numbers;          use TripDobl_Complex_Numbers;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Complex_Series_VecVecs;
with TripDobl_CSeries_Poly_Systems;
with TripDobl_CSeries_Poly_SysFun;
with TripDobl_CSeries_Jaco_Matrices;

package TripDobl_Newton_Matrix_Series is

-- DESCRIPTION :
--   This package contains the application of Newton's method to compute
--   truncated power series approximations for a space curve.
--   The working precision is triple double precision.
--   The solving of the linear systems uses matrix series rather than
--   series matrices, applying linearization on the Jacobian matrix.

-- ONE NEWTON STEP WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );

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
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

-- ONE NEWTON STEP WITH LU ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure LU_Newton_Step
              ( f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Performs one step with Newton's method on the square system p,
  --   starting at the series approximation x, 
  --   calculating with power series up to the given degree,
  --   using LU factorization to solve the linear system.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, LU_Newton_Step is silent;
  --   f        coefficient-parameter homotopy for evaluation;
  --   c        coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   degree   the degree at which to solve the linear system;
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; vrblvl : in integer32 := 0 );

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
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   x        updated approximation for the series solution;
  --   rcond    estimate for the inverse of the condition number,
  --            if close to zero, then the Jacobian matrix at x is
  --            ill conditioned and x may be wrongly updated.

-- ONE NEWTON STEP WITH QR :

  procedure QR_Newton_Step
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure QR_Newton_Step
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure QR_Newton_Step
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure QR_Newton_Step
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );

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
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

-- ONE NEWTON STEP WITH QR ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure QR_Newton_Step
              ( f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure QR_Newton_Step
              ( file : in file_type;
                f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Performs one step with Newton's method on the system p,
  --   starting at the series approximation x, 
  --   calculating with power series up to the given degree,
  --   using QR decomposition to solve the linear system.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, LU_Newton_Step is silent;
  --   f        coefficient-parameter homotopy for evaluation;
  --   c        coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   degree   the degree at which to solve the linear system;
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

-- ONE NEWTON STEP WITH SVD :

  procedure SVD_Newton_Step
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out triple_double;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Step
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out triple_double;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Step
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out triple_double;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Step
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out triple_double;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Performs one step with Newton's method on the system p,
  --   starting at the series approximation x, 
  --   calculating with power series up to the given degree,
  --   using the SVD and least squares to solve the linear system.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, LU_Newton_Step is silent;
  --   p        a polynomial system with series coefficients;
  --   jp       Jacobi matrix of the system p;
  --   degree   the degree at which to solve the linear system;
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   x        updated approximation for the series solution;
  --   info     return code of the SVD on the lead coefficient matrix;
  --   rcond    inverse condition number of the lead coeffficient matrix,
  --            computed with singular values.

-- ONE NEWTON STEP WITH ECHELON FORM :

  procedure Echelon_Newton_Step
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 );
  procedure Echelon_Newton_Step
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 );
  procedure Echelon_Newton_Step
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 );
  procedure Echelon_Newton_Step
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Performs one step with Newton's method on the system p,
  --   starting at the series approximation x, 
  --   calculating with power series up to the given degree,
  --   using a lower triangular echelon form to solve the linear system.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, LU_Newton_Step is silent;
  --   p        a polynomial system with series coefficients;
  --   jp       Jacobi matrix of the system p;
  --   degree   the degree at which to solve the linear system;
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   x        updated approximation for the series solution;
  --   det      determinant of the lower triangular echelon form,
  --            if nonzero, then the problem is regular,
  --            otherwise, the solution may be a formal Laurent series.

-- MANY NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );

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
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   degree   last degree of the computation;
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

-- MANY NEWTON STEPS WITH LU ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure LU_Newton_Steps
              ( f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does a number of Newton steps on the square system p,
  --   starting at x, doubling the degree after each step,
  --   with LU factorization on the Jacobian matrix,
  --   terminating if info /= 0 or if nbrit is reached.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, LU_Newton_Step is silent;
  --   f        coefficient-parameter homotopy for evaluation;
  --   c        coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   degree   the degree at start of the computations;
  --   nbrit    total number of Newton steps;
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   degree   last degree of the computation;
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

-- MANY NEWTON STEPS WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; vrblvl : in integer32 := 0 );

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
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   degree   last degree of the computation;
  --   x        updated approximation for the series solution;
  --   rcond    estimate for the inverse of the condition number,
  --            if close to zero, then the Jacobian matrix at x is
  --            ill conditioned and x may be wrongly updated.

-- MANY NEWTON STEPS WITH QR :

  procedure QR_Newton_Steps
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure QR_Newton_Steps
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does a number of Newton steps on the system p,
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
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   degree   last degree of the computation;
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

-- MANY NEWTON STEPS WITH QR ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure QR_Newton_Steps
              ( f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );
  procedure QR_Newton_Steps
              ( file : in file_type;
                f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does a number of Newton steps on the system p,
  --   starting at x, doubling the degree after each step,
  --   with QR decomposition on the Jacobian matrix,
  --   terminating if info /= 0 or if nbrit is reached.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, QR_Newton_Step is silent;
  --   f        coefficient-parameter homotopy for evaluation;
  --   c        coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   degree   the degree at start of the computations;
  --   maxdeg   the maximal degree of the series in the steps;
  --   nbrit    total number of Newton steps;
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   degree   last degree of the computation;
  --   x        updated approximation for the series solution;
  --   info     if zero, then the Jacobian matrix at x is regular,
  --            otherwise, info indicates the column at which the
  --            pivoting failed to find an invertible element.

-- MANY NEWTON STEPS WITH SVD :

  procedure SVD_Newton_Steps
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out triple_double;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Steps
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out triple_double;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Steps
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out triple_double;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Steps
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out triple_double;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does a number of Newton steps on the system p,
  --   starting at x, doubling the degree after each step,
  --   with SVD decomposition on the Jacobian matrix,
  --   terminating if info /= 0 or if nbrit is reached.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, QR_Newton_Step is silent;
  --   p        a polynomial system with series coefficients;
  --   jp       Jacobi matrix of the system p;
  --   degree   the degree at start of the computations;
  --   maxdeg   the maximal degree of the series in the steps;
  --   nbrit    total number of Newton steps;
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   degree   last degree of the computation;
  --   x        updated approximation for the series solution;
  --   info     return code of the SVD on the lead coefficient matrix;
  --   rcond    inverse condition number of the lead coeffficient matrix,
  --            computed with singular values.

-- MANY NEWTON STEPS WITH ECHELON FORM :

  procedure Echelon_Newton_Steps
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 );
  procedure Echelon_Newton_Steps
              ( p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 );
  procedure Echelon_Newton_Steps
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 );
  procedure Echelon_Newton_Steps
              ( file : in file_type;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Does a number of Newton steps on the square system p,
  --   starting at x, doubling the degree after each step,
  --   with an echelon form on the Jacobian matrix,
  --   terminating if nbrit is reached.

  -- ON ENTRY :
  --   file     for intermediate output: p(x) and the update dx,
  --            if omitted, LU_Newton_Step is silent;
  --   p        a polynomial system with series coefficients;
  --   jp       Jacobi matrix of the system p;
  --   degree   the degree at start of the computations;
  --   nbrit    total number of Newton steps;
  --   x        current approximation for the series solution;
  --   vrblvl   if > 0, then the name of the procedure is written.

  -- ON RETURN :
  --   degree   last degree of the computation;
  --   x        updated approximation for the series solution;
  --   det      dterminant of the echelon Hermite-Laurent matrix
  --            if nonzero, then there are no nonzero coefficients
  --            in the series with negative exponents, otherwise,
  --            the solution is a formal Laurent series.

end TripDobl_Newton_Matrix_Series;
