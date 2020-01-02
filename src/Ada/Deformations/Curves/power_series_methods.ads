with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_SysFun;
with Standard_CSeries_Jaco_Matrices;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_SysFun;
with DoblDobl_CSeries_Jaco_Matrices;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_SysFun;
with QuadDobl_CSeries_Jaco_Matrices;

package Power_Series_Methods is

-- DESCRIPTION :
--   A power series method applies Newton's method to a system
--   which has as coefficients power series.
--   The procedures in this packages split in several categories:
--   (1) whether to be verbose or not;
--   (2) which precision: double, double double, or quad double;
--   (3) whether to use LU, QR, SVD, or a general Echelon form;
--   (4) whether to work on one vector or on a vector of vectors of series;
--   Depending on the choices, there are 2x3x4x2 = 48 procedures.

-- NEWTON on ONE VECTOR of POWER SERIES :

  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying LU factorization to compute the Newton updates,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   s        leading coefficients for a power series solution;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        a power series solution to p, up to some order;
  --   info     return code of lufac on the Jacobian matrix.

-- LU NEWTON ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                f : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                f : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying LU factorization to compute the Newton updates,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   f        coefficient-parameter homotopy for faster evaluation;
  --   c        coefficient vectors of the homotopy;
  --   ejm      coefficient-parameter matrix of all derivatives;
  --   mlt      multiplication factors for the derivatives;
  --   s        leading coefficients for a power series solution;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        a power series solution to p, up to some order;
  --   info     return code of lufac on the Jacobian matrix.

  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                rcond : out double_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                rcond : out double_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying LU factorization to compute the Newton updates,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   s        leading coefficients for a power series solution;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        a power series solution to p, up to some order;
  --   rcond    estimate of the inverse condition number the Jacobian
  --            matrix computed by lufco.

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying QR decomposition to compute the Newton updates,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   s        leading coefficients for a power series solution;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        a power series solution to p, up to some order.

-- QR NEWTON ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                f : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                f : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying QR decomposition to compute the Newton updates,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   s        leading coefficients for a power series solution;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        a power series solution to p, up to some order.

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                rcond : out double_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                rcond : out double_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying Singular Value Decomposition to compute the Newton updates,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   s        leading coefficients for a power series solution;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        a power series solution to p, up to some order;
  --   rcond    the inverse condition number computed from the singular values,
  --            if 1.0 + rcond = 1.0, then the problem is singular.

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                det : out Standard_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                det : out Standard_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                det : out DoblDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                det : out DoblDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out QuadDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out QuadDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying Singular Value Decomposition to compute the Newton updates,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   s        leading coefficients for a power series solution;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        a power series solution to p, up to some order;
  --   rcond    the inverse condition number computed from the singular values,
  --            if 1.0 + rcond = 1.0, then the problem is singular.

-- NEWTON on a VECTOR of VECTOR of POWER SERIES :

  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the vector of power series in v.
  --   using LU factorization to compute the updates,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   v        leading coefficients for power series solutions;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   pause    if verbose and pause, then prompts the user to continue;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        updated power series solutions to p, up to some order.

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the vector of power series in v.
  --   using QR decomposition to compute the updates,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   v        leading coefficients for power series solutions;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   pause    if verbose and pause, then prompts the user to continue;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        updated power series solutions to p, up to some order.

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the vector of power series in v.
  --   using the singular value decomposition to compute the updates,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   v        leading coefficients for power series solutions;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   pause    if verbose and pause, prompts the user to continue;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        updated power series solutions to p, up to some order.

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the vector of power series in v.
  --   using the lower triangular echelon form to compute the updates,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   maxdeg   maximal degree of the series;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   v        leading coefficients for power series solutions;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   pause    if verbose and pause, prompts the user to continue;
  --   vrblvl   verbose level to indicate name of the procedure.

  -- ON RETURN :
  --   s        updated power series solutions to p, up to some order.

end Power_Series_Methods;
