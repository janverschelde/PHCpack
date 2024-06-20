with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
-- with PentDobl_Complex_Numbers;
-- with OctoDobl_Complex_Numbers;
-- with DecaDobl_Complex_Numbers;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with PentDobl_Complex_Series_Vectors;
with PentDobl_Complex_Series_VecVecs;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_Complex_Series_VecVecs;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_SysFun;
with Standard_CSeries_Jaco_Matrices;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_SysFun;
with DoblDobl_CSeries_Jaco_Matrices;
with TripDobl_CSeries_Poly_SysFun;
with TripDobl_CSeries_Jaco_Matrices;
with TripDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_SysFun;
with QuadDobl_CSeries_Jaco_Matrices;
with PentDobl_CSeries_Poly_Systems;
with PentDobl_CSeries_Poly_SysFun;
-- with PentDobl_CSeries_Jaco_Matrices;
with OctoDobl_CSeries_Poly_Systems;
with OctoDobl_CSeries_Poly_SysFun;
-- with OctoDobl_CSeries_Jaco_Matrices;
with DecaDobl_CSeries_Poly_Systems;
with DecaDobl_CSeries_Poly_SysFun;
-- with DecaDobl_CSeries_Jaco_Matrices;

package Power_Series_Methods is

-- DESCRIPTION :
--   A power series method applies Newton's method to a system
--   which has as coefficients power series.
--   The procedures in this packages split in several categories:
--   (1) verbose or not;
--   (2) seven levels of multiple double precision are supported;
--   (3) use LU, QR, SVD, or a general Echelon form;
--   (4) work on one vector or on a vector of vectors of series.

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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
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
  --   in double, double double, triple double, or quad double precision.

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
                f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                info : out integer32; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
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
  --   in double, double double, triple double, or quad double precision.

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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; verbose : in boolean := false;
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
  --   in double, double double, triple double, or quad double precision.

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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
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
  --   in double, double double, triple double, or quad double precision.

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
                f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
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
  --   in double, double double, triple double, or quad double precision.

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
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float;
                evp : out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float;
                evp : out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying Singular Value Decomposition to compute the Newton updates,
  --   in double precision.

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
  --            if 1.0 + rcond = 1.0, then the problem is singular;
  --   evp      if the problem is nonsingular, then evp contains
  --            the solution s evaluated at p.

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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                rcond : out triple_double; verbose : in boolean := false;
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
  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out PentDobl_Complex_Series_Vectors.Vector;
                rcond : out penta_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out PentDobl_Complex_Series_Vectors.Vector;
                rcond : out penta_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out OctoDobl_Complex_Series_Vectors.Vector;
                rcond : out octo_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out OctoDobl_Complex_Series_Vectors.Vector;
                rcond : out octo_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DecaDobl_Complex_Series_Vectors.Vector;
                rcond : out deca_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DecaDobl_Complex_Series_Vectors.Vector;
                rcond : out deca_double; verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying Singular Value Decomposition to compute the Newton updates,
  --   in double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out TripDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out TripDobl_Complex_Numbers.Complex_Number;
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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
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
  --   in double, double double, triple double, or quad double precision.

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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
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
              ( maxdeg,nbrit : in integer32;
                p : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in PentDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in OctoDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DecaDobl_Complex_Series_VecVecs.VecVec;
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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in PentDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in OctoDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );
  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DecaDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs Newton's method on the vector of power series in v.
  --   using the singular value decomposition to compute the updates,
  --   in double, double double, triple double, quad double,
  --   penta double, octo double, or deca double precision.

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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
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
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
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
  --   in double, double double, triple double, or quad double precision.

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
