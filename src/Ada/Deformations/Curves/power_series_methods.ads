with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_VecVecs;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Dense_Series_VecVecs;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_VecVecs;
with Standard_Series_Poly_Systems;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_Systems;

package Power_Series_Methods is

-- DESCRIPTION :
--   A power series method applies Newton's method to a system
--   which has as coefficients power series.
--   The procedures in this packages split in several categories:
--   (1) whether to be verbose or not;
--   (2) which precision: double, double double, or quad double;
--   (3) whether to use LU factorization, QR decomposition, or SVD;
--   (4) whether to work on one vector or on a vector of vectors of series;
--   Depending on the choices, there are 2x3x3x2 = 36 procedures.

-- NEWTON on ONE VECTOR of POWER SERIES :

  procedure Run_LU_Newton
              ( nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                s : in out Standard_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Run_LU_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                s : in out Standard_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Run_LU_Newton
              ( nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Run_LU_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Run_LU_Newton
              ( nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Run_LU_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying LU factorization to compute the Newton updates,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   s        leading coefficients for a power series solution;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output.

  -- ON RETURN :
  --   s        a power series solution to p, up to some order.

  procedure Run_QR_Newton
              ( nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                s : in out Standard_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Run_QR_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                s : in out Standard_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Run_QR_Newton
              ( nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Run_QR_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Run_QR_Newton
              ( nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );
  procedure Run_QR_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Dense_Series_Vectors.Vector;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying QR decomposition to compute the Newton updates,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   s        leading coefficients for a power series solution;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output.

  -- ON RETURN :
  --   s        a power series solution to p, up to some order.

  procedure Run_SVD_Newton
              ( nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                s : in out Standard_Dense_Series_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := false );
  procedure Run_SVD_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                s : in out Standard_Dense_Series_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := false );
  procedure Run_SVD_Newton
              ( nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Dense_Series_Vectors.Vector;
                rcond : out double_double; verbose : in boolean := false );
  procedure Run_SVD_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Dense_Series_Vectors.Vector;
                rcond : out double_double; verbose : in boolean := false );
  procedure Run_SVD_Newton
              ( nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Dense_Series_Vectors.Vector;
                rcond : out quad_double; verbose : in boolean := false );
  procedure Run_SVD_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Dense_Series_Vectors.Vector;
                rcond : out quad_double; verbose : in boolean := false );

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p,
  --   applying Singular Value Decomposition to compute the Newton updates,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   s        leading coefficients for a power series solution;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output.

  -- ON RETURN :
  --   s        a power series solution to p, up to some order;
  --   rcond    the inverse condition number computed from the singular values,
  --            if 1.0 + rcond = 1.0, then the problem is singular.

-- NEWTON on a VECTOR of VECTOR of POWER SERIES :

  procedure Run_LU_Newton
              ( nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false );
  procedure Run_LU_Newton
              ( nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false );
  procedure Run_LU_Newton
              ( nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false );
  procedure Run_LU_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Run_LU_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Run_LU_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs Newton's method on the vector of power series in v.
  --   using LU factorization to compute the updates,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   v        leading coefficients for power series solutions;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   pause   

  -- ON RETURN :
  --   s        updated power series solutions to p, up to some order.

  procedure Run_QR_Newton
              ( nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false );
  procedure Run_QR_Newton
              ( nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false );
  procedure Run_QR_Newton
              ( nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false );
  procedure Run_QR_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Run_QR_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Run_QR_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs Newton's method on the vector of power series in v.
  --   using QR decomposition to compute the updates,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   v        leading coefficients for power series solutions;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   pause   

  -- ON RETURN :
  --   s        updated power series solutions to p, up to some order.

  procedure Run_SVD_Newton
              ( nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false );
  procedure Run_SVD_Newton
              ( nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false );
  procedure Run_SVD_Newton
              ( nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false );
  procedure Run_SVD_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Run_SVD_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false );
  procedure Run_SVD_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false );

  -- DESCRIPTION :
  --   Runs Newton's method on the vector of power series in v.
  --   using the singular value decomposition to compute the updates,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output, to write results,
  --            if not provided, then output is written to screen;
  --   nbrit    number of new iterations;
  --   p        a polynomial system with series coefficients;
  --   v        leading coefficients for power series solutions;
  --   verbose  indicates if results of intermediate Newton steps
  --            need to be written to file or to standard output;
  --   pause   

  -- ON RETURN :
  --   s        updated power series solutions to p, up to some order.

end Power_Series_Methods;
