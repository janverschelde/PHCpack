with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Coefficient_Convolutions;

package Staggered_Newton_Convolutions is

-- DESCRIPTION :
--   Several steps with Newton's method on power series are performed
--   with convolutions on the coefficients and linearization.
--   The iterative methods are staggered in increasing degrees of the
--   power series, starting with degrees 1, 2, 4, 8, etc.
--   All Newton steps on coefficient convolution circuits
--   are done in standard double precision.
--   For LU factorizations on the lead matrix, there are 3 versions :
--   (1) using complex arithmetic for the LU factorization;
--   (2) on columns of real and imaginary parts, inlined complex;
--   (3) inlined and indexed, avoiding to recompute coefficients.

-- INDEXED NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure Indexed_LU_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; idxtoldx : out integer32;
                absdx : out double_float; fail : out boolean;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure Indexed_LU_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; idxtoldx : out integer32;
                absdx : out double_float; fail : out boolean;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf, in double precision,
  --   using inlined LU factorization to solve the linear series systems
  --   and the tolerance index is used to avoid the recomputation of
  --   already sufficiently accurate coefficient vectors of the solution.

  -- REQUIRED :
  --   rc'range = ic'range = 1..csr.dim and for all k in rc'range:
  --   rc(k)'range = ic(k)'range = 1..csr.dim.
  --   rv'range = iv'range = 1..csr.deg and for all k in rv'range:
  --   rv(k)'range = iv(k)'range = 1..csr.dim and for all i in rv(k)'range:
  --   rv(k)(i)'range = iv(k)(i)'range = 1..csr.dim.
  --   rb'range = ib'range = 0..csr.deg and for all k in rb'range:
  --   rb(k)'range = ib(k)'range = 1..csr.dim.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   rx       work space for the real parts of the coefficients;
  --   ix       work space for the imaginary parts of the coefficients;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on absdx, used as stop criterium,
  --            and for the computation of the tolerance index;
  --   rc       work space for real parts of the columns of A(0);
  --   ic       work space for imaginary parts of the columns of A(0);
  --   rv       work space for all real parts of all A(k), for k in 1..degree;
  --   iv       work space for all imag parts of all A(k), for k in 1..degree;
  --   rb       real parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ib       imaginary parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ry       allocated work space vector of range 1..s.dim;
  --   iy       allocated work space vector of range 1..s.dim;
  --   scale    if true, then the k-th component of the update dx
  --            is divided by k!, otherwise no scaling to dx is applied;
  --   verbose  if verbose, then the progress index will be written;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   nbrit    number of iterations done;
  --   idxtoldx is the highest index in the update with component less than
  --            the given tolerance tol, the tolerance index;
  --   absdx    absolute value of the update to the last coefficient;
  --   fail     true if absdx > tol after nbrit iterations;
  --   info     info from the LU factorization;
  --   ipvt     pivoting of the LU factorization on the lead matrix;
  --   rc       real parts of the output of lufac on A(0);
  --   ic       imaginary parts of the output of lufac on A(0);
  --   rb       rb(k) stores the real parts of the solution x(k);
  --   ib       ib(k) stores the imaginary parts of the solution x(k).

-- INDEXED NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure Indexed_LU_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; idxtoldx : out integer32;
                absdx : out double_float; fail : out boolean;
                rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure Indexed_LU_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; idxtoldx : out integer32;
                absdx : out double_float; fail : out boolean;
                rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf, in double precision,
  --   using inlined LU factorization to solve the linear series systems
  --   and the tolerance index is used to avoid the recomputation of
  --   already sufficiently accurate coefficient vectors of the solution.

  -- REQUIRED :
  --   rc'range = ic'range = 1..csr.dim and for all k in rc'range:
  --   rc(k)'range = ic(k)'range = 1..csr.dim.
  --   rv'range = iv'range = 1..csr.deg and for all k in rv'range:
  --   rv(k)'range = iv(k)'range = 1..csr.dim and for all i in rv(k)'range:
  --   rv(k)(i)'range = iv(k)(i)'range = 1..csr.dim.
  --   rb'range = ib'range = 0..csr.deg and for all k in rb'range:
  --   rb(k)'range = ib(k)'range = 1..csr.dim.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   rx       work space for the real parts of the coefficients;
  --   ix       work space for the imaginary parts of the coefficients;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on absdx, used as stop criterium,
  --            and for the computation of the tolerance index;
  --   rc       work space for real parts of the columns of A(0);
  --   ic       work space for imaginary parts of the columns of A(0);
  --   rv       work space for all real parts of all A(k), for k in 1..degree;
  --   iv       work space for all imag parts of all A(k), for k in 1..degree;
  --   rb       real parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ib       imaginary parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ry       allocated work space vector of range 1..s.dim;
  --   iy       allocated work space vector of range 1..s.dim;
  --   scale    if true, then the k-th component of the update dx
  --            is divided by k!, otherwise no scaling to dx is applied;
  --   verbose  if verbose, then the progress index will be written;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   nbrit    number of iterations done;
  --   idxtoldx is the highest index in the update with component less than
  --            the given tolerance tol, the tolerance index;
  --   absdx    absolute value of the update to the last coefficient;
  --   fail     true if absdx > tol after nbrit iterations;
  --   rcond    estimate for the inverse condition number of A(0);
  --   ipvt     pivoting of the LU factorization on the lead matrix;
  --   rc       real parts of the output of lufac on A(0);
  --   ic       imaginary parts of the output of lufac on A(0);
  --   rb       rb(k) stores the real parts of the solution x(k);
  --   ib       ib(k) stores the imaginary parts of the solution x(k).

-- INLINED NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure Inlined_LU_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure Inlined_LU_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf, in double precision,
  --   using inlined LU factorization to solve the linear series systems.

  -- REQUIRED :
  --   rc'range = ic'range = 1..csr.dim and for all k in rc'range:
  --   rc(k)'range = ic(k)'range = 1..csr.dim.
  --   rv'range = iv'range = 1..csr.deg and for all k in rv'range:
  --   rv(k)'range = iv(k)'range = 1..csr.dim and for all i in rv(k)'range:
  --   rv(k)(i)'range = iv(k)(i)'range = 1..csr.dim.
  --   rb'range = ib'range = 0..csr.deg and for all k in rb'range:
  --   rb(k)'range = ib(k)'range = 1..csr.dim.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   rx       work space for the real parts of the coefficients;
  --   ix       work space for the imaginary parts of the coefficients;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on absdx, used as stop criterium;
  --   rc       work space for real parts of the columns of A(0);
  --   ic       work space for imaginary parts of the columns of A(0);
  --   rv       work space for all real parts of all A(k), for k in 1..degree;
  --   iv       work space for all imag parts of all A(k), for k in 1..degree;
  --   rb       real parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ib       imaginary parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ry       allocated work space vector of range 1..s.dim;
  --   iy       allocated work space vector of range 1..s.dim;
  --   scale    if true, then the k-th component of the update dx
  --            is divided by k!, otherwise no scaling to dx is applied;
  --   verbose  if verbose, then the progress index will be written;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the update to the last coefficient;
  --   fail     true if absdx > tol after nbrit iterations;
  --   info     info from the LU factorization;
  --   ipvt     pivoting of the LU factorization on the lead matrix;
  --   rc       real parts of the output of lufac on A(0);
  --   ic       imaginary parts of the output of lufac on A(0);
  --   rb       rb(k) stores the real parts of the solution x(k);
  --   ib       ib(k) stores the imaginary parts of the solution x(k).

-- INLINED NEWTON STEPS WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure Inlined_LU_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure Inlined_LU_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf, in double precision,
  --   using inlined LU factorization to solve the linear series systems.

  -- REQUIRED :
  --   rc'range = ic'range = 1..csr.dim and for all k in rc'range:
  --   rc(k)'range = ic(k)'range = 1..csr.dim.
  --   rv'range = iv'range = 1..csr.deg and for all k in rv'range:
  --   rv(k)'range = iv(k)'range = 1..csr.dim and for all i in rv(k)'range:
  --   rv(k)(i)'range = iv(k)(i)'range = 1..csr.dim.
  --   rb'range = ib'range = 0..csr.deg and for all k in rb'range:
  --   rb(k)'range = ib(k)'range = 1..csr.dim.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   rx       work space for the real parts of the coefficients;
  --   ix       work space for the imaginary parts of the coefficients;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on absdx, used as stop criterium;
  --   rc       work space for real parts of the columns of A(0);
  --   ic       work space for imaginary parts of the columns of A(0);
  --   rv       work space for all real parts of all A(k), for k in 1..degree;
  --   iv       work space for all imag parts of all A(k), for k in 1..degree;
  --   rb       real parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ib       imaginary parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ry       allocated work space vector of range 1..s.dim;
  --   iy       allocated work space vector of range 1..s.dim;
  --   scale    if true, then the k-th component of the update dx
  --            is divided by k!, otherwise no scaling to dx is applied;
  --   verbose  if verbose, then the progress index will be written;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the update to the last coefficient;
  --   fail     true if absdx > tol after nbrit iterations;
  --   rcond    estimate for the inverse condition number of A(0);
  --   ipvt     pivoting of the LU factorization on the lead matrix;
  --   rc       real parts of the output of lufac on A(0);
  --   ic       imaginary parts of the output of lufac on A(0);
  --   rb       rb(k) stores the real parts of the solution x(k);
  --   ib       ib(k) stores the imaginary parts of the solution x(k).

-- NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf, in double precision,
  --   using LU factorization to solve the linear systems.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   rx       work space for the real parts of the coefficients;
  --   ix       work space for the imaginary parts of the coefficients;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on absdx, used as stop criterium;
  --   wrk      work space for the matrix series solver;
  --   scale    if true, then scaling is applied after the Newton steps,
  --            otherwise no scaling is applied;
  --   verbose  if verbose, then the progress index will be written;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the update to the last coefficient;
  --   fail     true if absdx > tol after nbrit iterations;
  --   info     info from the LU factorization;
  --   ipvt     pivoting of the LU factorization on the lead matrix.

-- NEWTON STEPS WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf, in double precision,
  --   using LU factorization to solve the linear systems,
  --   with an estimate for the condition number.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   s        system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   rx       work space for the real parts of the coefficients;
  --   ix       work space for the imaginary parts of the coefficients;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on absdx, used as stop criterium;
  --   wrk      work space for the matrix series solver;
  --   scale    if true, then scaling is applied after the Newton steps,
  --            otherwise no scaling is applied;
  --   verbose  if verbose, then the progress index will be written;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the update to the last coefficient;
  --   fail     true if absdx > tol after nbrit iterations;
  --   rcond    estimate for the inverse condition number;
  --   ipvt     pivoting of the LU factorization on the lead matrix.

-- NEWTON STEPS WITH QR ON COEFFICIENT CONVOLUTIONS :

  procedure QR_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure QR_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf, in double precision,
  --   solving the linear systems in the least squares sense with QR.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   dx       work space for the update to scf, delinearized;
  --   xd       work space for the update to scf, in linearized format;
  --   rx       work space for the real parts of the coefficients;
  --   ix       work space for the imaginary parts of the coefficients;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on absdx, used as stop criterium;
  --   qraux    information to recover the orthogonal part;
  --   w1       work space vector of range 1..n, n = number of rows;
  --   w2       work space vector of range 1..n, n = number of rows;
  --   w3       work space vector of range 1..n, n = number of rows;
  --   w4       work space vector of range 1..n, n = number of rows;
  --   w5       work space vector of range 1..n, n = number of rows.
  --   wrk      work space for the matrix series solver;
  --   scale    if true, then scaling is applied after the Newton steps,
  --            otherwise no scaling is applied;
  --   verbose  if verbose, then the progress index will be written;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   dx       delinearized update to the coefficients;
  --   xd       update to the coefficients, in linearized form;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the update to the last coefficient;
  --   fail     true if absdx > tol after nbrit iterations;
  --   ipvt     pivoting of the LU factorization on the lead matrix.

-- NEWTON STEPS WITH SVD ON COEFFICIENT CONVOLUTIONS :

  procedure SVD_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                svl : out Standard_Complex_Vectors.Vector;
                U,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                svl : out Standard_Complex_Vectors.Vector;
                U,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf, in double precision,
  --   solving the linear systems with singular value decompositions.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   dx       work space for the update to scf, delinearized;
  --   xd       work space for the update to scf, in linearized format;
  --   rx       work space for the real parts of the coefficients;
  --   ix       work space for the imaginary parts of the coefficients;
  --   maxit    maximum number of iterations;
  --   tol      tolerance on absdx, used as stop criterium;
  --   ewrk     work space allocated for the SVD of the lead A(0);
  --   wrkv     work space vector for the next coefficient computation;
  --   scale    if true, then scaling is applied after the Newton steps,
  --            otherwise no scaling is applied;
  --   verbose  if verbose, then the progress index will be written;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   dx       delinearized update to the coefficients;
  --   xd       update to the coefficients, in linearized form;
  --   nbrit    number of iterations done;
  --   absdx    absolute value of the update to the last coefficient;
  --   fail     true if absdx > tol after nbrit iterations;
  --   svl      vector of range 1..mm, where mm = min(n+1,p),
  --            where n = A(0)'last(1) and p = A(0)'last(2),
  --            the first min(n,p) entries of s contain the singular values
  --            of x arranged in descending order of magnitude;
  --   U        matrix with n rows and k columns, 
  --            if joba = 1, then k = n, if joba >= 2 then k = min(n,p),
  --            u contains the matrix of left singular vectors,
  --            u is not referenced if joba = 0, if n <= p or if joba > 2,
  --            then u may be identified with x in the subroutine call;
  --   V        matrix with p rows and p columns,
  --            v contains the matrix of right singular vectors,
  --            v is not referenced if jobb = 0, if p <= n, then v may be
  --            identified with x in the subroutine call;
  --   info     the singular values (and their corresponding singular vectors)
  --            s(info+1),s(info+2),...,s(m) are correct (here m=min(n,p)),
  --            thus if info = 0, all the singular values and their vectors
  --            are correct, in any event, the matrix b = ctrans(u)*x*v is
  --            the bidiagonal matrix with the elements of s on its diagonal
  --            and the elements of e on its super diagonal (ctrans(u) is the
  --            conjugate-transpose of u), thus the singular values of x 
  --            and b are the same;
  --   rcond    estimate for the inverse of the condition number,
  --            if close to zero, then the Jacobian matrix at scf is
  --            ill conditioned and scf may be wrongly updated.

end Staggered_Newton_Convolutions;
