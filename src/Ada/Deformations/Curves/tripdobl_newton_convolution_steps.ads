with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Standard_Integer_Vectors;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_Matrices;
with TripDobl_Speelpenning_Convolutions;

package TripDobl_Newton_Convolution_Steps is

-- DESCRIPTION :
--   Several steps with Newton's method on power series are performed
--   with convolutions on the coefficients and linearization,
--   in triple double precision arithmetic.

-- NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( csr : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out triple_double;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type; 
                csr : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out triple_double;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf,
  --   using LU factorization to solve the linear systems.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
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
              ( csr : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out triple_double;
                fail : out boolean; rcond : out triple_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector;
		scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Steps
              ( file : in file_type; 
                csr : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out triple_double;
                fail : out boolean; rcond : out triple_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector;
		scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf,
  --   using LU factorization to solve the linear systems,
  --   with an estimate for the condition number.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   s        system of convolution circuits;
  --   scf      vector of coefficients of power series;
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

-- NEWTON STEPS WITH QR :

  procedure QR_Newton_Steps
              ( csr : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out triple_double;
                fail : out boolean;
                qraux : out TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector;
		scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure QR_Newton_Steps
              ( file : in file_type;
                csr : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out triple_double;
                fail : out boolean;
                qraux : out TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector;
		scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf,
  --   solving the linear systems in the least squares sense with QR.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   dx       work space for the update to scf, delinearized;
  --   xd       work space for the update to scf, in linearized format;
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

-- NEWTON STEPS WITH SVD :

  procedure SVD_Newton_Steps
              ( csr : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out triple_double;
                fail : out boolean;
                svl : out TripDobl_Complex_Vectors.Vector;
                U,V : out TripDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out triple_double;
                ewrk : in TripDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in TripDobl_Complex_Vectors.Link_to_Vector;
		scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Steps
              ( file : in file_type;
                csr : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out triple_double;
                fail : out boolean;
                svl : out TripDobl_Complex_Vectors.Vector;
                U,V : out TripDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out triple_double;
                ewrk : in TripDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in TripDobl_Complex_Vectors.Link_to_Vector;
		scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies several Newton steps on the convolution circuits csr,
  --   departing from the series coefficients in scf,
  --   solving the linear systems with singular value decompositions.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   csr      system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   dx       work space for the update to scf, delinearized;
  --   xd       work space for the update to scf, in linearized format;
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

end TripDobl_Newton_Convolution_Steps;
