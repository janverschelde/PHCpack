with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Speelpenning_Convolutions;

package DoblDobl_Newton_Convolutions is

-- DESCRIPTION :
--   Newton's method on power series computed with convolution circuits.
--   The computational procedures in this package take their work space
--   from the input arguments and do not allocate for thread safety.
--   All computations are performed in double double precision.

-- AUXILIARIES :

  function Series_Coefficients
             ( v : DoblDobl_Complex_Vectors.Vector;
               d : integer32 )
             return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the coefficients of a vector of series,
  --   with leading coefficients the components of a solution v.
  --   The series are of degree d.

  procedure Minus ( v : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Flips the sign of all coefficients in v.

  procedure Minus ( deg : in integer32;
                    v : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Flips the sign of all coefficients in v,
  --   up to the given degree deg.

  procedure Power_Divide
	      ( x : in DoblDobl_Complex_VecVecs.VecVec;
                f : in double_double );

  -- DESCRIPTION :
  --   Divides all numbers in x(k) by f**k, for k > 0.

  -- REQUIRED :
  --   x contains the linearized representation of a vector of power series,
  --   that is: x(k) contains all coefficients with t^k, for all components
  --   of the vector of power series.

  procedure Update ( x,y : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Adds to all coefficients of x the corresponding coefficient of y.

  procedure Update ( deg : in integer32;
                     x,y : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Adds to all coefficients of x the corresponding coefficient of y,
  --   up to the given degree deg.

  -- REQUIRED : x'range = y'range, and for all k in x'range
  --   x(k)'range = y(k)'range.

  function Max ( v : DoblDobl_Complex_Vectors.Link_to_Vector )
               return double_double;

  -- DESCRIPTION :
  --   Returns the largest absolute value over all values in v.

  function Max ( deg : integer32;
                 v : DoblDobl_Complex_Vectors.Link_to_Vector )
               return double_double;

  -- DESCRIPTION :
  --   Returns the largest absolute value over all values in v,
  --   up to the last index, as defined by the value of deg.

  function Max ( v : DoblDobl_Complex_VecVecs.VecVec ) return double_double;

  -- DESCRIPTION :
  --   Returns the largest absolute value over all values in v.

  function Max ( deg : integer32;
                 v : DoblDobl_Complex_VecVecs.VecVec ) return double_double;

  -- DESCRIPTION :
  --   Returns the largest absolute value over all values in v(k),
  --   for the numbers in v(k) ranging up to the given degree deg.

  procedure MaxIdx ( v : in DoblDobl_Complex_VecVecs.VecVec;
                     tol : in double_float;
                     maxval : out double_double; idx : out integer32 );

  -- DESCRIPTION :
  --   Returns in idx the highest index in v for which 
  --   maxval = Max(v(idx)) <= tol.
  --   If idx < v'first, then already Max(v(v'first)) > tol,
  --   otherwise for all k from v'first to idx, Max(v(k)) <= tol.

  -- REQUIRED :
  --   The vector v stores the updates to the coefficients of a power series
  --   in linearized form, that is: v(k) stores all coefficients with t^k,
  --   for k in v'range.

  procedure MaxIdx ( deg : in integer32;
                     v : in DoblDobl_Complex_VecVecs.VecVec;
                     tol : in double_float;
                     maxval : out double_double; idx : out integer32 );

  -- DESCRIPTION :
  --   Returns in idx the highest index in v for which 
  --   maxval = Max(v(idx)) <= tol.
  --   If idx < v'first, then already Max(v(v'first)) > tol,
  --   otherwise for all k from v'first to idx, Max(v(k)) <= tol.
  --   Only coefficients with terms of degree <= deg are considered.

  -- REQUIRED :
  --   The vector v stores the updates to the coefficients of a power series
  --   in linearized form, that is: v(k) stores all coefficients with t^k,
  --   for k in v'range.

-- ONE NEWTON STEP WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
		scaledx : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
		scaledx : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies one Newton step on the convolution circuits c,
  --   departing from the series coefficients in s,
  --   in double, double double, or quad double precision,
  --   using LU factorization to solve the linear systems.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   s        system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   wrk      work space for the matrix series solver;
  --   scaledx  if true, then the k-th component of the update dx
  --            is divided by k!, otherwise no scaling to dx is applied;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   absdx    the absolute value of the maximal component of the update dx;
  --   info     info from the LU factorization;
  --   ipvt     pivoting of the LU factorization on the lead matrix.

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                absdx,rcond : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                absdx,rcond : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies one Newton step on the convolution circuits c,
  --   departing from the series coefficients in s,
  --   using LU factorization to solve the linear systems.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   s        system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   wrk      work space for the matrix series solver;
  --   scaledx  if true, then the k-th component of the update dx
  --            is divided by k!, otherwise no scaling to dx is applied;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   absdx    the absolute value of the maximal component of the update dx;
  --   rcond    estimate for the inverse of the condition number,
  --            if close to zero, then the Jacobian matrix at scf is
  --            ill conditioned and scf may be wrongly updated.
  --   ipvt     pivoting of the LU factorization on the lead matrix.

-- ONE NEWTON STEP WITH QR :

  procedure QR_Newton_Step
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure QR_Newton_Step
              ( file : in file_type;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies one Newton step on the convolution circuits c,
  --   departing from the series coefficients in s,
  --   solving the linear system in the least squares sense with QR.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   s        system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   dx       work space for the update to scf, delinearized;
  --   xd       work space for the update to scf, in linearized format;
  --   qraux    information to recover the orthogonal part;
  --   w1       work space vector of range 1..n, n = number of rows;
  --   w2       work space vector of range 1..n, n = number of rows;
  --   w3       work space vector of range 1..n, n = number of rows;
  --   w4       work space vector of range 1..n, n = number of rows;
  --   w5       work space vector of range 1..n, n = number of rows.
  --   wrk      work space for the matrix series solver;
  --   scaledx  if true, then the k-th component of the update dx
  --            is divided by k!, otherwise no scaling to dx is applied;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   dx       delinearized update to the coefficients;
  --   xd       update to the coefficients, in linearized form;
  --   absdx    the absolute value of the maximal component of the update dx;
  --   ipvt     pivoting of the LU factorization on the lead matrix.

-- ONE NEWTON STEP WITH SVD :

  procedure SVD_Newton_Step
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double;
                svl : out DoblDobl_Complex_Vectors.Vector;
                U,V : out DoblDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_double;
                ewrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Step
              ( file : in file_type;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double;
                svl : out DoblDobl_Complex_Vectors.Vector;
                U,V : out DoblDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_double;
                ewrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies one Newton step on the convolution circuits c,
  --   departing from the series coefficients in s,
  --   solving the linear system with SVD.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   s        system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   dx       work space for the update to scf, delinearized;
  --   xd       work space for the update to scf, in linearized format;
  --   ewrk     work space allocated for the SVD of the lead A(0);
  --   wrkv     work space vector for the next coefficient computation;
  --   scaledx  if true, then the k-th component of the update dx
  --            is divided by k!, otherwise no scaling to dx is applied;
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   absdx    the absolute value of the maximal component of the update dx;
  --   dx       delinearized update to the coefficients;
  --   xd       update to the coefficients, in linearized form;
  --   svl      vector of range 1..min(n+1,p), where n = A(0)'last(1)
  --            and p = A(0)'last(2),
  --            the first min(n,p) entries of s contain the singular values
  --            of x arranged in descending order of magnitude;
  --   U        matrix with n rows and n columns;
  --   V        matrix with p rows and p columns;
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

end DoblDobl_Newton_Convolutions;
