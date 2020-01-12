with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Newton_Convolutions is

-- DESCRIPTION :
--   Newton's method on power series computed with convolution circuits.
--   The computational procedures in this package take their work space
--   from the input arguments and do not allocate for thread safety.

-- AUXILIARIES :

  function Series_Coefficients
             ( v : Standard_Complex_Vectors.Vector;
               d : integer32 )
             return Standard_Complex_VecVecs.VecVec;
  function Series_Coefficients
             ( v : DoblDobl_Complex_Vectors.Vector;
               d : integer32 )
             return DoblDobl_Complex_VecVecs.VecVec;
  function Series_Coefficients
             ( v : QuadDobl_Complex_Vectors.Vector;
               d : integer32 )
             return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the coefficients of a vector of series,
  --   with leading coefficients the components of a solution v.
  --   The series are of degree d.

  procedure Minus ( v : in Standard_Complex_VecVecs.VecVec );
  procedure Minus ( v : in DoblDobl_Complex_VecVecs.VecVec );
  procedure Minus ( v : in QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Flips the sign of all coefficients in v.

  procedure Update ( x,y : in Standard_Complex_VecVecs.VecVec );
  procedure Update ( x,y : in DoblDobl_Complex_VecVecs.VecVec );
  procedure Update ( x,y : in QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Adds to all coefficients of x the corresponding coefficient of y.

  -- REQUIRED : x'range = y'range, and for all k in x'range
  --   x(k)'range = y(k)'range.

  function Max ( v : Standard_Complex_Vectors.Link_to_Vector )
               return double_float;
  function Max ( v : DoblDobl_Complex_Vectors.Link_to_Vector )
               return double_double;
  function Max ( v : QuadDobl_Complex_Vectors.Link_to_Vector )
               return quad_double;

  -- DESCRIPTION :
  --   Returns the largest absolute value over all values in v.

  function Max ( v : Standard_Complex_VecVecs.VecVec ) return double_float;
  function Max ( v : DoblDobl_Complex_VecVecs.VecVec ) return double_double;
  function Max ( v : QuadDobl_Complex_VecVecs.VecVec ) return quad_double;

  -- DESCRIPTION :
  --   Returns the largest absolute value over all values in v.

-- ONE NEWTON STEP WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
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
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   absdx    the absolute value of the maximal component of the update dx;
  --   info     info from the LU factorization;
  --   ipvt     pivoting of the LU factorization on the lead matrix.

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                absdx,rcond : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                absdx,rcond : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                absdx,rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure LU_Newton_Step
              ( file : in file_type;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                absdx,rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
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
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure QR_Newton_Step
              ( file : in file_type;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure QR_Newton_Step
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
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
                vrblvl : in integer32 := 0 );
  procedure QR_Newton_Step
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double;
                qraux : out QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure QR_Newton_Step
              ( file : in file_type;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double;
                qraux : out QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies one Newton step on the convolution circuits c,
  --   departing from the series coefficients in s,
  --   in double, double double, or quad double precision,
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
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   dx       delinearized update to the coefficients;
  --   xd       update to the coefficients, in linearized form;
  --   absdx    the absolute value of the maximal component of the update dx;
  --   ipvt     pivoting of the LU factorization on the lead matrix.

-- ONE NEWTON STEP WITH SVD :

  procedure SVD_Newton_Step
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float;
                svl : out Standard_Complex_Vectors.Vector;
                U,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Step
              ( file : in file_type;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float;
                svl : out Standard_Complex_Vectors.Vector;
                U,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Step
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double;
                svl : out DoblDobl_Complex_Vectors.Vector;
                U,V : out DoblDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_double;
                ewrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in DoblDobl_Complex_Vectors.Link_to_Vector;
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
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Step
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double;
                svl : out QuadDobl_Complex_Vectors.Vector;
                U,V : out QuadDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out quad_double;
                ewrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in QuadDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );
  procedure SVD_Newton_Step
              ( file : in file_type;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double;
                svl : out QuadDobl_Complex_Vectors.Vector;
                U,V : out QuadDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out quad_double;
                ewrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in QuadDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies one Newton step on the convolution circuits c,
  --   departing from the series coefficients in s,
  --   in double, double double, or quad double precision,
  --   solving the linear system with SVD.

  -- ON ENTRY :
  --   file     if provided, the intermediate coefficient vectors
  --            are written to file, otherwise the procedure is silent;
  --   s        system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   dx       work space for the update to scf, delinearized;
  --   xd       work space for the update to scf, in linearized format;
  --   ewrk     work space allocated for the SVD of the lead A(0);
  --   wrkv     work space vector for the next coefficient computation.
  --   vrblvl   if positive, the name of the procedure is written to screen.

  -- ON RETURN :
  --   scf      updated coefficients of the series solution;
  --   absdx    the absolute value of the maximal component of the update dx;
  --   dx       delinearized update to the coefficients;
  --   xd       update to the coefficients, in linearized form;
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

end Newton_Convolutions;
