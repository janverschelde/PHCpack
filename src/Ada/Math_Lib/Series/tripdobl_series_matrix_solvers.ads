with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Triple_Double_Numbers;               use Triple_Double_Numbers;
with TripDobl_Complex_Numbers;            use TripDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_VecMats;
with TripDobl_Complex_Vector_Series;
with TripDobl_Complex_Matrix_Series;

package TripDobl_Series_Matrix_Solvers is

-- DESCRIPTION :
--   Given a linear system of truncated power series, linearization with
--   matrix series and vector series solves the linear system.

  procedure Solve_Lead_by_lufac
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                a0lu : out TripDobl_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies LU factorization and back substitution to compute the
  --   constant coefficient of the solution series to A*x = b.

  -- REQUIRED :
  --   A.deg >= 0 and b.deg >= 0.  Moreover, the system is square.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   a0lu     LU factorization of A.cff(0);
  --   ipvt     pivoting information on the LU factorization of A.cff(0);
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac;
  --   x        x.cff(0) is the constant coefficient of the solution
  --            provided info = 0,
  --            if info /= 0, then x.cff(0) is undefined.

  procedure Solve_Lead_by_lufco
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                a0lu : out TripDobl_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out triple_double;
                x : out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies LU factorization and back substitution to compute the
  --   constant coefficient of the solution series to A*x = b.
  --   An estimate for the inverse of the condition number is returned.

  -- REQUIRED :
  --   A.deg >= 0 and b.deg >= 0.  Moreover, the system is square.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   a0lu     LU factorization of A.cff(0);
  --   ipvt     pivoting information on the LU factorization of A.cff(0);
  --   rcond    computed by lufco, if 1.0 + rcond = 1.0, then the lead
  --            coefficient matrix should be considered as singular;
  --   x        x.cff(0) is the constant coefficient of the solution
  --            provided 1.0 + rcond /= 1.0,
  --            if 1.0 + rcond = 1.0, then x.cff(0) is undefined.

  procedure Solve_Lead_by_QRLS
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                a0qr : out TripDobl_Complex_Matrices.Matrix;
                qraux : out TripDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies QR decomposition and least squares solving to compute the
  --   constant coefficient of the solution series A*x = b.

  -- REQUIRED :
  --   A.deg >= 0 and b.deg >= 0.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   a0qr     the QR decomposition of the leading coefficient of A;
  --   qraux    information to recover the orthogonal part;
  --   ipvt     pivoting information if that was requested;
  --   info     is zero of nonsingular, otherwise, a nonzero info
  --            indicates a singular matrix.
  --   x        x.cff(0) is the constant coefficient of the solution,
  --            if info /= 0, then x.cff(0) may be inaccurate.

  procedure Solve_Lead_by_SVD
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                S : out TripDobl_Complex_Vectors.Vector;
                U,V : out TripDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out triple_double;
                x : out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies Singular Value Decomposition (SVD) to compute the
  --   constant coefficient of the solution series A*x = b.

  -- REQUIRED :
  --   A.deg >= 0 and b.deg >= 0.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   S        vector of range 1..mm, where mm = min(n+1,p),
  --            where n = A.cff(0)'last(1) and p = A.cff(0)'last(2),
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
  --   rcond    is the outcome of Inverse_Condition_Number(S),
  --            the inverse condition number based on the singular values;
  --   x        x.cff(0) is the constant coefficient of the solution,
  --            if info /= 0, then x.cff(0) may be inaccurate.

  procedure Solve_Next_by_lusolve
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                a0lu : in TripDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                idx : in integer32;
                x : in out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies LU back substitution, calling lusolve,
  --   on the linear system with a modified right hand side vector,
  --   using previously computed coefficient series vector of x.

  -- REQUIRED :
  --   All coefficients up to idx-1 are defined and A*x = b is square.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   a0lu     output of Solve_Lead_by_lufac, contains the LU
  --            factorization of the leading coefficient of A;
  --   ipvt     output of Solve_Lead_by_lufac, contains the pivoting
  --            information in the LU factorization of A.cff(0);
  --   idx      current coefficient index of x,
  --            x.cff(k) for k in range 0..idx-1 have been computed,
  --            x.cff(idx) will be computed;
  --   x        previously computed coefficients of the solution,
  --            at the very least, x.cff(0) must be defined.

  -- ON RETURN :
  --   x        computed coefficient at idx with respect to input.

  procedure Solve_Next_by_QRLS
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                a0qr : in TripDobl_Complex_Matrices.Matrix;
                qraux : in TripDobl_Complex_Vectors.Vector;
                idx : in integer32; info : out integer32;
                x : in out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies QR decomposition and least squares solving to compute the
  --   next coefficient of the solution series A*x = b.

  -- REQUIRED :
  --   All coefficients in x up to idx-1 are defined.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   a0qr     the QR decomposition of the leading coefficient of A,
  --            as output of Solve_Lead_by_QRLS;
  --   qraux    information to recover the orthogonal part,
  --            as output of Solve_Lead_by_QRLS;
  --   idx      current coefficient index of x,
  --            x.cff(k) for k in range 0..idx-1 have been computed,
  --            x.cff(idx) will be computed;
  --   x        previously computed coefficients of the solution,
  --            at the very least, x.cff(0) must be defined.

  -- ON RETURN :
  --   info     is zero of nonsingular, otherwise, a nonzero info
  --            indicates a singular matrix;
  --   x        computed coefficient at idx with respect to input.

  procedure Solve_Next_by_SVD
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                S : in TripDobl_Complex_Vectors.Vector;
                U,V : in TripDobl_Complex_Matrices.Matrix;
                idx : in integer32;
                x : in out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies the SVD decomposition of the lead coefficient of A
  --   to compute the next coefficient of the solution series A*x = b.

  -- REQUIRED :
  --   All coefficients up to x.deg are defined.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   S        singular values of the lead coefficient of A,
  --            computed by Solve_Lead_by_SVD;
  --   U,V      computed by Solve_Lead_by_SVD;
  --   idx      current coefficient index of x,
  --            x.cff(k) for k in range 0..idx-1 have been computed,
  --            x.cff(idx) will be computed;
  --   x        previously computed coefficients of the solution,
  --            at the very least, x.cff(0) must be defined.

  -- ON RETURN :
  --   x        computed coefficient at idx with respect to input.

  procedure Solve_by_lufac
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                info : out integer32;
                x : out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using LU factorization on the
  --   leading coefficient matrix of A, without condition number estimate.

  -- REQUIRED :
  --   A.deg >= 0 and b.deg >= 0.  Moreover, the system is square.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   info     if info /= 0, then the lead coefficient matrix of A
  --            was deemed singular and x is undefined,
  --            if info = 0, then the system is regular;
  --   x        all coefficients of the solution series up to b.deg,
  --            provided info = 0.

  procedure Solve_by_lufco
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                rcond : out triple_double;
                x : out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using LU factorization on the
  --   leading coefficient matrix of A, with condition number estimate.

  -- REQUIRED :
  --   A.deg >= 0 and b.deg >= 0.  Moreover, the system is square.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   rcond    estimate for the inverse condition number,
  --            if 1.0 + rcond = 1.0, then the lead coefficient matrix of A
  --            was deemed singular and x is undefined,
  --            if 1.0 + rcond /= 1.0, then the system is regular;
  --   x        all coefficients of the solution series up to b.deg,
  --            provided 1.0 + rcond /= 1.0.

  procedure Solve_by_QRLS
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                info : out integer32;
                x : out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using the QR decomposition of the
  --   leading coefficient matrix of A for least squares solving.

  -- REQUIRED :
  --   A.deg >= 0 and b.deg >= 0.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   info     if info /= 0, then the lead coefficient matrix of A
  --            was deemed singular and x is undefined,
  --            if info = 0, then the system is regular;
  --   x        all coefficients of the solution series up to b.deg,
  --            provided info = 0.

  procedure Solve_by_SVD
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                info : out integer32; rcond : out triple_double;
                x : out TripDobl_Complex_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using the SVD of the
  --   leading coefficient matrix of A for least squares solving.

  -- REQUIRED :
  --   A.deg >= 0 and b.deg >= 0.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   info     see the output of Solve_Lead_by_SVD;
  --   rcond    inverse condition number computed from the singular
  --            values of the lead coefficient of A;
  --   x        all coefficients of the solution series up to b.deg,
  --            provided rcond /= 0.0.

  procedure Echelon_Solve
              ( A : in TripDobl_Complex_Matrix_Series.Matrix;
                b : in TripDobl_Complex_Vector_Series.Vector;
                det : out Complex_Number;
                xp : out TripDobl_Complex_Vector_Series.Vector;
                xn : out TripDobl_Complex_Vector_Series.Vector );
            
  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using the lower triangular
  --   echelon form of the Hermite-Laurent system defined by A and b.
  --   The solution may have negative exponents, maybe a Laurent series.

  -- REQUIRED : A.deg = b.deg >= 0.

  -- REQUIRED :
  --   A.deg >= 0 and b.deg >= 0.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   det      determinant of the lower triangular echelon form,
  --            if zero, then the solution may be a formal Laurent series
  --            as nonzero coefficients with negative exponents appear,
  --            if nonzero, then xn.deg = -1;
  --   xp       all coefficients of the solution series up to b.deg,
  --            provided det /= 0;
  --   xn       if xn.deg = -1, then there are no terms in the solution
  --            series with negative exponents,
  --            if xn.deg > 0, the xn.cff(k) stores the coefficient with
  --            t**(-k) in the Laurent series of the solution.

-- ON FLATTENED DATA STRUCTURES :

  procedure Solve_Lead_by_lufac
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );

  -- DESCRIPTION :
  --   Applies LU factorization and back substitution to compute the
  --   constant coefficient of the solution series to A*x = b.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A are square.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series.

  -- ON RETURN :
  --   A        A(0) contains the output of lufac on A(0);
  --   b        b(0) is the constant coefficient of the solution,
  --            provided info = 0, otherwise b(0) is unchanged;
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac.

  procedure Solve_Lead_by_lufco
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out triple_double );

  -- DESCRIPTION :
  --   Applies LU factorization and back substitution to compute the
  --   constant coefficient of the solution series to A*x = b.
  --   An estimate for the inverse of the condition number is returned.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A are square.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series.

  -- ON RETURN :
  --   A        A(0) contains the output of lufco on A(0);
  --   b        b(0) is the constant coefficient of the solution,
  --            provided 1.0 + rcond /= 1.0, otherwise b(0) is unchanged;
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   rcond    computed by lufco, if 1.0 + rcond = 1.0, then the lead
  --            coefficient matrix should be considered as singular.

  procedure Solve_Lead_by_QRLS
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x0 : in TripDobl_Complex_Vectors.Link_to_Vector;
                qraux : out TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );

  -- DESCRIPTION :
  --   Applies QR decomposition and least squares solving to compute the
  --   constant coefficient of the solution series A*x = b.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A have the same dimensions.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series;
  --   x0       space allocated for the leading coefficient of the solution;
  --   w1       work space vector of range 1..n, n = number of rows;
  --   w2       work space vector of range 1..n, n = number of rows;
  --   w3       work space vector of range 1..n, n = number of rows;
  --   w4       work space vector of range 1..n, n = number of rows;
  --   w5       work space vector of range 1..n, n = number of rows.

  -- ON RETURN :
  --   A        A(0) contains QR decomposition of A(0);
  --   x0       the constant coefficient of the solution,
  --            if info /= 0, then b(0) may be inaccurate.
  --   qraux    information to recover the orthogonal part;
  --   ipvt     pivoting information if that was requested;
  --   info     is zero of nonsingular, otherwise, a nonzero info
  --            indicates a singular matrix.

  procedure Solve_Lead_by_SVD
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x0 : in TripDobl_Complex_Vectors.Link_to_Vector;
                S : out TripDobl_Complex_Vectors.Vector;
                U,V : out TripDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out triple_double;
                ewrk : in TripDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in TripDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies the SVD decomposition of the lead coefficient of A
  --   to compute the next coefficient of the solution series A*x = b.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A have the same dimensions.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series;
  --   x0       space allocated for the leading coefficient of the solution;
  --   ewrk     allocated work space of range 1..p, p = A(0)'last(2);
  --   wrkv     allocated work space of range 1..n, n = A(0)'last(1).

  -- ON RETURN :
  --   A        A(0) used as work space to compute the SVD;
  --   S        vector of range 1..mm, where mm = min(n+1,p),
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
  --   rcond    is the outcome of Inverse_Condition_Number(S),
  --            the inverse condition number based on the singular values;
  --   x        the constant coefficient of the solution.

  procedure Matrix_Vector_Multiply
              ( A : in TripDobl_Complex_Matrices.Link_to_Matrix;
                x,y : in TripDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION : computes y = A*x.

  -- REQUIRED : y'range = A'range(1) and x'range = A'range(2).

  procedure Subtract ( x,y : in TripDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION : x := x - y.

  -- REQUIRED : x'range = y'range.

  procedure Solve_Next_by_lusolve
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                idx : in integer32;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies LU back substitution, calling lusolve,
  --   on the linear system with a modified right hand side vector,
  --   using previously computed coefficient series vector of x.

  -- REQUIRED :
  --   All coefficients up to idx-1 are defined and A*x = b is square.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series,
  --            A(0) contains the output of lufac or lufco;
  --   b        the right hand side coefficients of a vector series,
  --            b(k) for k in 0..idx-1 contains the solutions;
  --   ipvt     output of Solve_Lead_by_lufac, contains the pivoting
  --            information in the LU factorization of A(0);
  --   idx      current coefficient index of the solution,
  --            b(k) for k in range 0..idx-1 have been computed,
  --            b(idx) will be computed;
  --   wrk      work vector, allocated of range at least A(0)'range(1).

  -- ON RETURN :
  --   b        computed coefficient at idx with respect to input.

  procedure Solve_Next_by_QRLS
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                qraux : in TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_Vectors.Vector;
                idx : in integer32; info : out integer32;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies QR decomposition and least squares solving to compute the
  --   next coefficient of the solution series A*x = b.

  -- REQUIRED :
  --   All coefficients in x up to idx-1 are defined.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series,
  --            A(0) contains the QR decomposition of the lead matrix,
  --            obtained as output of Solve_Lead_by_QRLS;
  --   b        the right hand side coefficients of a vector series;
  --   x        x(k) for k in 0..idx-1 contains the solutions;
  --   qraux    information to recover the orthogonal part,
  --            as output of Solve_Lead_by_QRLS;
  --   w1       work space vector of range 1..n, n = number of rows;
  --   w2       work space vector of range 1..n, n = number of rows;
  --   w3       work space vector of range 1..n, n = number of rows;
  --   w4       work space vector of range 1..n, n = number of rows;
  --   w5       work space vector of range 1..n, n = number of rows.
  --   idx      current coefficient index of x,
  --            b(k) for k in range 0..idx-1 have been computed,
  --            b(idx) will be computed;
  --   wrk      allocated work space for the matrix-by-vector computations.

  -- ON RETURN :
  --   info     is zero of nonsingular, otherwise, a nonzero info
  --            indicates a singular matrix;
  --   b        updated right hand side after back substitution;
  --   x        computed coefficient at idx with respect to input.

  procedure Solve_Next_by_SVD
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                S : in TripDobl_Complex_Vectors.Vector;
                U,V : in TripDobl_Complex_Matrices.Matrix;
                idx : in integer32;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies the SVD decomposition of the lead coefficient of A
  --   to compute the next coefficient of the solution series A*x = b.

  -- REQUIRED :
  --   All coefficients in x up to idx-1 are defined.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series,
  --   b        the right hand side as a vector series;
  --   S        singular values of the lead coefficient of A,
  --            computed by Solve_Lead_by_SVD;
  --   U,V      computed by Solve_Lead_by_SVD;
  --   idx      current coefficient index of x,
  --            x(k) for k in range 0..idx-1 have been computed,
  --            x(idx) will be computed;
  --   x        previously computed coefficients of the solution,
  --            at the very least, x(0) must be defined;
  --   wrk      work space allocated for the right hand side computations.

  -- ON RETURN :
  --   b        updated right hand side after back substitution;
  --   x        computed coefficient at idx with respect to input.

  procedure Solve_by_lufac
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using LU factorization on the
  --   leading coefficient matrix of A, without condition number estimate.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A are square.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series;
  --   wrk      work vector, allocated of range at least A(0)'range(1).

  -- ON RETURN :
  --   A        A(0) contains the output of lufac on A(0);
  --   b        b contains the coefficients of the solution series x,
  --            provided info = 0, otherwise b is unchanged;
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac.

  procedure Solve_by_lufco
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out triple_double;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using LU factorization on the
  --   leading coefficient matrix of A, with condition number estimate.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A are square.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series;
  --   wrk      work vector, allocated of range at least A(0)'range(1).

  -- ON RETURN :
  --   A        A(0) contains the output of lufac on A(0);
  --   b        b contains the coefficients of the solution series x,
  --            provided info = 0, otherwise b is unchanged;
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   rcond    computed by lufco, if 1.0 + rcond = 1.0, then the lead
  --            coefficient matrix should be considered as singular.

  procedure Solve_by_QRLS
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                qraux : out TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in TripDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using the QR decomposition of the
  --   leading coefficient matrix of A for least squares solving.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A have the same dimension.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series;
  --   x        space allocated for the solution series;
  --   w1       work space vector of range 1..n, n = number of rows;
  --   w2       work space vector of range 1..n, n = number of rows;
  --   w3       work space vector of range 1..n, n = number of rows;
  --   w4       work space vector of range 1..n, n = number of rows;
  --   w5       work space vector of range 1..n, n = number of rows.
  --   wrk      work vector, allocated of range at least A(0)'range(1).

  -- ON RETURN :
  --   A        A(0) contains the output of QRD on A(0);
  --   b        modified right hand side vectors after back substitution;
  --   x        contains the coefficients of the solution series x,
  --            provided info = 0,
  --   qraux    information to recover the orthogonal part;
  --   ipvt     pivoting information if that was requested;
  --   info     is zero of nonsingular, otherwise, a nonzero info
  --            indicates a singular matrix.

  procedure Solve_by_SVD
              ( A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                S : out TripDobl_Complex_Vectors.Vector;
                U,V : out TripDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out triple_double;
                ewrk : in TripDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in TripDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using the SVD of the
  --   leading coefficient matrix of A for least squares solving.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A have the same dimension.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series;
  --   x        space allocated for the solution series;
  --   ewrk     work space allocated for the SVD of the lead A(0);
  --   wrkv     work space vector for the next coefficient computation.

  -- ON RETURN :
  --   A        A(0) modified as work space in SVD computation;
  --   b        modified right hand side vectors after back substitution;
  --   S,U,V    see the output of Solve_Lead_by_SVD;
  --   info     see the output of Solve_Lead_by_SVD;
  --   rcond    inverse condition number computed from the singular
  --            values of the lead coefficient of A;
  --   x        all coefficients of the solution series up to b'last,
  --            provided rcond /= 0.0.

end TripDobl_Series_Matrix_Solvers;
