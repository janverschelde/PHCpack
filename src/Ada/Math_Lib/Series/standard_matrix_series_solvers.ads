with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Dense_Vector_Series;
with Standard_Dense_Matrix_Series;

package Standard_Matrix_Series_Solvers is

-- DESCRIPTION :
--   Given a linear system of truncated power series, linearization with
--   matrix series and vector series solves the linear system.

  procedure Solve_Lead_by_lufac
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                a0lu : out Standard_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out Standard_Dense_Vector_Series.Vector );

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
  --            and x.deg = 0, provided info = 0,
  --            if info /= 0, then x.deg = -1 and x.cff(0) is undefined.

  procedure Solve_Lead_by_lufco
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                a0lu : out Standard_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float;
                x : out Standard_Dense_Vector_Series.Vector );

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
  --            and x.deg = 0, provided 1.0 + rcond /= 1.0,
  --            if 1.0 + rcond = 1.0, then x.deg = -1
  --            and x.cff(0) is undefined.

  procedure Solve_Lead_by_QRLS
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                a0qr : out Standard_Complex_Matrices.Matrix;
                qraux : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out Standard_Dense_Vector_Series.Vector );

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

  procedure Solve_Next_by_lusolve
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                a0lu : in Standard_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                x : in out Standard_Dense_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies LU back substitution, calling lusolve,
  --   on the linear system with a modified right hand side vector,
  --   using previously computed coefficient series vector of x.

  -- REQUIRED :
  --   All coefficients up to x.deg are defined and A*x = b is square.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   a0lu     output of Solve_Lead_by_lufac, contains the LU
  --            factorization of the leading coefficient of A;
  --   ipvt     output of Solve_Lead_by_lufac, contains the pivoting
  --            information in the LU factorization of A.cff(0);
  --   x        previously computed coefficients of the solution,
  --            at the very least, x.cff(0) must be defined.

  -- ON RETURN :
  --   x        computed coefficient at x.deg+1 with respect to input.

  procedure Solve_Next_by_QRLS
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                a0qr : in Standard_Complex_Matrices.Matrix;
                qraux : in Standard_Complex_Vectors.Vector;
                ipvt : in Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : in out Standard_Dense_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies QR decomposition and least squares solving to compute the
  --   next coefficient of the solution series A*x = b.

  -- REQUIRED :
  --   All coefficients up to x.deg are defined.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.
  --   a0qr     the QR decomposition of the leading coefficient of A,
  --            as output of Solve_Lead_by_QRLS;
  --   qraux    information to recover the orthogonal part,
  --            as output of Solve_Lead_by_QRLS;
  --   ipvt     pivoting information if that was requested,
  --            as output of Solve_Lead_by_QRLS;

  -- ON RETURN :
  --   info     is zero of nonsingular, otherwise, a nonzero info
  --            indicates a singular matrix;
  --   x        computed coefficient at x.deg+1 with respect to input.

  procedure Solve_by_lufac
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                info : out integer32;
                x : out Standard_Dense_Vector_Series.Vector );

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
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                rcond : out double_float;
                x : out Standard_Dense_Vector_Series.Vector );

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
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                info : out integer32;
                x : out Standard_Dense_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using QR decomposition of the
  --   leading coefficient matrix of A and least squares solving.

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

end Standard_Matrix_Series_Solvers;
